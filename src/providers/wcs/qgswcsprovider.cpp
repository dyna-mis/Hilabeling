/***************************************************************************
  qgswcsprovider.cpp  -  QGIS Data provider for
                         OGC Web Coverage Service layers
                             -------------------
    begin                : 2 July, 2012
    copyright            : (C) 2012 by Radim Blazek
    email                : radim dot blazek at gmail.com

    Based on qgswmsprovider.cpp written by Brendan Morley.

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgslogger.h"
#include "qgswcsprovider.h"
#include "qgscoordinatetransform.h"
#include "qgsdatasourceuri.h"
#include "qgsrasteridentifyresult.h"
#include "qgsrasterlayer.h"
#include "qgsrectangle.h"
#include "qgscoordinatereferencesystem.h"
#include "qgsnetworkaccessmanager.h"
#include "qgsnetworkreplyparser.h"
#include "qgsmessagelog.h"
#include "qgsexception.h"

#ifdef HAVE_GUI
#include "qgswcssourceselect.h"
#include "qgssourceselectprovider.h"
#endif

#include <QNetworkRequest>
#include <QNetworkReply>
#include <QNetworkProxy>

#include <QUrl>
#include <QEventLoop>
#include <QFile>

#ifdef QGISDEBUG
#include <QDir>
#endif

#include "gdalwarper.h"
#include "ogr_srs_api.h"
#include "cpl_conv.h"
#include "cpl_string.h"

#define ERR(message) QGS_ERROR_MESSAGE(message,"WCS provider")
#define SRVERR(message) QGS_ERROR_MESSAGE(message,"WCS server")
#define QGS_ERROR(message) QgsError(message,"WCS provider")

static QString WCS_KEY = QStringLiteral( "wcs" );
static QString WCS_DESCRIPTION = QStringLiteral( "OGC Web Coverage Service version 1.0/1.1 data provider" );

static QString DEFAULT_LATLON_CRS = QStringLiteral( "CRS:84" );

// TODO: colortable - use common baseclass with gdal, mapserver does not support http://trac.osgeo.org/mapserver/ticket/1671

QgsWcsProvider::QgsWcsProvider( const QString &uri, const ProviderOptions &options )
  : QgsRasterDataProvider( uri, options )
  , mCachedViewExtent( 0 )
{
  QgsDebugMsg( "constructing with uri '" + mHttpUri + "'." );

  mValid = false;
  mCachedMemFilename = QStringLiteral( "/vsimem/qgis/wcs/%0.dat" ).arg( reinterpret_cast<std::uintptr_t>( this ) );

  if ( !parseUri( uri ) ) return;

  // GetCapabilities and DescribeCoverage
  // TODO(?): do only DescribeCoverage to avoid one request
  // We need to get at least server version, which is not set in of URI (if not part of url)
  // and probably also rangeSet

  QgsDataSourceUri capabilitiesUri;
  capabilitiesUri.setEncodedUri( uri );
  // remove non relevant params
  capabilitiesUri.removeParam( QStringLiteral( "identifier" ) );
  capabilitiesUri.removeParam( QStringLiteral( "crs" ) );
  capabilitiesUri.removeParam( QStringLiteral( "format" ) );
  // TODO: check if successful (add return to capabilities)
  mCapabilities.setUri( capabilitiesUri );

  // 1.0 get additional coverage info
  if ( !mCapabilities.describeCoverage( mIdentifier ) )
  {
    appendError( ERR( tr( "Cannot describe coverage" ) ) );
    return;
  }

  mCoverageSummary = mCapabilities.coverage( mIdentifier );
  if ( !mCoverageSummary.valid )
  {
    // Should not happen if describeCoverage() did not fail
    appendError( ERR( tr( "Coverage not found" ) ) );
    return;
  }

  // It may happen that format is empty (e.g. uri created in python script),
  // in that casei select one from available formats
  if ( mFormat.isEmpty() )
  {
    // TIFF is known by GDAL
    mFormat = mCoverageSummary.supportedFormat.filter( QStringLiteral( "tif" ), Qt::CaseInsensitive ).value( 0 );
  }
  if ( mFormat.isEmpty() )
  {
    // Take the first if TIFF was not found
    mFormat = mCoverageSummary.supportedFormat.value( 0 );
  }

  // We cannot continue without format, it is required
  if ( mFormat.isEmpty() ) return;

  // It could happen (usually not with current QgsWCSSourceSelect if at least
  // one CRS is available) that crs is not set in uri, in that case we
  // use the native, if available or WGS84 or the first supported
  if ( mCoverageCrs.isEmpty() )
  {
    QgsDebugMsg( "nativeCrs = " + mCoverageSummary.nativeCrs );
    QgsDebugMsg( "supportedCrs = " + mCoverageSummary.supportedCrs.join( "," ) );
    if ( !mCoverageSummary.nativeCrs.isEmpty() )
    {
      setCoverageCrs( mCoverageSummary.nativeCrs );
    }
    else if ( mCoverageSummary.supportedCrs.contains( QStringLiteral( "EPSG:4326" ), Qt::CaseInsensitive ) )
    {
      setCoverageCrs( QStringLiteral( "EPSG:4326" ) );
    }
    else if ( !mCoverageSummary.supportedCrs.isEmpty() )
    {
      setCoverageCrs( mCoverageSummary.supportedCrs.value( 0 ) );
    }
  }
  QgsDebugMsg( "mCoverageCrs = " + mCoverageCrs );

  // It may happen that coverage CRS is not given or it is unknown
  // in that case we continue without CRS and user is asked for it
  //if ( mCoverageCrs.isEmpty() ) return;

  // Native size
  mWidth = mCoverageSummary.width;
  mHeight = mCoverageSummary.height;
  mHasSize = mCoverageSummary.hasSize;

  QgsDebugMsg( QStringLiteral( "mWidth = %1 mHeight = %2" ).arg( mWidth ).arg( mHeight ) );

  // TODO: Consider if/how to recalculate mWidth, mHeight if non native CRS is used

  if ( !calculateExtent() )
  {
    appendError( ERR( tr( "Cannot calculate extent" ) ) );
    return;
  }

  // Get small piece of coverage to find GDAL data type and number of bands
  int bandNo = 0; // All bands
  int width;
  int height;
  QString crs;
  QgsRectangle box; // box to use to calc extent
  // Prefer native CRS
  if ( !mCoverageSummary.nativeCrs.isEmpty() &&
       !mCoverageSummary.nativeBoundingBox.isEmpty() &&
       mCoverageSummary.supportedCrs.contains( mCoverageSummary.nativeCrs ) &&
       mHasSize )
  {
    box = mCoverageSummary.nativeBoundingBox;
    width = mWidth;
    height = mHeight;
    crs = mCoverageSummary.nativeCrs;
  }
  else
  {
    box = mCoverageExtent;
    if ( mHasSize )
    {
      width = mWidth;
      height = mHeight;
    }
    else
    {
      // just a number to get smaller piece of coverage
      width = 1000;
      height = 1000;
    }
  }
  double xRes = box.width() / width;
  double yRes = box.height() / height;
  QgsPointXY p = box.center();

  // width and height different to recognize rotation
  int requestWidth = 6;
  int requestHeight = 3;

  // extent to be used for test request
  double halfWidth = xRes * ( requestWidth / 2. );
  double halfHeight = yRes * ( requestHeight / 2. );

  // Using non native CRS (if we don't know which is native) it could easily happen,
  // that a small part of bbox in request CRS near margin falls outside
  // coverage native bbox and server reports error => take a piece from center
  QgsRectangle extent = QgsRectangle( p.x() - halfWidth, p.y() - halfHeight, p.x() + halfWidth, p.y() + halfHeight );

  getCache( bandNo, extent, requestWidth, requestHeight, crs );

  if ( !mCachedGdalDataset )
  {
    setError( mCachedError );
    appendError( ERR( tr( "Cannot get test dataset." ) ) );
    return;
  }

  mBandCount = GDALGetRasterCount( mCachedGdalDataset.get() );
  QgsDebugMsg( QStringLiteral( "mBandCount = %1" ).arg( mBandCount ) );

  // Check for server particularities (bbox, rotation)
  int responseWidth = GDALGetRasterXSize( mCachedGdalDataset.get() );
  int responseHeight = GDALGetRasterYSize( mCachedGdalDataset.get() );

  QgsDebugMsg( QStringLiteral( "requestWidth = %1 requestHeight = %2 responseWidth = %3 responseHeight = %4)" ).arg( requestWidth ).arg( requestHeight ).arg( responseWidth ).arg( responseHeight ) );
  // GeoServer and ArcGIS are using for 1.1 box "pixel" edges
  // Mapserver is using pixel centers according to 1.1. specification
  if ( ( responseWidth == requestWidth - 1 && responseHeight == requestHeight - 1 ) ||
       ( responseWidth == requestHeight - 1 && responseHeight == requestWidth - 1 ) )
  {
    mFixBox = true;
    QgsDebugMsg( QStringLiteral( "Test response size is smaller by pixel, using mFixBox" ) );
  }
  // Geoserver is giving rotated raster for geographic CRS - switched axis,
  // Geoserver developers argue that changed axis order applies also to
  // returned raster, that is exaggerated IMO but we have to handle that.
  if ( ( responseWidth == requestHeight && responseHeight == requestWidth ) ||
       ( responseWidth == requestHeight - 1 && responseHeight == requestWidth - 1 ) )
  {
    mFixRotate = true;
    QgsDebugMsg( QStringLiteral( "Test response is rotated, using mFixRotate" ) );
  }

  // Get types
  // TODO: we are using the same data types like GDAL (not wider like GDAL provider)
  // with expectation to replace 'no data' values by NaN
  mSrcGdalDataType.reserve( mBandCount );
  mGdalDataType.reserve( mBandCount );
  for ( int i = 1; i <= mBandCount; i++ )
  {
    GDALRasterBandH gdalBand = GDALGetRasterBand( mCachedGdalDataset.get(), i );
    GDALDataType myGdalDataType = GDALGetRasterDataType( gdalBand );

    QgsDebugMsg( QStringLiteral( "myGdalDataType[%1] = %2" ).arg( i - 1 ).arg( myGdalDataType ) );
    mSrcGdalDataType.append( myGdalDataType );

    // UMN Mapserver does not automatically set null value, METADATA wcs_rangeset_nullvalue must be used
    // http://lists.osgeo.org/pipermail/mapserver-users/2010-April/065328.html

    // TODO: This could be shared with GDAL provider
    int isValid = false;
    double myNoDataValue = GDALGetRasterNoDataValue( gdalBand, &isValid );
    if ( isValid )
    {
      QgsDebugMsg( QStringLiteral( "GDALGetRasterNoDataValue = %1" ).arg( myNoDataValue ) );
      myNoDataValue = QgsRaster::representableValue( myNoDataValue, dataTypeFromGdal( myGdalDataType ) );
      mSrcNoDataValue.append( myNoDataValue );
      mSrcHasNoDataValue.append( true );
      mUseSrcNoDataValue.append( true );
    }
    else
    {
      mSrcNoDataValue.append( std::numeric_limits<double>::quiet_NaN() );
      mSrcHasNoDataValue.append( false );
      mUseSrcNoDataValue.append( false );
    }
    // It may happen that nodata value given by GDAL is wrong and it has to be
    // disabled by user, in that case we need another value to be used for nodata
    // (for reprojection for example) -> always internaly represent as wider type
    // with mInternalNoDataValue in reserve.
    // No retyping, no internal values for now
#if 0
    int myInternalGdalDataType = myGdalDataType;
    double myInternalNoDataValue;
    switch ( srcDataType( i ) )
    {
      case Qgis::Byte:
        myInternalNoDataValue = -32768.0;
        myInternalGdalDataType = GDT_Int16;
        break;
      case Qgis::Int16:
        myInternalNoDataValue = -2147483648.0;
        myInternalGdalDataType = GDT_Int32;
        break;
      case Qgis::UInt16:
        myInternalNoDataValue = -2147483648.0;
        myInternalGdalDataType = GDT_Int32;
        break;
      case Qgis::Int32:
        // We believe that such values is no used in real data
        myInternalNoDataValue = -2147483648.0;
        break;
      case Qgis::UInt32:
        // We believe that such values is no used in real data
        myInternalNoDataValue = 4294967295.0;
        break;
      default: // Float32, Float64
        //myNoDataValue = std::numeric_limits<int>::max();
        // NaN should work well
        myInternalNoDataValue = std::numeric_limits<double>::quiet_NaN();
    }
    mGdalDataType.append( myInternalGdalDataType );
    mInternalNoDataValue.append( myInternalNoDataValue );
    QgsDebugMsg( QStringLiteral( "mInternalNoDataValue[%1] = %2" ).arg( i - 1 ).arg( mInternalNoDataValue[i - 1] ) );
#endif
    mGdalDataType.append( myGdalDataType );

#if 0
    // TODO: what to do if null values from DescribeCoverage differ?
    if ( !mCoverageSummary.nullValues.contains( myNoDataValue ) )
    {
      QgsDebugMsg( QStringLiteral( "noDataValue %1 is missing in nullValues from CoverageDescription" ).arg( myNoDataValue ) );
    }
#endif

    QgsDebugMsg( QStringLiteral( "mSrcGdalDataType[%1] = %2" ).arg( i - 1 ).arg( mSrcGdalDataType.at( i - 1 ) ) );
    QgsDebugMsg( QStringLiteral( "mGdalDataType[%1] = %2" ).arg( i - 1 ).arg( mGdalDataType.at( i - 1 ) ) );
    QgsDebugMsg( QStringLiteral( "mSrcNoDataValue[%1] = %2" ).arg( i - 1 ).arg( mSrcNoDataValue.at( i - 1 ) ) );

    // Create and store color table
    // TODO: never tested because mapserver (6.0.3) does not support color tables
    mColorTables.append( QgsGdalProviderBase::colorTable( mCachedGdalDataset.get(), i ) );
  }

  clearCache();

  // Block size is used for for statistics
  // TODO: How to find maximum block size supported by server?
  if ( mHasSize )
  {
    // This is taken from GDAL, how they come to these numbers?
    if ( mWidth > 1800 ) mXBlockSize = 1024;
    else mXBlockSize = mWidth;

    if ( mHeight > 900 ) mYBlockSize = 512;
    else mYBlockSize = mHeight;
  }

  mValid = true;
  QgsDebugMsg( QStringLiteral( "Constructed OK, provider valid." ) );
}

bool QgsWcsProvider::parseUri( const QString &uriString )
{

  QgsDebugMsg( "uriString = " + uriString );
  QgsDataSourceUri uri;
  uri.setEncodedUri( uriString );

  mMaxWidth = 0;
  mMaxHeight = 0;

  mHttpUri = uri.param( QStringLiteral( "url" ) );
  mBaseUrl = prepareUri( mHttpUri );
  QgsDebugMsg( "mBaseUrl = " + mBaseUrl );

  mIgnoreGetCoverageUrl = uri.hasParam( QStringLiteral( "IgnoreGetMapUrl" ) );
  mIgnoreAxisOrientation = uri.hasParam( QStringLiteral( "IgnoreAxisOrientation" ) ); // must be before parsing!
  mInvertAxisOrientation = uri.hasParam( QStringLiteral( "InvertAxisOrientation" ) ); // must be before parsing!

  mAuth.mUserName = uri.param( QStringLiteral( "username" ) );
  QgsDebugMsg( "set username to " + mAuth.mUserName );

  mAuth.mPassword = uri.param( QStringLiteral( "password" ) );
  QgsDebugMsg( "set password to " + mAuth.mPassword );

  if ( uri.hasParam( QStringLiteral( "authcfg" ) ) )
  {
    mAuth.mAuthCfg = uri.param( QStringLiteral( "authcfg" ) );
  }
  QgsDebugMsg( "set authcfg to " + mAuth.mAuthCfg );

  mIdentifier = uri.param( QStringLiteral( "identifier" ) );

  mTime = uri.param( QStringLiteral( "time" ) );

  setFormat( uri.param( QStringLiteral( "format" ) ) );

  if ( !uri.param( QStringLiteral( "crs" ) ).isEmpty() )
  {
    setCoverageCrs( uri.param( QStringLiteral( "crs" ) ) );
  }

  QString cache = uri.param( QStringLiteral( "cache" ) );
  if ( !cache.isEmpty() )
  {
    mCacheLoadControl = QgsNetworkAccessManager::cacheLoadControlFromName( cache );
  }
  QgsDebugMsg( QStringLiteral( "mCacheLoadControl = %1" ).arg( mCacheLoadControl ) );

  return true;
}

QString QgsWcsProvider::prepareUri( QString uri ) const
{
  if ( !uri.contains( '?' ) )
  {
    uri.append( '?' );
  }
  else if ( uri.right( 1 ) != QLatin1String( "?" ) && uri.right( 1 ) != QLatin1String( "&" ) )
  {
    uri.append( '&' );
  }

  return uri;
}

QgsWcsProvider::~QgsWcsProvider()
{
  QgsDebugMsg( QStringLiteral( "deconstructing." ) );

  // Dispose of any cached image as created by draw()
  clearCache();
}

QgsWcsProvider *QgsWcsProvider::clone() const
{
  QgsDataProvider::ProviderOptions providerOptions;
  QgsWcsProvider *provider = new QgsWcsProvider( dataSourceUri(), providerOptions );
  provider->copyBaseSettings( *this );
  return provider;
}

QString QgsWcsProvider::baseUrl() const
{
  return mBaseUrl;
}

QString QgsWcsProvider::format() const
{
  return mFormat;
}

void QgsWcsProvider::setFormat( QString const &format )
{
  QgsDebugMsg( "Setting format to " + format + '.' );
  mFormat = format;
}


void QgsWcsProvider::setCoverageCrs( QString const &crs )
{
  QgsDebugMsg( "Setting coverage CRS to " + crs + '.' );

  if ( crs != mCoverageCrs && !crs.isEmpty() )
  {
    // delete old coordinate transform as it is no longer valid
    mCoordinateTransform = QgsCoordinateTransform();

    mExtentDirty = true;

    mCoverageCrs = crs;

    mCrs = QgsCoordinateReferenceSystem::fromOgcWmsCrs( mCoverageCrs );
  }
}

void QgsWcsProvider::setQueryItem( QUrl &url, const QString &item, const QString &value ) const
{
  url.removeQueryItem( item );
  url.addQueryItem( item, value );
}

void QgsWcsProvider::readBlock( int bandNo, QgsRectangle  const &viewExtent, int pixelWidth, int pixelHeight, void *block, QgsRasterBlockFeedback *feedback )
{
  // TODO: set block to null values, move that to function and call only if fails
  memset( block, 0, pixelWidth * pixelHeight * QgsRasterBlock::typeSize( dataType( bandNo ) ) );

  // Requested extent must at least partially overlap coverage extent, otherwise
  // server gives error. QGIS usually does not request blocks outside raster extent
  // (higher level checks) but it is better to do check here as well
  if ( !viewExtent.intersects( mCoverageExtent ) )
  {
    return;
  }

  // Can we reuse the previously cached coverage?
  if ( !mCachedGdalDataset ||
       mCachedViewExtent != viewExtent ||
       mCachedViewWidth != pixelWidth ||
       mCachedViewHeight != pixelHeight )
  {
    getCache( bandNo, viewExtent, pixelWidth, pixelHeight, QString(), feedback );
  }

  if ( mCachedGdalDataset )
  {
    // It may happen (Geoserver) that if requested BBOX is larger than coverage
    // extent, the returned data cover intersection of requested BBOX and coverage
    // extent scaled to requested WIDTH/HEIGHT => check extent
    // Unfortunately if the received raster does not have a CRS, the extent is the raster size
    // and in that case it cannot be used to verify extent
    QgsCoordinateReferenceSystem cacheCrs;
    if ( !cacheCrs.createFromWkt( GDALGetProjectionRef( mCachedGdalDataset.get() ) ) &&
         !cacheCrs.createFromWkt( GDALGetGCPProjection( mCachedGdalDataset.get() ) ) )
    {
      QgsDebugMsg( QStringLiteral( "Cached does not have CRS" ) );
    }
    QgsDebugMsg( "Cache CRS: " + cacheCrs.authid() + ' ' + cacheCrs.description() );

    QgsRectangle cacheExtent = QgsGdalProviderBase::extent( mCachedGdalDataset.get() );
    QgsDebugMsg( "viewExtent = " + viewExtent.toString() );
    QgsDebugMsg( "cacheExtent = " + cacheExtent.toString() );
    // TODO: check also rotated
    if ( cacheCrs.isValid() && !mFixRotate )
    {
      // using qgsDoubleNear is too precise, example accetable difference:
      // 179.9999999306699863 x 179.9999999306700431
      if ( !qgsDoubleNearSig( cacheExtent.xMinimum(), viewExtent.xMinimum(), 10 ) ||
           !qgsDoubleNearSig( cacheExtent.yMinimum(), viewExtent.yMinimum(), 10 ) ||
           !qgsDoubleNearSig( cacheExtent.xMaximum(), viewExtent.xMaximum(), 10 ) ||
           !qgsDoubleNearSig( cacheExtent.yMaximum(), viewExtent.yMaximum(), 10 ) )
      {
        QgsDebugMsg( QStringLiteral( "cacheExtent and viewExtent differ" ) );
        QgsMessageLog::logMessage( tr( "Received coverage has wrong extent %1 (expected %2)" ).arg( cacheExtent.toString(), viewExtent.toString() ), tr( "WCS" ) );
        // We are doing all possible to avoid this situation,
        // If it happens, it would be possible to rescale the portion we get
        // to only part of the data block, but it is better to left it
        // blank, so that the problem may be discovered in its origin.
        return;
      }
    }

    int width = GDALGetRasterXSize( mCachedGdalDataset.get() );
    int height = GDALGetRasterYSize( mCachedGdalDataset.get() );
    QgsDebugMsg( QStringLiteral( "cached data width = %1 height = %2 (expected %3 x %4)" ).arg( width ).arg( height ).arg( pixelWidth ).arg( pixelHeight ) );

    GDALRasterBandH gdalBand = GDALGetRasterBand( mCachedGdalDataset.get(), bandNo );
    // TODO: check type?, check band count?
    if ( mFixRotate && width == pixelHeight && height == pixelWidth )
    {
      // Rotate counter clockwise
      // If GridOffsets With GeoServer,
      QgsDebugMsg( QStringLiteral( "Rotating raster" ) );
      int pixelSize = QgsRasterBlock::typeSize( dataType( bandNo ) );
      QgsDebugMsg( QStringLiteral( "pixelSize = %1" ).arg( pixelSize ) );
      int size = width * height * pixelSize;
      void *tmpData = malloc( size );
      if ( ! tmpData )
      {
        QgsDebugMsg( QStringLiteral( "Couldn't allocate memory of %1 bytes" ).arg( size ) );
        return;
      }
      if ( GDALRasterIO( gdalBand, GF_Read, 0, 0, width, height, tmpData, width, height, ( GDALDataType ) mGdalDataType.at( bandNo - 1 ), 0, 0 ) != CE_None )
      {
        QgsDebugMsg( QStringLiteral( "Raster IO Error" ) );
      }
      for ( int i = 0; i < pixelHeight; i++ )
      {
        for ( int j = 0; j < pixelWidth; j++ )
        {
          int destIndex = pixelSize * ( i * pixelWidth + j );
          int srcIndex = pixelSize * ( j * width + ( width - i - 1 ) );
          memcpy( ( char * )block + destIndex, ( char * )tmpData + srcIndex, pixelSize );
        }
      }
      free( tmpData );
    }
    else if ( width == pixelWidth && height == pixelHeight )
    {
      if ( GDALRasterIO( gdalBand, GF_Read, 0, 0, pixelWidth, pixelHeight, block, pixelWidth, pixelHeight, ( GDALDataType ) mGdalDataType.at( bandNo - 1 ), 0, 0 ) != CE_None )
      {
        QgsDebugMsg( QStringLiteral( "Raster IO Error" ) );
      }
      else
      {
        QgsDebugMsg( QStringLiteral( "Block read OK" ) );
      }
    }
    else
    {
      // This should not happen, but it is better to give distorted result + warning
      if ( GDALRasterIO( gdalBand, GF_Read, 0, 0, width, height, block, pixelWidth, pixelHeight, ( GDALDataType ) mGdalDataType.at( bandNo - 1 ), 0, 0 ) != CE_None )
      {
        QgsDebugMsg( QStringLiteral( "Raster IO Error" ) );
      }
      QgsMessageLog::logMessage( tr( "Received coverage has wrong size %1 x %2 (expected %3 x %4)" ).arg( width ).arg( height ).arg( pixelWidth ).arg( pixelHeight ), tr( "WCS" ) );
    }
  }
}

void QgsWcsProvider::getCache( int bandNo, QgsRectangle  const &viewExtent, int pixelWidth, int pixelHeight, QString crs, QgsRasterBlockFeedback *feedback ) const
{
  Q_UNUSED( bandNo );
  // delete cached data
  clearCache();

  if ( crs.isEmpty() )
  {
    crs = mCoverageCrs;
  }

  mCachedViewExtent = viewExtent;
  mCachedViewWidth = pixelWidth;
  mCachedViewHeight = pixelHeight;

  // --------------- BOUNDING BOX --------------------------------
  //according to the WCS spec for 1.1, some CRS have inverted axis
  // box:
  //  1.0.0: minx,miny,maxx,maxy
  //  1.1.0, 1.1.2: OGC 07-067r5 (WCS 1.1.2) referes to OGC 06-121r3 which says
  //  "The number of axes included, and the order of these axes, shall be as specified
  //  by the referenced CRS." That means inverted for geographic.
  bool changeXY = false;
  if ( !mIgnoreAxisOrientation && ( mCapabilities.version().startsWith( QLatin1String( "1.1" ) ) ) )
  {
    //create CRS from string
    QgsCoordinateReferenceSystem srs = QgsCoordinateReferenceSystem::fromOgcWmsCrs( crs );
    if ( srs.isValid() && srs.hasAxisInverted() )
    {
      changeXY = true;
    }
  }

  if ( mInvertAxisOrientation ) changeXY = !changeXY;

  double xRes = viewExtent.width() / pixelWidth;
  double yRes = viewExtent.height() / pixelHeight;
  QgsRectangle extent = viewExtent;
  // WCS 1.1 grid is using grid points (surrounded by sample spaces) and
  // "The spatial extent of a grid coverage extends only as far as the outermost
  // grid points contained in the bounding box. It does NOT include any area
  // (partial or whole grid cells or sample spaces) beyond those grid points."
  // Mapserver and GDAL are using bbox defined by grid points, i.e. shrunk
  // by 1 pixel, but Geoserver and ArcGIS are using full bbox including
  // the space around edge grid points.
  if ( mCapabilities.version().startsWith( QLatin1String( "1.1" ) ) && !mFixBox )
  {
    // shrink the extent to border cells centers by half cell size
    extent = QgsRectangle( viewExtent.xMinimum() + xRes / 2., viewExtent.yMinimum() + yRes / 2., viewExtent.xMaximum() - xRes / 2., viewExtent.yMaximum() - yRes / 2. );
  }

  // Bounding box in WCS format (Warning: does not work with scientific notation)
  QString bbox = QString( changeXY ? "%2,%1,%4,%3" : "%1,%2,%3,%4" )
                 .arg( qgsDoubleToString( extent.xMinimum() ),
                       qgsDoubleToString( extent.yMinimum() ),
                       qgsDoubleToString( extent.xMaximum() ),
                       qgsDoubleToString( extent.yMaximum() ) );

  QUrl url( mIgnoreGetCoverageUrl ? mBaseUrl : mCapabilities.getCoverageUrl() );

  // Version 1.0.0, 1.1.0, 1.1.2
  setQueryItem( url, QStringLiteral( "SERVICE" ), QStringLiteral( "WCS" ) );
  setQueryItem( url, QStringLiteral( "VERSION" ), mCapabilities.version() );
  setQueryItem( url, QStringLiteral( "REQUEST" ), QStringLiteral( "GetCoverage" ) );
  setQueryItem( url, QStringLiteral( "FORMAT" ), mFormat );

  // Version 1.0.0
  if ( mCapabilities.version().startsWith( QLatin1String( "1.0" ) ) )
  {
    setQueryItem( url, QStringLiteral( "COVERAGE" ), mIdentifier );
    if ( !mTime.isEmpty() )
    {
      // It seems that Mmapserver (6.0.3) WCS 1.1 completely ignores
      // TemporalDomain. Some code (copy-pasted from 1.0) is commented in
      // msWCSDescribeCoverage_CoverageDescription11() and GetCoverage
      // TimeSequence param is not supported at all. If a coverage is defined
      // with timeposition in mapfile, the result of GetCoverage is empty
      // raster (all values 0).
      setQueryItem( url, QStringLiteral( "TIME" ), mTime );
    }
    setQueryItem( url, QStringLiteral( "BBOX" ), bbox );
    setQueryItem( url, QStringLiteral( "CRS" ), crs ); // request BBOX CRS
    setQueryItem( url, QStringLiteral( "RESPONSE_CRS" ), crs ); // response CRS
    setQueryItem( url, QStringLiteral( "WIDTH" ), QString::number( pixelWidth ) );
    setQueryItem( url, QStringLiteral( "HEIGHT" ), QString::number( pixelHeight ) );
  }

  // Version 1.1.0, 1.1.2
  if ( mCapabilities.version().startsWith( QLatin1String( "1.1" ) ) )
  {
    setQueryItem( url, QStringLiteral( "IDENTIFIER" ), mIdentifier );
    QString crsUrn = QStringLiteral( "urn:ogc:def:crs:%1::%2" ).arg( crs.split( ':' ).value( 0 ), crs.split( ':' ).value( 1 ) );
    bbox += ',' + crsUrn;

    if ( !mTime.isEmpty() )
    {
      setQueryItem( url, QStringLiteral( "TIMESEQUENCE" ), mTime );
    }

    setQueryItem( url, QStringLiteral( "BOUNDINGBOX" ), bbox );

    //  Example:
    //   GridBaseCRS=urn:ogc:def:crs:SG:6.6:32618
    //   GridType=urn:ogc:def:method:CS:1.1:2dGridIn2dCrs
    //   GridCS=urn:ogc:def:cs:OGC:0Grid2dSquareCS
    //   GridOrigin=0,0
    //   GridOffsets=0.0707,-0.0707,0.1414,0.1414&

    setQueryItem( url, QStringLiteral( "GRIDBASECRS" ), crsUrn ); // response CRS

    setQueryItem( url, QStringLiteral( "GRIDCS" ), QStringLiteral( "urn:ogc:def:cs:OGC:0.0:Grid2dSquareCS" ) );

    setQueryItem( url, QStringLiteral( "GRIDTYPE" ), QStringLiteral( "urn:ogc:def:method:WCS:1.1:2dSimpleGrid" ) );

    // GridOrigin is BBOX minx, maxy
    // Note: shifting origin to cell center (not really necessary nor making sense)
    // does not work with Mapserver 6.0.3
    // Mapserver 6.0.3 does not work with origin on yMinimum (lower left)
    // Geoserver works OK with yMinimum (lower left)
    QString gridOrigin = QString( changeXY ? "%2,%1" : "%1,%2" )
                         .arg( qgsDoubleToString( extent.xMinimum() ),
                               qgsDoubleToString( extent.yMaximum() ) );
    setQueryItem( url, QStringLiteral( "GRIDORIGIN" ), gridOrigin );

    // GridOffsets WCS 1.1:
    // GridType urn:ogc:def:method:WCS:1.1:2dSimpleGrid : 2 values
    // GridType urn:ogc:def:method:WCS:1.1:2dGridIn2dCrs : 4 values, the center two of these offsets will be zero when the GridCRS is not rotated or skewed in the GridBaseCRS.
    // The yOff must be negative because origin is in upper corner
    // WCS 1.1.2: BaseX = origin(1) + offsets(1) * GridX
    //            BaseY = origin(2) + offsets(2) * GridY
    // otherwise GeoServer gives mirrored result, however
    // Mapserver works OK with both positive and negative
    // OTOH, the yOff offset must not be negative with mFixRotate and Geoserver 2.1-SNAPSHOT
    // but it must be negative with GeoServer 2.1.3 and mFixRotate. I am not sure
    // at this moment 100% -> disabling positive yOff for now - TODO: try other servers
    //double yOff = mFixRotate ? yRes : -yRes; // this was working with some servers I think
    double yOff = -yRes;
    QString gridOffsets = QString( changeXY ? "%2,%1" : "%1,%2" )
                          //setQueryItem( url, "GRIDTYPE", "urn:ogc:def:method:WCS:1.1:2dGridIn2dCrs" );
                          //QString gridOffsets = QString( changeXY ? "%2,0,0,%1" : "%1,0,0,%2" )
                          .arg( qgsDoubleToString( xRes ),
                                qgsDoubleToString( yOff ) );
    setQueryItem( url, QStringLiteral( "GRIDOFFSETS" ), gridOffsets );
  }

  QgsDebugMsg( QStringLiteral( "GetCoverage: %1" ).arg( url.toString() ) );

  // cache some details for if the user wants to do an identifyAsHtml() later

  //mGetFeatureInfoUrlBase = mIgnoreGetFeatureInfoUrl ? mBaseUrl : getFeatureInfoUrl();

  emit statusChanged( tr( "Getting map via WCS." ) );

  QgsWcsDownloadHandler handler( url, mAuth, mCacheLoadControl, mCachedData, mCapabilities.version(), mCachedError, feedback );
  handler.blockingDownload();

  QgsDebugMsg( QStringLiteral( "%1 bytes received" ).arg( mCachedData.size() ) );
  if ( mCachedData.isEmpty() )
  {
    QgsMessageLog::logMessage( tr( "No data received" ), tr( "WCS" ) );
    clearCache();
    return;
  }

#if 0
  QFile myFile( "/tmp/qgiswcscache.dat" );
  if ( myFile.open( QIODevice::WriteOnly ) )
  {
    QDataStream myStream( &myFile );
    myStream.writeRawData( mCachedData.data(), mCachedData.size() );
    myFile.close();
  }
#endif

  mCachedMemFile = VSIFileFromMemBuffer( mCachedMemFilename.toUtf8().constData(),
                                         ( GByte * )mCachedData.data(),
                                         ( vsi_l_offset )mCachedData.size(),
                                         FALSE );

  if ( !mCachedMemFile )
  {
    QgsMessageLog::logMessage( tr( "Cannot create memory file" ), tr( "WCS" ) );
    clearCache();
    return;
  }
  QgsDebugMsg( QStringLiteral( "Memory file created" ) );

  CPLErrorReset();
  mCachedGdalDataset.reset( GDALOpen( mCachedMemFilename.toUtf8().constData(), GA_ReadOnly ) );
  if ( !mCachedGdalDataset )
  {
    QgsMessageLog::logMessage( QString::fromUtf8( CPLGetLastErrorMsg() ), tr( "WCS" ) );
    clearCache();
    return;
  }
  QgsDebugMsg( QStringLiteral( "Dataset opened" ) );

}

// For stats only, maybe change QgsRasterDataProvider::bandStatistics() to
// use standard readBlock with extent
void QgsWcsProvider::readBlock( int bandNo, int xBlock, int yBlock, void *block )
{

  QgsDebugMsg( QStringLiteral( "xBlock = %1 yBlock = %2" ).arg( xBlock ).arg( yBlock ) );

  if ( !mHasSize ) return;

  double xRes = mCoverageExtent.width() / mWidth;
  double yRes = mCoverageExtent.height() / mHeight;

  // blocks on edges may run out of extent, that should not be problem (at least for
  // stats - there is a check for it)
  double xMin = mCoverageExtent.xMinimum() + xRes * xBlock * mXBlockSize;
  double xMax = xMin + xRes * mXBlockSize;
  double yMax = mCoverageExtent.yMaximum() - yRes * yBlock * mYBlockSize;
  double yMin = yMax - yRes * mYBlockSize;
  //QgsDebugMsg( QStringLiteral("yMin = %1 yMax = %2").arg(yMin).arg(yMax) );

  QgsRectangle extent( xMin, yMin, xMax, yMax );

  readBlock( bandNo, extent, mXBlockSize, mYBlockSize, block, nullptr );
}


// This could be shared with GDAL provider
Qgis::DataType QgsWcsProvider::sourceDataType( int bandNo ) const
{
  if ( bandNo < 0 || bandNo > mSrcGdalDataType.size() )
  {
    return Qgis::UnknownDataType;
  }

  return dataTypeFromGdal( mSrcGdalDataType[bandNo - 1] );
}

Qgis::DataType QgsWcsProvider::dataType( int bandNo ) const
{
  if ( bandNo < 0 || bandNo > mGdalDataType.size() )
  {
    return Qgis::UnknownDataType;
  }

  return dataTypeFromGdal( mGdalDataType[bandNo - 1] );
}

int QgsWcsProvider::bandCount() const
{
  return mBandCount;
}

// this is only called once when statistics are calculated
// TODO
int QgsWcsProvider::xBlockSize() const
{
  return mXBlockSize;
}
int QgsWcsProvider::yBlockSize() const
{
  return mYBlockSize;
}

int QgsWcsProvider::xSize() const { return mWidth; }
int QgsWcsProvider::ySize() const { return mHeight; }

void QgsWcsProvider::clearCache() const
{
  if ( mCachedGdalDataset )
  {
    QgsDebugMsg( QStringLiteral( "Close mCachedGdalDataset" ) );
    mCachedGdalDataset.reset();
    QgsDebugMsg( QStringLiteral( "Closed" ) );
  }
  if ( mCachedMemFile )
  {
    QgsDebugMsg( QStringLiteral( "Close mCachedMemFile" ) );
    VSIFCloseL( mCachedMemFile );
    mCachedMemFile = nullptr;
    QgsDebugMsg( QStringLiteral( "Closed" ) );
  }
  QgsDebugMsg( QStringLiteral( "Clear mCachedData" ) );
  mCachedData.clear();
  mCachedError.clear();
  QgsDebugMsg( QStringLiteral( "Cleared" ) );
}

QList<QgsColorRampShader::ColorRampItem> QgsWcsProvider::colorTable( int bandNumber )const
{
  return mColorTables.value( bandNumber - 1 );
}

int QgsWcsProvider::colorInterpretation( int bandNo ) const
{
  GDALRasterBandH myGdalBand = GDALGetRasterBand( mCachedGdalDataset.get(), bandNo );
  return colorInterpretationFromGdal( GDALGetRasterColorInterpretation( myGdalBand ) );
}


bool QgsWcsProvider::parseServiceExceptionReportDom( QByteArray const &xml, const QString &wcsVersion, QString &errorTitle, QString &errorText )
{

#ifdef QGISDEBUG
  //test the content of the QByteArray
  QString responsestring( xml );
  QgsDebugMsg( "received the following data: " + responsestring );
#endif

  // Convert completed document into a Dom
  QDomDocument doc;
  QString errorMsg;
  int errorLine;
  int errorColumn;
  bool contentSuccess = doc.setContent( xml, false, &errorMsg, &errorLine, &errorColumn );

  if ( !contentSuccess )
  {
    errorTitle = tr( "Dom Exception" );
    errorText = tr( "Could not get WCS Service Exception at %1 at line %2 column %3\n\nResponse was:\n\n%4" )
                .arg( errorMsg )
                .arg( errorLine )
                .arg( errorColumn )
                .arg( QString( xml ) );

    QgsLogger::debug( "Dom Exception: " + errorText );

    return false;
  }

  QDomElement docElem = doc.documentElement();

  // TODO: Assert the docElem.tagName() is
  //  ServiceExceptionReport // 1.0
  //  ows:ExceptionReport  // 1.1

  //QString version = docElem.attribute("version");

  QDomElement e;
  if ( wcsVersion.startsWith( QLatin1String( "1.0" ) ) )
  {
    e = QgsWcsCapabilities::domElement( docElem, QStringLiteral( "ServiceException" ) );
  }
  else // 1.1
  {
    e = QgsWcsCapabilities::domElement( docElem, QStringLiteral( "Exception" ) );
  }
  parseServiceException( e, wcsVersion, errorTitle, errorText );

  QgsDebugMsg( QStringLiteral( "exiting." ) );

  return true;
}

void QgsWcsProvider::parseServiceException( QDomElement const &e, const QString &wcsVersion, QString &errorTitle, QString &errorText )
{

  errorTitle = tr( "Service Exception" );

  QMap<QString, QString> exceptions;

  // Some codes are in both 1.0 and 1.1, in that case the 'meaning of code' is
  // taken from 1.0 specification

  // set up friendly descriptions for the service exception
  // 1.0
  exceptions[QStringLiteral( "InvalidFormat" )] = tr( "Request contains a format not offered by the server." );
  exceptions[QStringLiteral( "CoverageNotDefined" )] = tr( "Request is for a Coverage not offered by the service instance." );
  exceptions[QStringLiteral( "CurrentUpdateSequence" )] = tr( "Value of (optional) UpdateSequence parameter in GetCapabilities request is equal to current value of service metadata update sequence number." );
  exceptions[QStringLiteral( "InvalidUpdateSequence" )] = tr( "Value of (optional) UpdateSequence parameter in GetCapabilities request is greater than current value of service metadata update sequence number." );
  // 1.0, 1.1
  exceptions[QStringLiteral( "MissingParameterValue" )] = tr( "Request does not include a parameter value, and the server instance did not declare a default value for that dimension." );
  exceptions[QStringLiteral( "InvalidParameterValue" )] = tr( "Request contains an invalid parameter value." );
  // 1.1
  exceptions[QStringLiteral( "NoApplicableCode" )] = tr( "No other exceptionCode specified by this service and server applies to this exception." );
  exceptions[QStringLiteral( "UnsupportedCombination" )] = tr( "Operation request contains an output CRS that can not be used within the output format." );
  exceptions[QStringLiteral( "NotEnoughStorage" )] = tr( "Operation request specifies to \"store\" the result, but not enough storage is available to do this." );

  QString seCode;
  QString seText;
  if ( wcsVersion.startsWith( QLatin1String( "1.0" ) ) )
  {
    seCode = e.attribute( QStringLiteral( "code" ) );
    seText = e.text();
  }
  else
  {
    QStringList codes;
    seCode = e.attribute( QStringLiteral( "exceptionCode" ) );
    // UMN Mapserver (6.0.3) has messed/switched 'locator' and 'exceptionCode'
    if ( ! exceptions.contains( seCode ) )
    {
      seCode = e.attribute( QStringLiteral( "locator" ) );
      if ( ! exceptions.contains( seCode ) )
      {
        seCode.clear();
      }
    }
    seText = QgsWcsCapabilities::firstChildText( e, QStringLiteral( "ExceptionText" ) );
  }

  if ( seCode.isEmpty() )
  {
    errorText = tr( "(No error code was reported)" );
  }
  else if ( exceptions.contains( seCode ) )
  {
    errorText = exceptions.value( seCode );
  }
  else
  {
    errorText = seCode + ' ' + tr( "(Unknown error code)" );
  }

  errorText += '\n' + tr( "The WCS vendor also reported: " );
  errorText += seText;

  QgsMessageLog::logMessage( tr( "composed error message '%1'." ).arg( errorText ), tr( "WCS" ) );
  QgsDebugMsg( QStringLiteral( "exiting." ) );
}



QgsRectangle QgsWcsProvider::extent() const
{
  if ( mExtentDirty )
  {
    if ( calculateExtent() )
    {
      mExtentDirty = false;
    }
  }

  return mCoverageExtent;
}

bool QgsWcsProvider::isValid() const
{
  return mValid;
}


QString QgsWcsProvider::wcsVersion()
{
  return mCapabilities.version();
}

bool QgsWcsProvider::calculateExtent() const
{

  // Make sure we know what extents are available
  if ( !mCoverageSummary.described )
  {
    return false;
  }

  // Prefer to use extent from capabilities / coverage description because
  // transformation from WGS84 enlarges the extent
  mCoverageExtent = mCoverageSummary.boundingBoxes.value( mCoverageCrs );
  QgsDebugMsg( "mCoverageCrs = " + mCoverageCrs + " mCoverageExtent = " + mCoverageExtent.toString() );
  if ( !mCoverageExtent.isEmpty() && mCoverageExtent.isFinite() )
  {
    QgsDebugMsg( "mCoverageExtent = " + mCoverageExtent.toString() );
    //return true;
  }
  else
  {
    // Set up the coordinate transform from the WCS standard CRS:84 bounding
    // box to the user's selected CRS
    if ( !mCoordinateTransform.isValid() )
    {
      QgsCoordinateReferenceSystem qgisSrsSource = QgsCoordinateReferenceSystem::fromOgcWmsCrs( QStringLiteral( "EPSG:4326" ) );
      QgsCoordinateReferenceSystem qgisSrsDest = QgsCoordinateReferenceSystem::fromOgcWmsCrs( mCoverageCrs );

      //QgsDebugMsg( "qgisSrsSource: " + qgisSrsSource.toWkt() );
      //QgsDebugMsg( "qgisSrsDest: " + qgisSrsDest.toWkt() );

      mCoordinateTransform = QgsCoordinateTransform( qgisSrsSource, qgisSrsDest, transformContext() );
    }

    QgsDebugMsg( "mCoverageSummary.wgs84BoundingBox= " + mCoverageSummary.wgs84BoundingBox.toString() );

    // Convert to the user's CRS as required
    try
    {
      mCoverageExtent = mCoordinateTransform.transformBoundingBox( mCoverageSummary.wgs84BoundingBox, QgsCoordinateTransform::ForwardTransform );
    }
    catch ( QgsCsException &cse )
    {
      Q_UNUSED( cse );
      return false;
    }

    //make sure extent does not contain 'inf' or 'nan'
    if ( !mCoverageExtent.isFinite() )
    {
      return false;
    }
  }

  QgsDebugMsg( "mCoverageExtent = " + mCoverageExtent.toString() );

  // It may happen (GeoServer) that extent reported in spatialDomain.Envelope is larger
  // than the coverage. Then if that larger BBOX is requested, the server returns
  // request BBOX intersected with coverage box scaled to requested WIDTH and HEIGHT.
  // GDAL WCS client does not suffer from this probably because it probably takes
  // extent from lonLatEnvelope (it probably does not calculate it from
  // spatialDomain.RectifiedGrid because calculated value is slightly different).

  // To get the true extent (it can also be smaller than real if reported Envelope is
  // than real smaller, but smaller is safer because data cannot be shifted) we make
  // request of the whole extent cut the extent from spatialDomain.Envelope if
  // necessary

  getCache( 1, mCoverageExtent, 10, 10 );
  if ( mCachedGdalDataset )
  {
    QgsRectangle cacheExtent = QgsGdalProviderBase::extent( mCachedGdalDataset.get() );
    QgsDebugMsg( "mCoverageExtent = " + mCoverageExtent.toString() );
    QgsDebugMsg( "cacheExtent = " + cacheExtent.toString() );
    QgsCoordinateReferenceSystem cacheCrs;
    if ( !cacheCrs.createFromWkt( GDALGetProjectionRef( mCachedGdalDataset.get() ) ) &&
         !cacheCrs.createFromWkt( GDALGetGCPProjection( mCachedGdalDataset.get() ) ) )
    {
      QgsDebugMsg( QStringLiteral( "Cached does not have CRS" ) );
    }
    QgsDebugMsg( "Cache CRS: " + cacheCrs.authid() + ' ' + cacheCrs.description() );

    // We can only verify extent if CRS is set
    // If dataset comes rotated, GDAL probably cuts latitude extend, disable
    // extent check for rotated, TODO: verify
    if ( cacheCrs.isValid() && !mFixRotate )
    {
      if ( !qgsDoubleNearSig( cacheExtent.xMinimum(), mCoverageExtent.xMinimum(), 10 ) ||
           !qgsDoubleNearSig( cacheExtent.yMinimum(), mCoverageExtent.yMinimum(), 10 ) ||
           !qgsDoubleNearSig( cacheExtent.xMaximum(), mCoverageExtent.xMaximum(), 10 ) ||
           !qgsDoubleNearSig( cacheExtent.yMaximum(), mCoverageExtent.yMaximum(), 10 ) )
      {
        QgsDebugMsg( QStringLiteral( "cacheExtent and mCoverageExtent differ, mCoverageExtent cut to cacheExtent" ) );
        mCoverageExtent = cacheExtent;
      }
    }
  }
  else
  {
    // Unfortunately it may also happen that a server (cubewerx.com) does not have
    // overviews and it is not able to respond for the whole extent within timeout.
    // It returns timeout error.
    // In that case (if request failed) we do not report error to allow working
    // with such servers on smaller portions of extent
    // (http://lists.osgeo.org/pipermail/qgis-developer/2013-January/024019.html)

    // Unfortunately even if we get over this 10x10 check, QGIS also requests
    // 32x32 thumbnail where it is waiting for another timeout

    QgsDebugMsg( QStringLiteral( "Cannot get cache to verify extent" ) );
    QgsMessageLog::logMessage( tr( "Cannot verify coverage full extent: %1" ).arg( mCachedError.message() ), tr( "WCS" ) );
  }

  return true;
}


int QgsWcsProvider::capabilities() const
{
  int capability = NoCapabilities;
  capability |= QgsRasterDataProvider::Identify;
  capability |= QgsRasterDataProvider::IdentifyValue;

  if ( mHasSize )
  {
    capability |= QgsRasterDataProvider::Size;
  }

  return capability;
}

QString QgsWcsProvider::coverageMetadata( const QgsWcsCoverageSummary &coverage )
{
  QString metadata;

  // Use a nested table
  metadata += QLatin1String( "<tr><td>" );
  metadata += QLatin1String( "<table width=\"100%\">" );

  // Table header
  metadata += QLatin1String( "<tr><th class=\"strong\">" );
  metadata += tr( "Property" );
  metadata += QLatin1String( "</th>" );
  metadata += QLatin1String( "<th class=\"strong\">" );
  metadata += tr( "Value" );
  metadata += QLatin1String( "</th></tr>" );

  metadata += htmlRow( tr( "Name (identifier)" ), coverage.identifier );
  metadata += htmlRow( tr( "Title" ), coverage.title );
  metadata += htmlRow( tr( "Abstract" ), coverage.abstract );
#if 0
  // We don't have size, nativeCrs, nativeBoundingBox etc. until describe coverage which would be heavy for all coverages
  metadata += htmlRow( tr( "Fixed Width" ), QString::number( coverage.width ) );
  metadata += htmlRow( tr( "Fixed Height" ), QString::number( coverage.height ) );
  metadata += htmlRow( tr( "Native CRS" ), coverage.nativeCrs );
  metadata += htmlRow( tr( "Native Bounding Box" ), coverage.nativeBoundingBox.toString() );
#endif

  metadata += htmlRow( tr( "WGS 84 Bounding Box" ), coverage.wgs84BoundingBox.toString() );

  // Layer Coordinate Reference Systems
  // TODO(?): supportedCrs and supportedFormat are not available in 1.0
  // until coverage is described - it would be confusing to show it only if available
#if 0
  for ( int j = 0; j < std::min( coverage.supportedCrs.size(), 10 ); j++ )
  {
    metadata += htmlRow( tr( "Available in CRS" ), coverage.supportedCrs.value( j ) );
  }

  if ( coverage.supportedCrs.size() > 10 )
  {
    metadata += htmlRow( tr( "Available in CRS" ), tr( "(and %n more)", "crs", coverage.supportedCrs.size() - 10 ) );
  }

  for ( int j = 0; j < std::min( coverage.supportedFormat.size(), 10 ); j++ )
  {
    metadata += htmlRow( tr( "Available in format" ), coverage.supportedFormat.value( j ) );
  }

  if ( coverage.supportedFormat.size() > 10 )
  {
    metadata += htmlRow( tr( "Available in format" ), tr( "(and %n more)", "crs", coverage.supportedFormat.size() - 10 ) );
  }
#endif

  // Close the nested table
  metadata += QLatin1String( "</table>" );
  metadata += QLatin1String( "</td></tr>" );

  return metadata;
}

QString QgsWcsProvider::htmlMetadata()
{
  QString metadata;
  metadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "WCS Info" ) + QStringLiteral( "</td><td><div>" );

  metadata += QLatin1String( "</a>&nbsp;<a href=\"#coverages\">" );
  metadata += tr( "Coverages" );
  metadata += QLatin1String( "</a>" );

#if 0
  // TODO
  metadata += "<a href=\"#cachestats\">";
  metadata += tr( "Cache Stats" );
  metadata += "</a> ";
#endif

  metadata += QLatin1String( "<br /><table class=\"tabular-view\">" );  // Nested table 1

  // Server Properties section
  metadata += QLatin1String( "<tr><th class=\"strong\"><a name=\"serverproperties\"></a>" );
  metadata += tr( "Server Properties" );
  metadata += QLatin1String( "</th></tr>" );

  // Use a nested table
  metadata += QLatin1String( "<tr><td>" );
  metadata += QLatin1String( "<table width=\"100%\">" );  // Nested table 2

  // Table header
  metadata += QLatin1String( "<tr><th class=\"strong\">" );
  metadata += tr( "Property" );
  metadata += QLatin1String( "</th>" );
  metadata += QLatin1String( "<th class=\"strong\">" );
  metadata += tr( "Value" );
  metadata += QLatin1String( "</th></tr>" );

  metadata += htmlRow( ( "WCS Version" ), mCapabilities.version() );
  metadata += htmlRow( tr( "Title" ), mCapabilities.capabilities().title );
  metadata += htmlRow( tr( "Abstract" ), mCapabilities.capabilities().abstract );
#if 0
  // TODO: probably apply stylesheet in QgsWcsCapabilities and save as html
  metadata += htmlRow( tr( "Keywords" ), mCapabilities.service.keywordList.join( "<br />" ) );
  metadata += htmlRow( tr( "Online Resource" ), "-" );
  metadata += htmlRow( tr( "Contact Person" ),
                       mCapabilities.service.contactInformation.contactPersonPrimary.contactPerson
                       + "<br />" + mCapabilities.service.contactInformation.contactPosition;
                       + "<br />" + mCapabilities.service.contactInformation.contactPersonPrimary.contactOrganization );
  metadata += htmlRow( tr( "Fees" ), mCapabilities.service.fees );
  metadata += htmlRow( tr( "Access Constraints" ), mCapabilities.service.accessConstraints );
  metadata += htmlRow( tr( "Image Formats" ), mCapabilities.capability.request.getMap.format.join( "<br />" ) );
  metadata += htmlRow( tr( "GetCapabilitiesUrl" ), mBaseUrl );
#endif
  metadata += htmlRow( tr( "Get Coverage Url" ), mCapabilities.getCoverageUrl() + ( mIgnoreGetCoverageUrl ? tr( "&nbsp;<font color=\"red\">(advertised but ignored)</font>" ) : QString() ) );

  // Close the nested table
  metadata += QLatin1String( "</table>" );  // End nested table 2
  metadata += QLatin1String( "</td></tr>" );

  // Coverage properties
  metadata += QLatin1String( "<tr><th class=\"strong\"><a name=\"coverages\"></a>" );
  metadata += tr( "Coverages" );
  metadata += QLatin1String( "</th></tr>" );

  // Dialog takes too long to open if there are too many coverages (1000 for example)
  int count = 0;
  const auto constCoverages = mCapabilities.coverages();
  for ( const QgsWcsCoverageSummary &c : constCoverages )
  {
    metadata += coverageMetadata( c );
    count++;
    if ( count >= 100 ) break;
  }
  metadata += QLatin1String( "</table>" );
  if ( count < mCapabilities.coverages().size() )
  {
    metadata += tr( "And %1 more coverages" ).arg( mCapabilities.coverages().size() - count );
  }

  metadata += QStringLiteral( "</table></div></td></tr>\n" );  // End nested table 1
  return metadata;
}

QString QgsWcsProvider::htmlCell( const QString &text )
{
  return "<td>" + text + "</td>";
}

QString QgsWcsProvider:: htmlRow( const QString &text1, const QString &text2 )
{
  return "<tr>" + htmlCell( text1 ) +  htmlCell( text2 ) + "</tr>";
}

QgsRasterIdentifyResult QgsWcsProvider::identify( const QgsPointXY &point, QgsRaster::IdentifyFormat format, const QgsRectangle &boundingBox, int width, int height, int /*dpi*/ )
{
  QgsDebugMsg( QStringLiteral( "thePoint = %1 %2" ).arg( point.x(), 0, 'g', 10 ).arg( point.y(), 0, 'g', 10 ) );
  QgsDebugMsg( QStringLiteral( "theWidth = %1 height = %2" ).arg( width ).arg( height ) );
  QgsDebugMsg( "boundingBox = " + boundingBox.toString() );
  QMap<int, QVariant> results;

  if ( format != QgsRaster::IdentifyFormatValue )
  {
    return QgsRasterIdentifyResult( QGS_ERROR( tr( "Format not supported" ) ) );
  }

  if ( !extent().contains( point ) )
  {
    // Outside the raster
    for ( int i = 1; i <= bandCount(); i++ )
    {
      results.insert( i, QVariant() );
    }
    return QgsRasterIdentifyResult( QgsRaster::IdentifyFormatValue, results );
  }

  QgsRectangle finalExtent = boundingBox;
  int maxSize = 2000;
  int cacheSize = 1000; // tile cache size if context is not defined or small
  double relResTol = 0.1; // relative resolution tolerance (10%)

  // TODO: We are using cacheSize x cacheSize if context is not defined
  // or box is too small (in pixels). That is necessary to allow effective
  // map querying (valuetool for example). OTOH, it means reading of large
  // unnecessary amount of data for each identify() call if query points are too
  // distant (out of cache) for example when called from scripts

  // if context size is to large we have to cut it, in that case caching big
  // big part does not make sense
  if ( finalExtent.isEmpty() || width == 0 || height == 0 ||
       width > maxSize || height > maxSize )
  {
    // context missing, use a small area around the point and highest resolution if known

    double xRes, yRes;
    if ( mHasSize )
    {
      xRes = mCoverageExtent.width()  / mWidth;
      yRes = mCoverageExtent.height()  / mHeight;
    }
    else
    {
      // set resolution approximately to 1mm
      switch ( mCrs.mapUnits() )
      {
        case QgsUnitTypes::DistanceMeters:
          xRes = 0.001;
          break;
        case QgsUnitTypes::DistanceFeet:
          xRes = 0.003;
          break;
        case QgsUnitTypes::DistanceDegrees:
          // max length of degree of latitude on pole is 111694 m
          xRes = 1e-8;
          break;
        default:
          xRes = 0.001; // expecting meters
      }
      yRes = xRes;
    }

    width = cacheSize;
    height = cacheSize;

    finalExtent = QgsRectangle( point.x() - xRes * width / 2,
                                point.y() - yRes * height / 2,
                                point.x() + xRes * width / 2,
                                point.y() + yRes * height / 2 );

    double xResDiff = std::fabs( mCachedViewExtent.width() / mCachedViewWidth - xRes );
    double yResDiff = std::fabs( mCachedViewExtent.height() / mCachedViewHeight - yRes );

    if ( !mCachedGdalDataset ||
         !mCachedViewExtent.contains( point ) ||
         mCachedViewWidth == 0 || mCachedViewHeight == 0 ||
         xResDiff / xRes > relResTol ||
         yResDiff / yRes > relResTol )
    {
      getCache( 1, finalExtent, width, height );
    }
  }
  else
  {
    // Use context -> effective caching (usually, if context is constant)
    QgsDebugMsg( QStringLiteral( "Using context extent and resolution" ) );
    // To use the cache it is sufficient to have point within cache and
    // similar resolution
    double xRes = finalExtent.width() / width;
    double yRes = finalExtent.height() / height;
    QgsDebugMsg( QStringLiteral( "width = %1 height = %2 xRes = %3 yRes = %4" ).arg( finalExtent.width() ).arg( finalExtent.height() ).arg( xRes ).arg( yRes ) );
    double xResDiff = std::fabs( mCachedViewExtent.width() / mCachedViewWidth - xRes );
    double yResDiff = std::fabs( mCachedViewExtent.height() / mCachedViewHeight - yRes );
    QgsDebugMsg( QStringLiteral( "xRes diff = %1 yRes diff = %2 relative xResDiff = %3 relative yResDiff = %4" ).arg( xResDiff ).arg( yResDiff ).arg( xResDiff / xRes ).arg( yResDiff / yRes ) );
    if ( !mCachedGdalDataset ||
         !mCachedViewExtent.contains( point ) ||
         xResDiff / xRes > relResTol ||
         yResDiff / yRes > relResTol )
    {
      // Identify map tool is now using fake context 1x1 pixel to get point/line
      // features from WMS. In such case we enlarge the extent to get data cached
      // for next queries around.
      // BTW: UMN Mapserver (6.0.3) seems to be buggy with 1x1 pixels request (returns 'no data' value
      if ( width < cacheSize )
      {
        int buffer = ( cacheSize - width ) / 2;
        width += 2 * buffer;
        finalExtent.setXMinimum( finalExtent.xMinimum() - xRes * buffer );
        finalExtent.setXMaximum( finalExtent.xMaximum() + xRes * buffer );
      }
      if ( height < cacheSize )
      {
        int buffer = ( cacheSize - height ) / 2;
        height += 2 * buffer;
        finalExtent.setYMinimum( finalExtent.yMinimum() - yRes * buffer );
        finalExtent.setYMaximum( finalExtent.yMaximum() + yRes * buffer );
      }

      getCache( 1, finalExtent, width, height );
    }
  }

  if ( !mCachedGdalDataset ||
       !mCachedViewExtent.contains( point ) )
  {
    return QgsRasterIdentifyResult( QGS_ERROR( tr( "Read data error" ) ) );
  }

  double x = point.x();
  double y = point.y();

  // Calculate the row / column where the point falls
  double xRes = mCachedViewExtent.width() / mCachedViewWidth;
  double yRes = mCachedViewExtent.height() / mCachedViewHeight;

  // Offset, not the cell index -> flor
  int col = ( int ) std::floor( ( x - mCachedViewExtent.xMinimum() ) / xRes );
  int row = ( int ) std::floor( ( mCachedViewExtent.yMaximum() - y ) / yRes );

  QgsDebugMsg( "row = " + QString::number( row ) + " col = " + QString::number( col ) );

  for ( int i = 1; i <= GDALGetRasterCount( mCachedGdalDataset.get() ); i++ )
  {
    GDALRasterBandH gdalBand = GDALGetRasterBand( mCachedGdalDataset.get(), i );

    double value;
    CPLErr err = GDALRasterIO( gdalBand, GF_Read, col, row, 1, 1,
                               &value, 1, 1, GDT_Float64, 0, 0 );

    if ( err != CPLE_None )
    {
      QgsLogger::warning( "RasterIO error: " + QString::fromUtf8( CPLGetLastErrorMsg() ) );
      return QgsRasterIdentifyResult( QGS_ERROR( tr( "RasterIO error: " ) + QString::fromUtf8( CPLGetLastErrorMsg() ) ) );
    }

    // Apply no data and user no data
    if ( ( sourceHasNoDataValue( i ) && useSourceNoDataValue( i ) &&
           ( std::isnan( value ) || qgsDoubleNear( value, sourceNoDataValue( i ) ) ) ) ||
         ( QgsRasterRange::contains( value, userNoDataValues( i ) ) ) )
    {
      results.insert( i, QVariant() );
    }
    else
    {
      results.insert( i, value );
    }
  }

  return QgsRasterIdentifyResult( QgsRaster::IdentifyFormatValue, results );
}

QgsCoordinateReferenceSystem QgsWcsProvider::crs() const
{
  return mCrs;
}

QString QgsWcsProvider::lastErrorTitle()
{
  return mErrorCaption;
}

QString QgsWcsProvider::lastError()
{
  QgsDebugMsg( "returning '" + mError  + "'." );
  return mError;
}

QString QgsWcsProvider::lastErrorFormat()
{
  return mErrorFormat;
}

QString  QgsWcsProvider::name() const
{
  return WCS_KEY;
}

QString  QgsWcsProvider::description() const
{
  return WCS_DESCRIPTION;
}

void QgsWcsProvider::reloadData()
{
  clearCache();
}

QString QgsWcsProvider::nodeAttribute( const QDomElement &e, const QString &name, const QString &defValue )
{
  if ( e.hasAttribute( name ) )
    return e.attribute( name );

  QDomNamedNodeMap map( e.attributes() );
  for ( int i = 0; i < map.size(); i++ )
  {
    QDomAttr attr( map.item( i ).toElement().toAttr() );
    if ( attr.name().compare( name, Qt::CaseInsensitive ) == 0 )
      return attr.value();
  }

  return defValue;
}

QMap<QString, QString> QgsWcsProvider::supportedMimes()
{
  QMap<QString, QString> mimes;
  GDALAllRegister();

  QgsDebugMsg( QStringLiteral( "GDAL drivers cont %1" ).arg( GDALGetDriverCount() ) );
  for ( int i = 0; i < GDALGetDriverCount(); ++i )
  {
    GDALDriverH driver = GDALGetDriver( i );
    Q_CHECK_PTR( driver );

    if ( !driver )
    {
      QgsLogger::warning( "unable to get driver " + QString::number( i ) );
      continue;
    }

    QString desc = GDALGetDescription( driver );

    QString mimeType = GDALGetMetadataItem( driver, "DMD_MIMETYPE", "" );

    if ( mimeType.isEmpty() ) continue;

    desc = desc.isEmpty() ? mimeType : desc;

    QgsDebugMsg( "add GDAL format " + mimeType + ' ' + desc );

    mimes[mimeType] = desc;
  }
  return mimes;
}

QGISEXTERN QgsWcsProvider *classFactory( const QString *uri, const QgsDataProvider::ProviderOptions &options )
{
  return new QgsWcsProvider( *uri, options );
}

QGISEXTERN QString providerKey()
{
  return WCS_KEY;
}

QGISEXTERN QString description()
{
  return WCS_DESCRIPTION;
}

QGISEXTERN bool isProvider()
{
  return true;
}


// ----------

int QgsWcsDownloadHandler::sErrors = 0;

QgsWcsDownloadHandler::QgsWcsDownloadHandler( const QUrl &url, QgsWcsAuthorization &auth, QNetworkRequest::CacheLoadControl cacheLoadControl, QByteArray &cachedData, const QString &wcsVersion, QgsError &cachedError, QgsRasterBlockFeedback *feedback )
  : mAuth( auth )
  , mEventLoop( new QEventLoop )
  , mCachedData( cachedData )
  , mWcsVersion( wcsVersion )
  , mCachedError( cachedError )
  , mFeedback( feedback )
{
  if ( feedback )
  {
    connect( feedback, &QgsFeedback::canceled, this, &QgsWcsDownloadHandler::canceled, Qt::QueuedConnection );

    // rendering could have been canceled before we started to listen to canceled() signal
    // so let's check before doing the download and maybe quit prematurely
    if ( feedback->isCanceled() )
      return;
  }

  QNetworkRequest request( url );
  QgsSetRequestInitiatorClass( request, QStringLiteral( "QgsWcsDownloadHandler" ) );
  if ( !mAuth.setAuthorization( request ) )
  {
    QgsMessageLog::logMessage( tr( "Network request update failed for authentication config" ),
                               tr( "WCS" ) );
    return;
  }
  request.setAttribute( QNetworkRequest::CacheSaveControlAttribute, true );
  request.setAttribute( QNetworkRequest::CacheLoadControlAttribute, cacheLoadControl );

  mCacheReply = QgsNetworkAccessManager::instance()->get( request );
  if ( !mAuth.setAuthorizationReply( mCacheReply ) )
  {
    mCacheReply->deleteLater();
    mCacheReply = nullptr;
    QgsMessageLog::logMessage( tr( "Network reply update failed for authentication config" ),
                               tr( "WCS" ) );
    finish();
    return;
  }
  connect( mCacheReply, &QNetworkReply::finished, this, &QgsWcsDownloadHandler::cacheReplyFinished );
  connect( mCacheReply, &QNetworkReply::downloadProgress, this, &QgsWcsDownloadHandler::cacheReplyProgress );
}

QgsWcsDownloadHandler::~QgsWcsDownloadHandler()
{
  delete mEventLoop;
}

void QgsWcsDownloadHandler::blockingDownload()
{
  if ( mFeedback && mFeedback->isCanceled() )
    return; // nothing to do

  mEventLoop->exec( QEventLoop::ExcludeUserInputEvents );

  Q_ASSERT( !mCacheReply );
}

void QgsWcsDownloadHandler::cacheReplyFinished()
{
  QgsDebugMsg( QStringLiteral( "mCacheReply->error() = %1" ).arg( mCacheReply->error() ) );
  if ( mCacheReply->error() == QNetworkReply::NoError )
  {
    QVariant redirect = mCacheReply->attribute( QNetworkRequest::RedirectionTargetAttribute );
    if ( !redirect.isNull() )
    {
      mCacheReply->deleteLater();

      QgsDebugMsg( QStringLiteral( "redirected getmap: %1" ).arg( redirect.toString() ) );
      QNetworkRequest request( redirect.toUrl() );
      QgsSetRequestInitiatorClass( request, QStringLiteral( "QgsWcsDownloadHandler" ) );
      if ( !mAuth.setAuthorization( request ) )
      {
        QgsMessageLog::logMessage( tr( "Network request update failed for authentication config" ),
                                   tr( "WCS" ) );
        return;
      }
      mCacheReply = QgsNetworkAccessManager::instance()->get( request );
      if ( !mAuth.setAuthorizationReply( mCacheReply ) )
      {
        mCacheReply->deleteLater();
        mCacheReply = nullptr;
        QgsMessageLog::logMessage( tr( "Network reply update failed for authentication config" ),
                                   tr( "WCS" ) );
        finish();
        return;
      }
      connect( mCacheReply, &QNetworkReply::finished, this, &QgsWcsDownloadHandler::cacheReplyFinished );
      connect( mCacheReply, &QNetworkReply::downloadProgress, this, &QgsWcsDownloadHandler::cacheReplyProgress );

      return;
    }

    QVariant status = mCacheReply->attribute( QNetworkRequest::HttpStatusCodeAttribute );
    QgsDebugMsg( QStringLiteral( "status = %1" ).arg( status.toInt() ) );
    if ( !status.isNull() && status.toInt() >= 400 )
    {
      QVariant phrase = mCacheReply->attribute( QNetworkRequest::HttpReasonPhraseAttribute );

      QgsMessageLog::logMessage( tr( "Map request error (Status: %1; Reason phrase: %2; URL: %3)" )
                                 .arg( status.toInt() )
                                 .arg( phrase.toString(),
                                       mCacheReply->url().toString() ), tr( "WCS" ) );

      mCacheReply->deleteLater();
      mCacheReply = nullptr;

      finish();
      return;
    }

    // Read response

    QString contentType = mCacheReply->header( QNetworkRequest::ContentTypeHeader ).toString();
    QgsDebugMsg( "contentType: " + contentType );

    // Exception
    // Content type examples: text/xml
    //                        application/vnd.ogc.se_xml;charset=UTF-8
    //                        application/xml
    if ( contentType.startsWith( QLatin1String( "text/" ), Qt::CaseInsensitive ) ||
         contentType.compare( QLatin1String( "application/xml" ), Qt::CaseInsensitive ) == 0 ||
         contentType.startsWith( QLatin1String( "application/vnd.ogc.se_xml" ), Qt::CaseInsensitive ) )
    {
      QString errorTitle, errorText;
      QByteArray text = mCacheReply->readAll();
      if ( ( contentType.compare( QLatin1String( "text/xml" ), Qt::CaseInsensitive ) == 0 ||
             contentType.compare( QLatin1String( "application/xml" ), Qt::CaseInsensitive ) == 0 ||
             contentType.startsWith( QLatin1String( "application/vnd.ogc.se_xml" ), Qt::CaseInsensitive ) )
           && QgsWcsProvider::parseServiceExceptionReportDom( text, mWcsVersion, errorTitle, errorText ) )
      {
        mCachedError.append( SRVERR( tr( "Map request error:<br>Title: %1<br>Error: %2<br>URL: <a href='%3'>%3</a>)" )
                                     .arg( errorTitle, errorText,
                                           mCacheReply->url().toString() ) ) );
      }
      else
      {
        QgsMessageLog::logMessage( tr( "Map request error (Status: %1; Response: %2; URL: %3)" )
                                   .arg( status.toInt() )
                                   .arg( QString::fromUtf8( text ),
                                         mCacheReply->url().toString() ), tr( "WCS" ) );
      }

      mCacheReply->deleteLater();
      mCacheReply = nullptr;

      finish();
      return;
    }

    // WCS 1.1 gives response as multipart, e.g.
    if ( QgsNetworkReplyParser::isMultipart( mCacheReply ) )
    {
      QgsDebugMsg( QStringLiteral( "reply is multipart" ) );
      QgsNetworkReplyParser parser( mCacheReply );

      if ( !parser.isValid() )
      {
        QgsMessageLog::logMessage( tr( "Cannot parse multipart response: %1" ).arg( parser.error() ), tr( "WCS" ) );
        mCacheReply->deleteLater();
        mCacheReply = nullptr;

        finish();
        return;
      }

      if ( parser.parts() < 2 )
      {
        QgsMessageLog::logMessage( tr( "Expected 2 parts, %1 received" ).arg( parser.parts() ), tr( "WCS" ) );
        mCacheReply->deleteLater();
        mCacheReply = nullptr;

        finish();
        return;
      }
      else if ( parser.parts() > 2 )
      {
        // We will try the second one
        QgsMessageLog::logMessage( tr( "More than 2 parts (%1) received" ).arg( parser.parts() ), tr( "WCS" ) );
      }

      QString transferEncoding = parser.rawHeader( 1, QStringLiteral( "Content-Transfer-Encoding" ).toLatin1() );
      QgsDebugMsg( "transferEncoding = " + transferEncoding );

      // It may happen (GeoServer) that in part header is for example
      // Content-Type: image/tiff and Content-Transfer-Encoding: base64
      // but content is xml ExceptionReport which is not in base64
      QByteArray body = parser.body( 1 );
      if ( body.startsWith( "<?xml" ) )
      {
        QString errorTitle, errorText;
        if ( QgsWcsProvider::parseServiceExceptionReportDom( body, mWcsVersion, errorTitle, errorText ) )
        {
          QgsMessageLog::logMessage( tr( "Map request error (Title: %1; Error: %2; URL: %3)" )
                                     .arg( errorTitle, errorText,
                                           mCacheReply->url().toString() ), tr( "WCS" ) );
        }
        else
        {
          QgsMessageLog::logMessage( tr( "Map request error (Response: %1; URL: %2)" )
                                     .arg( QString::fromUtf8( body ),
                                           mCacheReply->url().toString() ), tr( "WCS" ) );
        }

        mCacheReply->deleteLater();
        mCacheReply = nullptr;

        finish();
        return;
      }

      if ( transferEncoding == QLatin1String( "binary" ) )
      {
        mCachedData = body;
      }
      else if ( transferEncoding == QLatin1String( "base64" ) )
      {
        mCachedData = QByteArray::fromBase64( body );
      }
      else
      {
        QgsMessageLog::logMessage( tr( "Content-Transfer-Encoding %1 not supported" ).arg( transferEncoding ), tr( "WCS" ) );
      }
    }
    else
    {
      // Mime types for WCS 1.0 response should be usually image/<format>
      // (image/tiff, image/png, image/jpeg, image/png; mode=8bit, etc.)
      // but other mime types (like application/*) may probably also appear

      // Create memory file
      mCachedData = mCacheReply->readAll();
    }

    mCacheReply->deleteLater();
    mCacheReply = nullptr;

    finish();
  }
  else
  {
    // report any errors except for the one we have caused by canceling the request
    if ( mCacheReply->error() != QNetworkReply::OperationCanceledError )
    {
      // Resend request if AlwaysCache
      QNetworkRequest request = mCacheReply->request();
      QgsSetRequestInitiatorClass( request, QStringLiteral( "QgsWcsDownloadHandler" ) );
      if ( request.attribute( QNetworkRequest::CacheLoadControlAttribute ).toInt() == QNetworkRequest::AlwaysCache )
      {
        QgsDebugMsg( QStringLiteral( "Resend request with PreferCache" ) );
        request.setAttribute( QNetworkRequest::CacheLoadControlAttribute, QNetworkRequest::PreferCache );

        mCacheReply->deleteLater();

        mCacheReply = QgsNetworkAccessManager::instance()->get( request );
        if ( !mAuth.setAuthorizationReply( mCacheReply ) )
        {
          mCacheReply->deleteLater();
          mCacheReply = nullptr;
          QgsMessageLog::logMessage( tr( "Network reply update failed for authentication config" ),
                                     tr( "WCS" ) );
          finish();
          return;
        }
        connect( mCacheReply, &QNetworkReply::finished, this, &QgsWcsDownloadHandler::cacheReplyFinished, Qt::DirectConnection );
        connect( mCacheReply, &QNetworkReply::downloadProgress, this, &QgsWcsDownloadHandler::cacheReplyProgress, Qt::DirectConnection );

        return;
      }

      sErrors++;
      if ( sErrors < 100 )
      {
        QgsMessageLog::logMessage( tr( "Map request failed [error: %1 url: %2]" ).arg( mCacheReply->errorString(), mCacheReply->url().toString() ), tr( "WCS" ) );
      }
      else if ( sErrors == 100 )
      {
        QgsMessageLog::logMessage( tr( "Not logging more than 100 request errors." ), tr( "WCS" ) );
      }
    }

    mCacheReply->deleteLater();
    mCacheReply = nullptr;

    finish();
  }
}

void QgsWcsDownloadHandler::cacheReplyProgress( qint64 bytesReceived, qint64 bytesTotal )
{
  Q_UNUSED( bytesReceived );
  Q_UNUSED( bytesTotal );
  QgsDebugMsgLevel( QStringLiteral( "%1 of %2 bytes of map downloaded." ).arg( bytesReceived ).arg( bytesTotal < 0 ? QString( "unknown number of" ) : QString::number( bytesTotal ) ), 3 );
}

void QgsWcsDownloadHandler::canceled()
{
  QgsDebugMsg( QStringLiteral( "Caught canceled() signal" ) );
  if ( mCacheReply )
  {
    QgsDebugMsg( QStringLiteral( "Aborting WCS network request" ) );
    mCacheReply->abort();
  }
}

#ifdef HAVE_GUI

//! Provider for WCS layers source select
class QgsWcsSourceSelectProvider : public QgsSourceSelectProvider
{
  public:

    QString providerKey() const override { return QStringLiteral( "wcs" ); }
    QString text() const override { return QObject::tr( "WCS" ); }
    int ordering() const override { return QgsSourceSelectProvider::OrderRemoteProvider + 20; }
    QIcon icon() const override { return QgsApplication::getThemeIcon( QStringLiteral( "/mActionAddWcsLayer.svg" ) ); }
    QgsAbstractDataSourceWidget *createDataSourceWidget( QWidget *parent = nullptr, Qt::WindowFlags fl = Qt::Widget, QgsProviderRegistry::WidgetMode widgetMode = QgsProviderRegistry::WidgetMode::Embedded ) const override
    {
      return new QgsWCSSourceSelect( parent, fl, widgetMode );
    }
};


QGISEXTERN QList<QgsSourceSelectProvider *> *sourceSelectProviders()
{
  QList<QgsSourceSelectProvider *> *providers = new QList<QgsSourceSelectProvider *>();

  *providers
      << new QgsWcsSourceSelectProvider;

  return providers;
}
#endif
