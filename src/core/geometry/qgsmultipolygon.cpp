/***************************************************************************
                        qgsmultipolygon.cpp
  -------------------------------------------------------------------
Date                 : 28 Oct 2014
Copyright            : (C) 2014 by Marco Hugentobler
email                : marco.hugentobler at sourcepole dot com
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgsmultipolygon.h"
#include "qgsapplication.h"
#include "qgsgeometryutils.h"
#include "qgssurface.h"
#include "qgslinestring.h"
#include "qgspolygon.h"
#include "qgscurvepolygon.h"
#include "qgsmultilinestring.h"

QgsMultiPolygon::QgsMultiPolygon()
{
  mWkbType = QgsWkbTypes::MultiPolygon;
}

QString QgsMultiPolygon::geometryType() const
{
  return QStringLiteral( "MultiPolygon" );
}

void QgsMultiPolygon::clear()
{
  QgsMultiSurface::clear();
  mWkbType = QgsWkbTypes::MultiPolygon;
}

QgsMultiPolygon *QgsMultiPolygon::createEmptyWithSameType() const
{
  auto result = qgis::make_unique< QgsMultiPolygon >();
  result->mWkbType = mWkbType;
  return result.release();
}

QgsMultiPolygon *QgsMultiPolygon::clone() const
{
  return new QgsMultiPolygon( *this );
}

bool QgsMultiPolygon::fromWkt( const QString &wkt )
{
  return fromCollectionWkt( wkt, QVector<QgsAbstractGeometry *>() << new QgsPolygon, QStringLiteral( "Polygon" ) );
}

QDomElement QgsMultiPolygon::asGml2( QDomDocument &doc, int precision, const QString &ns, const AxisOrder axisOrder ) const
{
  // GML2 does not support curves
  QDomElement elemMultiPolygon = doc.createElementNS( ns, QStringLiteral( "MultiPolygon" ) );

  if ( isEmpty() )
    return elemMultiPolygon;

  for ( const QgsAbstractGeometry *geom : mGeometries )
  {
    if ( qgsgeometry_cast<const QgsPolygon *>( geom ) )
    {
      QDomElement elemPolygonMember = doc.createElementNS( ns, QStringLiteral( "polygonMember" ) );
      elemPolygonMember.appendChild( geom->asGml2( doc, precision, ns, axisOrder ) );
      elemMultiPolygon.appendChild( elemPolygonMember );
    }
  }

  return elemMultiPolygon;
}

QDomElement QgsMultiPolygon::asGml3( QDomDocument &doc, int precision, const QString &ns, const QgsAbstractGeometry::AxisOrder axisOrder ) const
{
  QDomElement elemMultiSurface = doc.createElementNS( ns, QStringLiteral( "MultiPolygon" ) );

  if ( isEmpty() )
    return elemMultiSurface;

  for ( const QgsAbstractGeometry *geom : mGeometries )
  {
    if ( qgsgeometry_cast<const QgsPolygon *>( geom ) )
    {
      QDomElement elemSurfaceMember = doc.createElementNS( ns, QStringLiteral( "polygonMember" ) );
      elemSurfaceMember.appendChild( geom->asGml3( doc, precision, ns, axisOrder ) );
      elemMultiSurface.appendChild( elemSurfaceMember );
    }
  }

  return elemMultiSurface;
}

QString QgsMultiPolygon::asJson( int precision ) const
{
  // GeoJSON does not support curves
  QString json = QStringLiteral( "{\"type\": \"MultiPolygon\", \"coordinates\": [" );
  for ( const QgsAbstractGeometry *geom : mGeometries )
  {
    if ( qgsgeometry_cast<const QgsPolygon *>( geom ) )
    {
      json += '[';

      const QgsPolygon *polygon = static_cast<const QgsPolygon *>( geom );

      std::unique_ptr< QgsLineString > exteriorLineString( polygon->exteriorRing()->curveToLine() );
      QgsPointSequence exteriorPts;
      exteriorLineString->points( exteriorPts );
      json += QgsGeometryUtils::pointsToJSON( exteriorPts, precision ) + QLatin1String( ", " );

      std::unique_ptr< QgsLineString > interiorLineString;
      for ( int i = 0, n = polygon->numInteriorRings(); i < n; ++i )
      {
        interiorLineString.reset( polygon->interiorRing( i )->curveToLine() );
        QgsPointSequence interiorPts;
        interiorLineString->points( interiorPts );
        json += QgsGeometryUtils::pointsToJSON( interiorPts, precision ) + QLatin1String( ", " );
      }
      if ( json.endsWith( QLatin1String( ", " ) ) )
      {
        json.chop( 2 ); // Remove last ", "
      }

      json += QLatin1String( "], " );
    }
  }
  if ( json.endsWith( QLatin1String( ", " ) ) )
  {
    json.chop( 2 ); // Remove last ", "
  }
  json += QLatin1String( "] }" );
  return json;
}

bool QgsMultiPolygon::addGeometry( QgsAbstractGeometry *g )
{
  if ( !qgsgeometry_cast<QgsPolygon *>( g ) )
  {
    delete g;
    return false;
  }

  if ( mGeometries.empty() )
  {
    setZMTypeFromSubGeometry( g, QgsWkbTypes::MultiPolygon );
  }
  if ( is3D() && !g->is3D() )
    g->addZValue();
  else if ( !is3D() && g->is3D() )
    g->dropZValue();
  if ( isMeasure() && !g->isMeasure() )
    g->addMValue();
  else if ( !isMeasure() && g->isMeasure() )
    g->dropMValue();

  return QgsGeometryCollection::addGeometry( g ); // clazy:exclude=skipped-base-method
}

bool QgsMultiPolygon::insertGeometry( QgsAbstractGeometry *g, int index )
{
  if ( !g || !qgsgeometry_cast< QgsPolygon * >( g ) )
  {
    delete g;
    return false;
  }

  return QgsMultiSurface::insertGeometry( g, index );
}

QgsMultiSurface *QgsMultiPolygon::toCurveType() const
{
  QgsMultiSurface *multiSurface = new QgsMultiSurface();
  for ( int i = 0; i < mGeometries.size(); ++i )
  {
    multiSurface->addGeometry( mGeometries.at( i )->clone() );
  }
  return multiSurface;
}

QgsAbstractGeometry *QgsMultiPolygon::boundary() const
{
  std::unique_ptr< QgsMultiLineString > multiLine( new QgsMultiLineString() );
  for ( int i = 0; i < mGeometries.size(); ++i )
  {
    if ( QgsPolygon *polygon = qgsgeometry_cast<QgsPolygon *>( mGeometries.at( i ) ) )
    {
      QgsAbstractGeometry *polygonBoundary = polygon->boundary();

      if ( QgsLineString *lineStringBoundary = qgsgeometry_cast< QgsLineString * >( polygonBoundary ) )
      {
        multiLine->addGeometry( lineStringBoundary );
      }
      else if ( QgsMultiLineString *multiLineStringBoundary = qgsgeometry_cast< QgsMultiLineString * >( polygonBoundary ) )
      {
        for ( int j = 0; j < multiLineStringBoundary->numGeometries(); ++j )
        {
          multiLine->addGeometry( multiLineStringBoundary->geometryN( j )->clone() );
        }
        delete multiLineStringBoundary;
      }
      else
      {
        delete polygonBoundary;
      }
    }
  }
  if ( multiLine->numGeometries() == 0 )
  {
    return nullptr;
  }
  return multiLine.release();
}

bool QgsMultiPolygon::wktOmitChildType() const
{
  return true;
}
