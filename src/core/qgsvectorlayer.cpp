/***************************************************************************
                               qgsvectorlayer.cpp
                              --------------------
          begin                : Oct 29, 2003
          copyright            : (C) 2003 by Gary E.Sherman
          email                : sherman at mrcc.com

  This class implements a generic means to display vector layers. The features
  and attributes are read from the data store using a "data provider" plugin.
  QgsVectorLayer can be used with any data store for which an appropriate
  plugin is available.

***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <limits>

#include <QDir>
#include <QFile>
#include <QImage>
#include <QPainter>
#include <QPainterPath>
#include <QPolygonF>
#include <QProgressDialog>
#include <QString>
#include <QDomNode>
#include <QVector>
#include <QStringBuilder>
#include <QUrl>
#include <QUndoCommand>

#include "qgssettings.h"
#include "qgsvectorlayer.h"
#include "qgsactionmanager.h"
#include "qgis.h" //for globals
#include "qgsapplication.h"
#include "qgsclipper.h"
#include "qgsconditionalstyle.h"
#include "qgscoordinatereferencesystem.h"
#include "qgscoordinatetransform.h"
#include "qgsexception.h"
#include "qgscurve.h"
#include "qgsdatasourceuri.h"
#include "qgsexpressionfieldbuffer.h"
#include "qgsexpressionnodeimpl.h"
#include "qgsfeature.h"
#include "qgsfeaturerequest.h"
#include "qgsfields.h"
#include "qgsgeometry.h"
#include "qgslayermetadataformatter.h"
#include "qgslogger.h"
#include "qgsmaplayerlegend.h"
#include "qgsmaptopixel.h"
#include "qgsmessagelog.h"
#include "qgsogcutils.h"
#include "qgspainting.h"
#include "qgspointxy.h"
#include "qgsproject.h"
#include "qgsproviderregistry.h"
#include "qgsrectangle.h"
#include "qgsrelationmanager.h"
#include "qgsrendercontext.h"
#include "qgsvectordataprovider.h"
#include "qgsvectorlayereditbuffer.h"
#include "qgsvectorlayereditpassthrough.h"
#include "qgsvectorlayereditutils.h"
#include "qgsvectorlayerfeatureiterator.h"
#include "qgsvectorlayerjoinbuffer.h"
#include "qgsvectorlayerlabeling.h"
#include "qgsvectorlayerrenderer.h"
#include "qgsvectorlayerundocommand.h"
#include "qgsvectorlayerfeaturecounter.h"
#include "qgspoint.h"
#include "qgsrenderer.h"
#include "qgssymbollayer.h"
#include "qgssinglesymbolrenderer.h"
#include "qgsdiagramrenderer.h"
#include "qgsstyle.h"
#include "qgspallabeling.h"
#include "qgssimplifymethod.h"
#include "qgsexpressioncontext.h"
#include "qgsfeedback.h"
#include "qgsxmlutils.h"
#include "qgsunittypes.h"
#include "qgstaskmanager.h"
#include "qgstransaction.h"
#include "qgsauxiliarystorage.h"
#include "qgsgeometryoptions.h"
#include "qgsexpressioncontextutils.h"

#include "diagram/qgsdiagram.h"

#ifdef TESTPROVIDERLIB
#include <dlfcn.h>
#endif

typedef bool saveStyle_t(
  const QString &uri,
  const QString &qmlStyle,
  const QString &sldStyle,
  const QString &styleName,
  const QString &styleDescription,
  const QString &uiFileContent,
  bool useAsDefault,
  QString &errCause
);

typedef QString loadStyle_t(
  const QString &uri,
  QString &errCause
);

typedef int listStyles_t(
  const QString &uri,
  QStringList &ids,
  QStringList &names,
  QStringList &descriptions,
  QString &errCause
);

typedef QString getStyleById_t(
  const QString &uri,
  QString styleID,
  QString &errCause
);

typedef bool deleteStyleById_t(
  const QString &uri,
  QString styleID,
  QString &errCause
);


QgsVectorLayer::QgsVectorLayer( const QString &vectorLayerPath,
                                const QString &baseName,
                                const QString &providerKey,
                                const QgsVectorLayer::LayerOptions &options )
  : QgsMapLayer( QgsMapLayerType::VectorLayer, baseName, vectorLayerPath )
  , mAuxiliaryLayer( nullptr )
  , mAuxiliaryLayerKey( QString() )
  , mReadExtentFromXml( options.readExtentFromXml )
{
  if ( options.fallbackCrs.isValid() )
    setCrs( options.fallbackCrs, false );
  mWkbType = options.fallbackWkbType;

  setProviderType( providerKey );

  mGeometryOptions = qgis::make_unique<QgsGeometryOptions>();
  mActions = new QgsActionManager( this );
  mConditionalStyles = new QgsConditionalLayerStyles();

  mJoinBuffer = new QgsVectorLayerJoinBuffer( this );
  mJoinBuffer->setParent( this );
  connect( mJoinBuffer, &QgsVectorLayerJoinBuffer::joinedFieldsChanged, this, &QgsVectorLayer::onJoinedFieldsChanged );

  mExpressionFieldBuffer = new QgsExpressionFieldBuffer();
  // if we're given a provider type, try to create and bind one to this layer
  if ( !vectorLayerPath.isEmpty() && !mProviderKey.isEmpty() )
  {
    QgsDataProvider::ProviderOptions providerOptions { options.transformContext };
    setDataSource( vectorLayerPath, baseName, providerKey, providerOptions, options.loadDefaultStyle );
  }

  connect( this, &QgsVectorLayer::selectionChanged, this, [ = ] { emit repaintRequested(); } );
  connect( QgsProject::instance()->relationManager(), &QgsRelationManager::relationsLoaded, this, &QgsVectorLayer::onRelationsLoaded );

  connect( this, &QgsVectorLayer::subsetStringChanged, this, &QgsMapLayer::configChanged );

  // Default simplify drawing settings
  QgsSettings settings;
  mSimplifyMethod.setSimplifyHints( settings.flagValue( QStringLiteral( "qgis/simplifyDrawingHints" ), mSimplifyMethod.simplifyHints(), QgsSettings::NoSection ) );
  mSimplifyMethod.setSimplifyAlgorithm( settings.enumValue( QStringLiteral( "qgis/simplifyAlgorithm" ), mSimplifyMethod.simplifyAlgorithm() ) );
  mSimplifyMethod.setThreshold( settings.value( QStringLiteral( "qgis/simplifyDrawingTol" ), mSimplifyMethod.threshold() ).toFloat() );
  mSimplifyMethod.setForceLocalOptimization( settings.value( QStringLiteral( "qgis/simplifyLocal" ), mSimplifyMethod.forceLocalOptimization() ).toBool() );
  mSimplifyMethod.setMaximumScale( settings.value( QStringLiteral( "qgis/simplifyMaxScale" ), mSimplifyMethod.maximumScale() ).toFloat() );

} // QgsVectorLayer ctor


QgsVectorLayer::~QgsVectorLayer()
{
  emit willBeDeleted();

  mValid = false;

  delete mDataProvider;
  delete mEditBuffer;
  delete mJoinBuffer;
  delete mExpressionFieldBuffer;
  delete mLabeling;
  delete mDiagramLayerSettings;
  delete mDiagramRenderer;

  delete mActions;

  delete mRenderer;
  delete mConditionalStyles;

  if ( mFeatureCounter )
    mFeatureCounter->cancel();
}

QgsVectorLayer *QgsVectorLayer::clone() const
{
  QgsVectorLayer::LayerOptions options;
  if ( mDataProvider )
  {
    options.transformContext = mDataProvider->transformContext();
  }
  QgsVectorLayer *layer = new QgsVectorLayer( source(), name(), mProviderKey, options );
  QgsMapLayer::clone( layer );

  QList<QgsVectorLayerJoinInfo> joins = vectorJoins();
  const auto constJoins = joins;
  for ( const QgsVectorLayerJoinInfo &join : constJoins )
  {
    // do not copy join information for auxiliary layer
    if ( !auxiliaryLayer()
         || ( auxiliaryLayer() && auxiliaryLayer()->id() != join.joinLayerId() ) )
      layer->addJoin( join );
  }

  layer->setProviderEncoding( dataProvider()->encoding() );
  layer->setDisplayExpression( displayExpression() );
  layer->setMapTipTemplate( mapTipTemplate() );
  layer->setReadOnly( isReadOnly() );
  layer->selectByIds( selectedFeatureIds() );
  layer->setExcludeAttributesWms( excludeAttributesWms() );
  layer->setExcludeAttributesWfs( excludeAttributesWfs() );
  layer->setAttributeTableConfig( attributeTableConfig() );
  layer->setFeatureBlendMode( featureBlendMode() );
  layer->setOpacity( opacity() );
  layer->setReadExtentFromXml( readExtentFromXml() );

  const auto constActions = actions()->actions();
  for ( const QgsAction &action : constActions )
  {
    layer->actions()->addAction( action );
  }

  if ( renderer() )
  {
    layer->setRenderer( renderer()->clone() );
  }

  if ( labeling() )
  {
    layer->setLabeling( labeling()->clone() );
  }
  layer->setLabelsEnabled( labelsEnabled() );

  layer->setSimplifyMethod( simplifyMethod() );

  if ( diagramRenderer() )
  {
    layer->setDiagramRenderer( diagramRenderer()->clone() );
  }

  if ( diagramLayerSettings() )
  {
    layer->setDiagramLayerSettings( *diagramLayerSettings() );
  }

  for ( int i = 0; i < fields().count(); i++ )
  {
    layer->setFieldAlias( i, attributeAlias( i ) );
    layer->setEditorWidgetSetup( i, editorWidgetSetup( i ) );
    layer->setConstraintExpression( i, constraintExpression( i ), constraintDescription( i ) );
    layer->setDefaultValueDefinition( i, defaultValueDefinition( i ) );

    QMap< QgsFieldConstraints::Constraint, QgsFieldConstraints::ConstraintStrength> constraints = fieldConstraintsAndStrength( i );
    auto constraintIt = constraints.constBegin();
    for ( ; constraintIt != constraints.constEnd(); ++ constraintIt )
    {
      layer->setFieldConstraint( i, constraintIt.key(), constraintIt.value() );
    }

    if ( fields().fieldOrigin( i ) == QgsFields::OriginExpression )
    {
      layer->addExpressionField( expressionField( i ), fields().at( i ) );
    }
  }

  layer->setEditFormConfig( editFormConfig() );

  if ( auxiliaryLayer() )
    layer->setAuxiliaryLayer( auxiliaryLayer()->clone( layer ) );

  return layer;
}

QString QgsVectorLayer::storageType() const
{
  if ( mDataProvider )
  {
    return mDataProvider->storageType();
  }
  return QString();
}


QString QgsVectorLayer::capabilitiesString() const
{
  if ( mDataProvider )
  {
    return mDataProvider->capabilitiesString();
  }
  return QString();
}

QString QgsVectorLayer::dataComment() const
{
  if ( mDataProvider )
  {
    return mDataProvider->dataComment();
  }
  return QString();
}

QgsCoordinateReferenceSystem QgsVectorLayer::sourceCrs() const
{
  return crs();
}

QString QgsVectorLayer::sourceName() const
{
  return name();
}

void QgsVectorLayer::reload()
{
  if ( mDataProvider )
  {
    mDataProvider->reloadData();
    updateFields();
  }
}

QgsMapLayerRenderer *QgsVectorLayer::createMapRenderer( QgsRenderContext &rendererContext )
{
  return new QgsVectorLayerRenderer( this, rendererContext );
}


void QgsVectorLayer::drawVertexMarker( double x, double y, QPainter &p, QgsVectorLayer::VertexMarkerType type, int m )
{
  if ( type == QgsVectorLayer::SemiTransparentCircle )
  {
    p.setPen( QColor( 50, 100, 120, 200 ) );
    p.setBrush( QColor( 200, 200, 210, 120 ) );
    p.drawEllipse( x - m, y - m, m * 2 + 1, m * 2 + 1 );
  }
  else if ( type == QgsVectorLayer::Cross )
  {
    p.setPen( QColor( 255, 0, 0 ) );
    p.drawLine( x - m, y + m, x + m, y - m );
    p.drawLine( x - m, y - m, x + m, y + m );
  }
}

void QgsVectorLayer::select( QgsFeatureId fid )
{
  mSelectedFeatureIds.insert( fid );

  emit selectionChanged( QgsFeatureIds() << fid, QgsFeatureIds(), false );
}

void QgsVectorLayer::select( const QgsFeatureIds &featureIds )
{
  mSelectedFeatureIds.unite( featureIds );

  emit selectionChanged( featureIds, QgsFeatureIds(), false );
}

void QgsVectorLayer::deselect( const QgsFeatureId fid )
{
  mSelectedFeatureIds.remove( fid );

  emit selectionChanged( QgsFeatureIds(), QgsFeatureIds() << fid, false );
}

void QgsVectorLayer::deselect( const QgsFeatureIds &featureIds )
{
  mSelectedFeatureIds.subtract( featureIds );

  emit selectionChanged( QgsFeatureIds(), featureIds, false );
}

void QgsVectorLayer::selectByRect( QgsRectangle &rect, QgsVectorLayer::SelectBehavior behavior )
{
  // normalize the rectangle
  rect.normalize();

  QgsFeatureIds newSelection;

  QgsFeatureIterator features = getFeatures( QgsFeatureRequest()
                                .setFilterRect( rect )
                                .setFlags( QgsFeatureRequest::ExactIntersect | QgsFeatureRequest::NoGeometry )
                                .setNoAttributes() );

  QgsFeature feat;
  while ( features.nextFeature( feat ) )
  {
    newSelection << feat.id();
  }
  features.close();

  selectByIds( newSelection, behavior );
}

void QgsVectorLayer::selectByExpression( const QString &expression, QgsVectorLayer::SelectBehavior behavior )
{
  QgsFeatureIds newSelection;

  QgsExpressionContext context( QgsExpressionContextUtils::globalProjectLayerScopes( this ) );

  if ( behavior == SetSelection || behavior == AddToSelection )
  {
    QgsFeatureRequest request = QgsFeatureRequest().setFilterExpression( expression )
                                .setExpressionContext( context )
                                .setFlags( QgsFeatureRequest::NoGeometry )
                                .setNoAttributes();

    QgsFeatureIterator features = getFeatures( request );

    if ( behavior == AddToSelection )
    {
      newSelection = selectedFeatureIds();
    }
    QgsFeature feat;
    while ( features.nextFeature( feat ) )
    {
      newSelection << feat.id();
    }
    features.close();
  }
  else if ( behavior == IntersectSelection || behavior == RemoveFromSelection )
  {
    QgsExpression exp( expression );
    exp.prepare( &context );

    QgsFeatureIds oldSelection = selectedFeatureIds();
    QgsFeatureRequest request = QgsFeatureRequest().setFilterFids( oldSelection );

    //refine request
    if ( !exp.needsGeometry() )
      request.setFlags( QgsFeatureRequest::NoGeometry );
    request.setSubsetOfAttributes( exp.referencedColumns(), fields() );

    QgsFeatureIterator features = getFeatures( request );
    QgsFeature feat;
    while ( features.nextFeature( feat ) )
    {
      context.setFeature( feat );
      bool matches = exp.evaluate( &context ).toBool();

      if ( matches && behavior == IntersectSelection )
      {
        newSelection << feat.id();
      }
      else if ( !matches && behavior == RemoveFromSelection )
      {
        newSelection << feat.id();
      }
    }
  }

  selectByIds( newSelection );
}

void QgsVectorLayer::selectByIds( const QgsFeatureIds &ids, QgsVectorLayer::SelectBehavior behavior )
{
  QgsFeatureIds newSelection;

  switch ( behavior )
  {
    case SetSelection:
      newSelection = ids;
      break;

    case AddToSelection:
      newSelection = mSelectedFeatureIds + ids;
      break;

    case RemoveFromSelection:
      newSelection = mSelectedFeatureIds - ids;
      break;

    case IntersectSelection:
      newSelection = mSelectedFeatureIds.intersect( ids );
      break;
  }

  QgsFeatureIds deselectedFeatures = mSelectedFeatureIds - newSelection;
  mSelectedFeatureIds = newSelection;

  emit selectionChanged( newSelection, deselectedFeatures, true );
}

void QgsVectorLayer::modifySelection( const QgsFeatureIds &selectIds, const QgsFeatureIds &deselectIds )
{
  QgsFeatureIds intersectingIds = selectIds & deselectIds;
  if ( !intersectingIds.isEmpty() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Trying to select and deselect the same item at the same time. Unsure what to do. Selecting dubious items." ), 3 );
  }

  mSelectedFeatureIds -= deselectIds;
  mSelectedFeatureIds += selectIds;

  emit selectionChanged( selectIds, deselectIds - intersectingIds, false );
}

void QgsVectorLayer::invertSelection()
{
  QgsFeatureIds ids = allFeatureIds();
  ids.subtract( mSelectedFeatureIds );
  selectByIds( ids );
}

void QgsVectorLayer::selectAll()
{
  selectByIds( allFeatureIds() );
}

void QgsVectorLayer::invertSelectionInRectangle( QgsRectangle &rect )
{
  // normalize the rectangle
  rect.normalize();

  QgsFeatureIterator fit = getFeatures( QgsFeatureRequest()
                                        .setFilterRect( rect )
                                        .setFlags( QgsFeatureRequest::NoGeometry | QgsFeatureRequest::ExactIntersect )
                                        .setNoAttributes() );

  QgsFeatureIds selectIds;
  QgsFeatureIds deselectIds;

  QgsFeature fet;
  while ( fit.nextFeature( fet ) )
  {
    if ( mSelectedFeatureIds.contains( fet.id() ) )
    {
      deselectIds << fet.id();
    }
    else
    {
      selectIds << fet.id();
    }
  }

  modifySelection( selectIds, deselectIds );
}

void QgsVectorLayer::removeSelection()
{
  if ( mSelectedFeatureIds.isEmpty() )
    return;

  selectByIds( QgsFeatureIds() );
}

QgsVectorDataProvider *QgsVectorLayer::dataProvider()
{
  return mDataProvider;
}

const QgsVectorDataProvider *QgsVectorLayer::dataProvider() const
{
  return mDataProvider;
}

void QgsVectorLayer::setProviderEncoding( const QString &encoding )
{
  if ( mValid && mDataProvider && mDataProvider->encoding() != encoding )
  {
    mDataProvider->setEncoding( encoding );
    updateFields();
  }
}

void QgsVectorLayer::setDiagramRenderer( QgsDiagramRenderer *r )
{
  delete mDiagramRenderer;
  mDiagramRenderer = r;
  emit rendererChanged();
  emit styleChanged();
}

QgsWkbTypes::GeometryType QgsVectorLayer::geometryType() const
{
  return QgsWkbTypes::geometryType( mWkbType );
}

QgsWkbTypes::Type QgsVectorLayer::wkbType() const
{
  return mWkbType;
}

QgsRectangle QgsVectorLayer::boundingBoxOfSelected() const
{
  if ( !mValid || !isSpatial() || mSelectedFeatureIds.isEmpty() ) //no selected features
  {
    return QgsRectangle( 0, 0, 0, 0 );
  }

  QgsRectangle r, retval;
  retval.setMinimal();

  QgsFeature fet;
  if ( mDataProvider->capabilities() & QgsVectorDataProvider::SelectAtId )
  {
    QgsFeatureIterator fit = getFeatures( QgsFeatureRequest()
                                          .setFilterFids( mSelectedFeatureIds )
                                          .setNoAttributes() );

    while ( fit.nextFeature( fet ) )
    {
      if ( !fet.hasGeometry() )
        continue;
      r = fet.geometry().boundingBox();
      retval.combineExtentWith( r );
    }
  }
  else
  {
    QgsFeatureIterator fit = getFeatures( QgsFeatureRequest()
                                          .setNoAttributes() );

    while ( fit.nextFeature( fet ) )
    {
      if ( mSelectedFeatureIds.contains( fet.id() ) )
      {
        if ( fet.hasGeometry() )
        {
          r = fet.geometry().boundingBox();
          retval.combineExtentWith( r );
        }
      }
    }
  }

  if ( retval.width() == 0.0 || retval.height() == 0.0 )
  {
    // If all of the features are at the one point, buffer the
    // rectangle a bit. If they are all at zero, do something a bit
    // more crude.

    if ( retval.xMinimum() == 0.0 && retval.xMaximum() == 0.0 &&
         retval.yMinimum() == 0.0 && retval.yMaximum() == 0.0 )
    {
      retval.set( -1.0, -1.0, 1.0, 1.0 );
    }
  }

  return retval;
}

bool QgsVectorLayer::labelsEnabled() const
{
  return mLabelsEnabled && static_cast< bool >( mLabeling );
}

void QgsVectorLayer::setLabelsEnabled( bool enabled )
{
  mLabelsEnabled = enabled;
}

bool QgsVectorLayer::diagramsEnabled() const
{
  if ( !mDiagramRenderer || !mDiagramLayerSettings )
    return false;

  QList<QgsDiagramSettings> settingList = mDiagramRenderer->diagramSettings();
  if ( !settingList.isEmpty() )
  {
    return settingList.at( 0 ).enabled;
  }
  return false;
}

long QgsVectorLayer::featureCount( const QString &legendKey ) const
{
  if ( !mSymbolFeatureCounted )
    return -1;

  return mSymbolFeatureCountMap.value( legendKey );
}



QgsVectorLayerFeatureCounter *QgsVectorLayer::countSymbolFeatures()
{
  if ( mSymbolFeatureCounted || mFeatureCounter )
    return mFeatureCounter;

  mSymbolFeatureCountMap.clear();

  if ( !mValid )
  {
    QgsDebugMsgLevel( QStringLiteral( "invoked with invalid layer" ), 3 );
    return mFeatureCounter;
  }
  if ( !mDataProvider )
  {
    QgsDebugMsgLevel( QStringLiteral( "invoked with null mDataProvider" ), 3 );
    return mFeatureCounter;
  }
  if ( !mRenderer )
  {
    QgsDebugMsgLevel( QStringLiteral( "invoked with null mRenderer" ), 3 );
    return mFeatureCounter;
  }

  if ( !mFeatureCounter )
  {
    mFeatureCounter = new QgsVectorLayerFeatureCounter( this );
    connect( mFeatureCounter, &QgsTask::taskCompleted, this, &QgsVectorLayer::onFeatureCounterCompleted );
    connect( mFeatureCounter, &QgsTask::taskTerminated, this, &QgsVectorLayer::onFeatureCounterTerminated );

    QgsApplication::taskManager()->addTask( mFeatureCounter );
  }

  return mFeatureCounter;
}

void QgsVectorLayer::updateExtents( bool force )
{
  // do not update extent by default when trust project option is activated
  if ( force || !mReadExtentFromXml || ( mReadExtentFromXml && mXmlExtent.isNull() ) )
    mValidExtent = false;
}

void QgsVectorLayer::setExtent( const QgsRectangle &r )
{
  QgsMapLayer::setExtent( r );
  mValidExtent = true;
}

void QgsVectorLayer::updateDefaultValues( QgsFeatureId fid, QgsFeature feature )
{
  if ( !mDefaultValueOnUpdateFields.isEmpty() )
  {
    if ( !feature.isValid() )
      feature = getFeature( fid );

    int size = mFields.size();
    for ( int idx : qgis::as_const( mDefaultValueOnUpdateFields ) )
    {
      if ( idx < 0 || idx >= size )
        continue;

      feature.setAttribute( idx, defaultValue( idx, feature ) );
      updateFeature( feature, true );
    }
  }
}

QgsRectangle QgsVectorLayer::extent() const
{
  QgsRectangle rect;
  rect.setMinimal();

  if ( !isSpatial() )
    return rect;


  if ( !mValidExtent && mLazyExtent && mDataProvider && !mDataProvider->hasMetadata() && mReadExtentFromXml && !mXmlExtent.isNull() )
  {
    mExtent = mXmlExtent;
    mValidExtent = true;
    mLazyExtent = false;
  }

  if ( !mValidExtent && mLazyExtent && mDataProvider && mDataProvider->isValid() )
  {
    // get the extent
    QgsRectangle mbr = mDataProvider->extent();

    // show the extent
    QgsDebugMsgLevel( QStringLiteral( "Extent of layer: %1" ).arg( mbr.toString() ), 3 );
    // store the extent
    mValidExtent = true;
    mExtent = mbr;

    mLazyExtent = false;
  }

  if ( mValidExtent )
    return QgsMapLayer::extent();

  if ( !mValid || !mDataProvider )
  {
    QgsDebugMsgLevel( QStringLiteral( "invoked with invalid layer or null mDataProvider" ), 3 );
    return rect;
  }

  if ( !mEditBuffer ||
       ( mEditBuffer->mDeletedFeatureIds.isEmpty() && mEditBuffer->mChangedGeometries.isEmpty() ) ||
       QgsDataSourceUri( mDataProvider->dataSourceUri() ).useEstimatedMetadata() )
  {
    mDataProvider->updateExtents();

    // get the extent of the layer from the provider
    // but only when there are some features already
    if ( mDataProvider->featureCount() != 0 )
    {
      QgsRectangle r = mDataProvider->extent();
      rect.combineExtentWith( r );
    }

    if ( mEditBuffer )
    {
      for ( QgsFeatureMap::const_iterator it = mEditBuffer->mAddedFeatures.constBegin(); it != mEditBuffer->mAddedFeatures.constEnd(); ++it )
      {
        if ( it->hasGeometry() )
        {
          QgsRectangle r = it->geometry().boundingBox();
          rect.combineExtentWith( r );
        }
      }
    }
  }
  else
  {
    QgsFeatureIterator fit = getFeatures( QgsFeatureRequest()
                                          .setNoAttributes() );

    QgsFeature fet;
    while ( fit.nextFeature( fet ) )
    {
      if ( fet.hasGeometry() && fet.geometry().type() != QgsWkbTypes::UnknownGeometry )
      {
        QgsRectangle bb = fet.geometry().boundingBox();
        rect.combineExtentWith( bb );
      }
    }
  }

  if ( rect.xMinimum() > rect.xMaximum() && rect.yMinimum() > rect.yMaximum() )
  {
    // special case when there are no features in provider nor any added
    rect = QgsRectangle(); // use rectangle with zero coordinates
  }

  mValidExtent = true;
  mExtent = rect;

  // Send this (hopefully) up the chain to the map canvas
  emit recalculateExtents();

  return rect;
}

QgsRectangle QgsVectorLayer::sourceExtent() const
{
  return extent();
}

QString QgsVectorLayer::subsetString() const
{
  if ( !mValid || !mDataProvider )
  {
    QgsDebugMsgLevel( QStringLiteral( "invoked with invalid layer or null mDataProvider" ), 3 );
    return QString();
  }
  return mDataProvider->subsetString();
}

bool QgsVectorLayer::setSubsetString( const QString &subset )
{
  if ( !mValid || !mDataProvider || mEditBuffer )
  {
    QgsDebugMsgLevel( QStringLiteral( "invoked with invalid layer or null mDataProvider or while editing" ), 3 );
    return false;
  }

  if ( subset == mDataProvider->subsetString() )
    return true;

  bool res = mDataProvider->setSubsetString( subset );

  // get the updated data source string from the provider
  mDataSource = mDataProvider->dataSourceUri();
  updateExtents();
  updateFields();

  if ( res )
  {
    emit subsetStringChanged();
    emit repaintRequested();
  }

  return res;
}

bool QgsVectorLayer::simplifyDrawingCanbeApplied( const QgsRenderContext &renderContext, QgsVectorSimplifyMethod::SimplifyHint simplifyHint ) const
{
  if ( mValid && mDataProvider && !mEditBuffer && ( isSpatial() && geometryType() != QgsWkbTypes::PointGeometry ) && ( mSimplifyMethod.simplifyHints() & simplifyHint ) && renderContext.useRenderingOptimization() )
  {
    double maximumSimplificationScale = mSimplifyMethod.maximumScale();

    // check maximum scale at which generalisation should be carried out
    return !( maximumSimplificationScale > 1 && renderContext.rendererScale() <= maximumSimplificationScale );
  }
  return false;
}

QgsConditionalLayerStyles *QgsVectorLayer::conditionalStyles() const
{
  return mConditionalStyles;
}

QgsFeatureIterator QgsVectorLayer::getFeatures( const QgsFeatureRequest &request ) const
{
  if ( !mValid || !mDataProvider )
    return QgsFeatureIterator();

  return QgsFeatureIterator( new QgsVectorLayerFeatureIterator( new QgsVectorLayerFeatureSource( this ), true, request ) );
}

QgsGeometry QgsVectorLayer::getGeometry( QgsFeatureId fid ) const
{
  QgsFeature feature;
  getFeatures( QgsFeatureRequest( fid ).setFlags( QgsFeatureRequest::SubsetOfAttributes ) ).nextFeature( feature );
  if ( feature.isValid() )
    return feature.geometry();
  else
    return QgsGeometry();
}

bool QgsVectorLayer::addFeature( QgsFeature &feature, Flags )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return false;


  if ( mGeometryOptions->isActive() )
  {
    QgsGeometry geom = feature.geometry();
    mGeometryOptions->apply( geom );
    feature.setGeometry( geom );
  }

  bool success = mEditBuffer->addFeature( feature );

  if ( success )
  {
    updateExtents();

    if ( mJoinBuffer->containsJoins() )
      success = mJoinBuffer->addFeature( feature );
  }

  return success;
}

bool QgsVectorLayer::updateFeature( QgsFeature &updatedFeature, bool skipDefaultValues )
{
  if ( !mEditBuffer || !mDataProvider )
  {
    return false;
  }

  QgsFeature currentFeature = getFeature( updatedFeature.id() );
  if ( currentFeature.isValid() )
  {
    bool hasChanged = false;
    bool hasError = false;

    if ( ( updatedFeature.hasGeometry() || currentFeature.hasGeometry() ) && !updatedFeature.geometry().equals( currentFeature.geometry() ) )
    {
      QgsGeometry geometry = updatedFeature.geometry();
      if ( changeGeometry( updatedFeature.id(), geometry, true ) )
      {
        hasChanged = true;
        updatedFeature.setGeometry( geometry );
      }
      else
      {
        QgsDebugMsgLevel( QStringLiteral( "geometry of feature %1 could not be changed." ).arg( updatedFeature.id() ), 3 );
      }
    }

    QgsAttributes fa = updatedFeature.attributes();
    QgsAttributes ca = currentFeature.attributes();

    for ( int attr = 0; attr < fa.count(); ++attr )
    {
      if ( fa.at( attr ) != ca.at( attr ) )
      {
        if ( changeAttributeValue( updatedFeature.id(), attr, fa.at( attr ), ca.at( attr ), true ) )
        {
          hasChanged = true;
        }
        else
        {
          QgsDebugMsgLevel( QStringLiteral( "attribute %1 of feature %2 could not be changed." ).arg( attr ).arg( updatedFeature.id() ), 3 );
          hasError = true;
        }
      }
    }
    if ( hasChanged && !mDefaultValueOnUpdateFields.isEmpty() && !skipDefaultValues )
      updateDefaultValues( updatedFeature.id(), updatedFeature );

    return !hasError;
  }
  else
  {
    QgsDebugMsgLevel( QStringLiteral( "feature %1 could not be retrieved" ).arg( updatedFeature.id() ), 3 );
    return false;
  }
}


bool QgsVectorLayer::insertVertex( double x, double y, QgsFeatureId atFeatureId, int beforeVertex )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return false;

  QgsVectorLayerEditUtils utils( this );
  bool result = utils.insertVertex( x, y, atFeatureId, beforeVertex );
  if ( result )
    updateExtents();
  return result;
}


bool QgsVectorLayer::insertVertex( const QgsPoint &point, QgsFeatureId atFeatureId, int beforeVertex )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return false;

  QgsVectorLayerEditUtils utils( this );
  bool result = utils.insertVertex( point, atFeatureId, beforeVertex );
  if ( result )
    updateExtents();
  return result;
}


bool QgsVectorLayer::moveVertex( double x, double y, QgsFeatureId atFeatureId, int atVertex )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return false;

  QgsVectorLayerEditUtils utils( this );
  bool result = utils.moveVertex( x, y, atFeatureId, atVertex );

  if ( result )
    updateExtents();
  return result;
}

bool QgsVectorLayer::moveVertex( const QgsPoint &p, QgsFeatureId atFeatureId, int atVertex )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return false;

  QgsVectorLayerEditUtils utils( this );
  bool result = utils.moveVertex( p, atFeatureId, atVertex );

  if ( result )
    updateExtents();
  return result;
}

QgsVectorLayer::EditResult QgsVectorLayer::deleteVertex( QgsFeatureId featureId, int vertex )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsVectorLayer::InvalidLayer;

  QgsVectorLayerEditUtils utils( this );
  EditResult result = utils.deleteVertex( featureId, vertex );

  if ( result == Success )
    updateExtents();
  return result;
}


bool QgsVectorLayer::deleteSelectedFeatures( int *deletedCount )
{
  if ( !mValid || !mDataProvider || !( mDataProvider->capabilities() & QgsVectorDataProvider::DeleteFeatures ) )
  {
    return false;
  }

  if ( !isEditable() )
  {
    return false;
  }

  int deleted = 0;
  int count = mSelectedFeatureIds.size();
  // Make a copy since deleteFeature modifies mSelectedFeatureIds
  QgsFeatureIds selectedFeatures( mSelectedFeatureIds );
  const auto constSelectedFeatures = selectedFeatures;
  for ( QgsFeatureId fid : constSelectedFeatures )
  {
    deleted += deleteFeature( fid );  // removes from selection
  }

  triggerRepaint();
  updateExtents();

  if ( deletedCount )
  {
    *deletedCount = deleted;
  }

  return deleted == count;
}

QgsGeometry::OperationResult QgsVectorLayer::addRing( const QVector<QgsPointXY> &ring, QgsFeatureId *featureId )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsGeometry::OperationResult::LayerNotEditable;

  QgsVectorLayerEditUtils utils( this );
  QgsGeometry::OperationResult result = QgsGeometry::OperationResult::AddRingNotInExistingFeature;

  //first try with selected features
  if ( !mSelectedFeatureIds.isEmpty() )
  {
    result = utils.addRing( ring, mSelectedFeatureIds, featureId );
  }

  if ( result != QgsGeometry::OperationResult::Success )
  {
    //try with all intersecting features
    result = utils.addRing( ring, QgsFeatureIds(), featureId );
  }

  return result;
}

QgsGeometry::OperationResult QgsVectorLayer::addRing( QgsCurve *ring, QgsFeatureId *featureId )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
  {
    delete ring;
    return QgsGeometry::OperationResult::LayerNotEditable;
  }

  if ( !ring )
  {
    return QgsGeometry::OperationResult::InvalidInputGeometryType;
  }

  if ( !ring->isClosed() )
  {
    delete ring;
    return QgsGeometry::OperationResult::AddRingNotClosed;
  }

  QgsVectorLayerEditUtils utils( this );
  QgsGeometry::OperationResult result = QgsGeometry::OperationResult::AddRingNotInExistingFeature;

  //first try with selected features
  if ( !mSelectedFeatureIds.isEmpty() )
  {
    result = utils.addRing( static_cast< QgsCurve * >( ring->clone() ), mSelectedFeatureIds, featureId );
  }

  if ( result != QgsGeometry::OperationResult::Success )
  {
    //try with all intersecting features
    result = utils.addRing( static_cast< QgsCurve * >( ring->clone() ), QgsFeatureIds(), featureId );
  }

  delete ring;
  return result;
}

QgsGeometry::OperationResult QgsVectorLayer::addPart( const QList<QgsPointXY> &points )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsGeometry::OperationResult::LayerNotEditable;

  //number of selected features must be 1

  if ( mSelectedFeatureIds.empty() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Number of selected features < 1" ), 3 );
    return QgsGeometry::OperationResult::SelectionIsEmpty;
  }
  else if ( mSelectedFeatureIds.size() > 1 )
  {
    QgsDebugMsgLevel( QStringLiteral( "Number of selected features > 1" ), 3 );
    return QgsGeometry::OperationResult::SelectionIsGreaterThanOne;
  }

  QgsVectorLayerEditUtils utils( this );
  QgsGeometry::OperationResult result = utils.addPart( points, *mSelectedFeatureIds.constBegin() );

  if ( result == QgsGeometry::OperationResult::Success )
    updateExtents();
  return result;
}

QgsGeometry::OperationResult QgsVectorLayer::addPart( const QgsPointSequence &points )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsGeometry::OperationResult::LayerNotEditable;

  //number of selected features must be 1

  if ( mSelectedFeatureIds.empty() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Number of selected features <1" ), 3 );
    return QgsGeometry::OperationResult::SelectionIsEmpty;
  }
  else if ( mSelectedFeatureIds.size() > 1 )
  {
    QgsDebugMsgLevel( QStringLiteral( "Number of selected features >1" ), 3 );
    return QgsGeometry::OperationResult::SelectionIsGreaterThanOne;
  }

  QgsVectorLayerEditUtils utils( this );
  QgsGeometry::OperationResult result = utils.addPart( points, *mSelectedFeatureIds.constBegin() );

  if ( result == QgsGeometry::OperationResult::Success )
    updateExtents();
  return result;
}

QgsGeometry::OperationResult QgsVectorLayer::addPart( QgsCurve *ring )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsGeometry::OperationResult::LayerNotEditable;

  //number of selected features must be 1

  if ( mSelectedFeatureIds.empty() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Number of selected features <1" ), 3 );
    return QgsGeometry::OperationResult::SelectionIsEmpty;
  }
  else if ( mSelectedFeatureIds.size() > 1 )
  {
    QgsDebugMsgLevel( QStringLiteral( "Number of selected features >1" ), 3 );
    return QgsGeometry::OperationResult::SelectionIsGreaterThanOne;
  }

  QgsVectorLayerEditUtils utils( this );
  QgsGeometry::OperationResult result = utils.addPart( ring, *mSelectedFeatureIds.constBegin() );

  if ( result == QgsGeometry::OperationResult::Success )
    updateExtents();
  return result;
}

int QgsVectorLayer::translateFeature( QgsFeatureId featureId, double dx, double dy )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsGeometry::OperationResult::LayerNotEditable;

  QgsVectorLayerEditUtils utils( this );
  int result = utils.translateFeature( featureId, dx, dy );

  if ( result == QgsGeometry::OperationResult::Success )
    updateExtents();
  return result;
}

QgsGeometry::OperationResult QgsVectorLayer::splitParts( const QVector<QgsPointXY> &splitLine, bool topologicalEditing )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsGeometry::OperationResult::LayerNotEditable;

  QgsVectorLayerEditUtils utils( this );
  return utils.splitParts( splitLine, topologicalEditing );
}

QgsGeometry::OperationResult QgsVectorLayer::splitFeatures( const QVector<QgsPointXY> &splitLine, bool topologicalEditing )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return QgsGeometry::OperationResult::LayerNotEditable;

  QgsVectorLayerEditUtils utils( this );
  return utils.splitFeatures( splitLine, topologicalEditing );
}

int QgsVectorLayer::addTopologicalPoints( const QgsGeometry &geom )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return -1;

  QgsVectorLayerEditUtils utils( this );
  return utils.addTopologicalPoints( geom );
}

int QgsVectorLayer::addTopologicalPoints( const QgsPointXY &p )
{
  if ( !mValid || !mEditBuffer || !mDataProvider )
    return -1;

  QgsVectorLayerEditUtils utils( this );
  return utils.addTopologicalPoints( p );
}

void QgsVectorLayer::setLabeling( QgsAbstractVectorLayerLabeling *labeling )
{
  if ( mLabeling == labeling )
    return;

  delete mLabeling;
  mLabeling = labeling;
}

bool QgsVectorLayer::startEditing()
{
  if ( !mValid || !mDataProvider )
  {
    return false;
  }

  // allow editing if provider supports any of the capabilities
  if ( !( mDataProvider->capabilities() & QgsVectorDataProvider::EditingCapabilities ) )
  {
    return false;
  }

  if ( mReadOnly )
  {
    return false;
  }

  if ( mEditBuffer )
  {
    // editing already underway
    return false;
  }

  emit beforeEditingStarted();

  mDataProvider->enterUpdateMode();

  if ( mDataProvider->transaction() )
  {
    mEditBuffer = new QgsVectorLayerEditPassthrough( this );

    connect( mDataProvider->transaction(), &QgsTransaction::dirtied, this, &QgsVectorLayer::onDirtyTransaction, Qt::UniqueConnection );
  }
  else
  {
    mEditBuffer = new QgsVectorLayerEditBuffer( this );
  }
  // forward signals
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::layerModified, this, &QgsVectorLayer::invalidateSymbolCountedFlag );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::layerModified, this, &QgsVectorLayer::layerModified ); // TODO[MD]: necessary?
  //connect( mEditBuffer, SIGNAL( layerModified() ), this, SLOT( triggerRepaint() ) ); // TODO[MD]: works well?
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::featureAdded, this, &QgsVectorLayer::featureAdded );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::featureDeleted, this, &QgsVectorLayer::onFeatureDeleted );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::geometryChanged, this, &QgsVectorLayer::geometryChanged );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::attributeValueChanged, this, &QgsVectorLayer::attributeValueChanged );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::attributeAdded, this, &QgsVectorLayer::attributeAdded );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::attributeDeleted, this, &QgsVectorLayer::attributeDeleted );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::committedAttributesDeleted, this, &QgsVectorLayer::committedAttributesDeleted );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::committedAttributesAdded, this, &QgsVectorLayer::committedAttributesAdded );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::committedFeaturesAdded, this, &QgsVectorLayer::committedFeaturesAdded );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::committedFeaturesRemoved, this, &QgsVectorLayer::committedFeaturesRemoved );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::committedAttributeValuesChanges, this, &QgsVectorLayer::committedAttributeValuesChanges );
  connect( mEditBuffer, &QgsVectorLayerEditBuffer::committedGeometriesChanges, this, &QgsVectorLayer::committedGeometriesChanges );

  updateFields();

  emit editingStarted();

  return true;
}

void QgsVectorLayer::setTransformContext( const QgsCoordinateTransformContext &transformContext )
{
  if ( mDataProvider )
    mDataProvider->setTransformContext( transformContext );
}

bool QgsVectorLayer::readXml( const QDomNode &layer_node, QgsReadWriteContext &context )
{
  QgsDebugMsgLevel( QStringLiteral( "Datasource in QgsVectorLayer::readXml: %1" ).arg( mDataSource.toLocal8Bit().data() ), 3 );

  //process provider key
  QDomNode pkeyNode = layer_node.namedItem( QStringLiteral( "provider" ) );

  if ( pkeyNode.isNull() )
  {
    mProviderKey.clear();
  }
  else
  {
    QDomElement pkeyElt = pkeyNode.toElement();
    mProviderKey = pkeyElt.text();
  }

  // determine type of vector layer
  if ( !mProviderKey.isNull() )
  {
    // if the provider string isn't empty, then we successfully
    // got the stored provider
  }
  else if ( mDataSource.contains( QLatin1String( "dbname=" ) ) )
  {
    mProviderKey = QStringLiteral( "postgres" );
  }
  else
  {
    mProviderKey = QStringLiteral( "ogr" );
  }

  QgsDataProvider::ProviderOptions options { context.transformContext() };
  if ( !setDataProvider( mProviderKey, options ) )
  {
    QgsDebugMsg( QStringLiteral( "Could not set data provider for layer %1" ).arg( publicSource() ) );
    const QDomElement elem = layer_node.toElement();

    // for invalid layer sources, we fallback to stored wkbType if available
    if ( elem.hasAttribute( QStringLiteral( "wkbType" ) ) )
      mWkbType = qgsEnumKeyToValue( elem.attribute( QStringLiteral( "wkbType" ) ), mWkbType );
  }

  QDomElement pkeyElem = pkeyNode.toElement();
  if ( !pkeyElem.isNull() )
  {
    QString encodingString = pkeyElem.attribute( QStringLiteral( "encoding" ) );
    if ( mDataProvider && !encodingString.isEmpty() )
    {
      mDataProvider->setEncoding( encodingString );
    }
  }

  // load vector joins - does not resolve references to layers yet
  mJoinBuffer->readXml( layer_node );

  updateFields();

  QString errorMsg;
  if ( !readSymbology( layer_node, errorMsg, context ) )
  {
    return false;
  }

  readStyleManager( layer_node );

  QDomNode depsNode = layer_node.namedItem( QStringLiteral( "dataDependencies" ) );
  QDomNodeList depsNodes = depsNode.childNodes();
  QSet<QgsMapLayerDependency> sources;
  for ( int i = 0; i < depsNodes.count(); i++ )
  {
    QString source = depsNodes.at( i ).toElement().attribute( QStringLiteral( "id" ) );
    sources << QgsMapLayerDependency( source );
  }
  setDependencies( sources );

  QgsMapLayerLegend *legend = QgsMapLayerLegend::defaultVectorLegend( this );
  QDomElement legendElem = layer_node.firstChildElement( QStringLiteral( "legend" ) );
  if ( !legendElem.isNull() )
    legend->readXml( legendElem, context );
  setLegend( legend );

  // read extent
  if ( mReadExtentFromXml )
  {
    QDomNode extentNode = layer_node.namedItem( QStringLiteral( "extent" ) );
    if ( !extentNode.isNull() )
    {
      mXmlExtent = QgsXmlUtils::readRectangle( extentNode.toElement() );
    }
  }

  // auxiliary layer
  const QDomNode asNode = layer_node.namedItem( QStringLiteral( "auxiliaryLayer" ) );
  const QDomElement asElem = asNode.toElement();
  if ( !asElem.isNull() )
  {
    mAuxiliaryLayerKey = asElem.attribute( QStringLiteral( "key" ) );
  }

  return mValid;               // should be true if read successfully

} // void QgsVectorLayer::readXml


void QgsVectorLayer::setDataSource( const QString &dataSource, const QString &baseName, const QString &provider, bool loadDefaultStyleFlag )
{
  QgsDataProvider::ProviderOptions options;
  setDataSource( dataSource, baseName, provider, options, loadDefaultStyleFlag );
}

void QgsVectorLayer::setDataSource( const QString &dataSource, const QString &baseName, const QString &provider, const QgsDataProvider::ProviderOptions &options, bool loadDefaultStyleFlag )
{
  QgsWkbTypes::GeometryType geomType = geometryType();

  mDataSource = dataSource;
  setName( baseName );
  setDataProvider( provider, options );

  if ( !mValid )
  {
    emit dataSourceChanged();
    return;
  }

  // Always set crs
  setCoordinateSystem();

  // reset style if loading default style, style is missing, or geometry type is has changed (and layer is valid)
  if ( !renderer() || !legend() || ( mValid && geomType != geometryType() ) || loadDefaultStyleFlag )
  {
    bool defaultLoadedFlag = false;

    if ( loadDefaultStyleFlag && isSpatial() && mDataProvider->capabilities() & QgsVectorDataProvider::CreateRenderer )
    {
      // first try to create a renderer directly from the data provider
      std::unique_ptr< QgsFeatureRenderer > defaultRenderer( mDataProvider->createRenderer() );
      if ( defaultRenderer )
      {
        defaultLoadedFlag = true;
        setRenderer( defaultRenderer.release() );
      }
    }

    // else check if there is a default style / propertysheet defined
    // for this layer and if so apply it
    if ( !defaultLoadedFlag && loadDefaultStyleFlag )
    {
      loadDefaultStyle( defaultLoadedFlag );
    }

    // if the default style failed to load or was disabled use some very basic defaults
    if ( !defaultLoadedFlag && isSpatial() )
    {
      // add single symbol renderer
      setRenderer( QgsFeatureRenderer::defaultRenderer( geometryType() ) );
    }

    setLegend( QgsMapLayerLegend::defaultVectorLegend( this ) );

    if ( mDataProvider->capabilities() & QgsVectorDataProvider::CreateLabeling )
    {
      std::unique_ptr< QgsAbstractVectorLayerLabeling > defaultLabeling( mDataProvider->createLabeling() );
      if ( defaultLabeling )
      {
        setLabeling( defaultLabeling.release() );
        setLabelsEnabled( true );
      }
    }
  }

  emit dataSourceChanged();
  emit repaintRequested();
}

QString QgsVectorLayer::loadDefaultStyle( bool &resultFlag )
{
  if ( isSpatial() && mDataProvider->capabilities() & QgsVectorDataProvider::CreateRenderer )
  {
    // first try to create a renderer directly from the data provider
    std::unique_ptr< QgsFeatureRenderer > defaultRenderer( mDataProvider->createRenderer() );
    if ( defaultRenderer )
    {
      resultFlag = true;
      setRenderer( defaultRenderer.release() );
      return QString();
    }
  }

  return QgsMapLayer::loadDefaultStyle( resultFlag );
}


bool QgsVectorLayer::setDataProvider( QString const &provider, const QgsDataProvider::ProviderOptions &options )
{
  mProviderKey = provider;
  delete mDataProvider;

  // For Postgres provider primary key unicity is tested at construction time,
  // so it has to be set before initializing the provider,
  // this manipulation is necessary to preserve default behavior when
  // "trust layer metadata" project level option is set and checkPrimaryKeyUnicity
  // was not explicitly passed in the uri
  if ( provider.compare( QLatin1String( "postgres" ) ) == 0 )
  {
    const QString checkUnicityKey { QStringLiteral( "checkPrimaryKeyUnicity" ) };
    QgsDataSourceUri uri( mDataSource );
    if ( ! uri.hasParam( checkUnicityKey ) )
    {
      uri.setParam( checkUnicityKey, mReadExtentFromXml ? "0" : "1" );
      mDataSource = uri.uri( false );
    }
  }

  mDataProvider = qobject_cast<QgsVectorDataProvider *>( QgsProviderRegistry::instance()->createProvider( provider, mDataSource, options ) );
  if ( !mDataProvider )
  {
    mValid = false;
    QgsDebugMsgLevel( QStringLiteral( "Unable to get data provider" ), 2 );
    return false;
  }

  mDataProvider->setParent( this );
  connect( mDataProvider, &QgsVectorDataProvider::raiseError, this, &QgsVectorLayer::raiseError );

  QgsDebugMsgLevel( QStringLiteral( "Instantiated the data provider plugin" ), 2 );

  mValid = mDataProvider->isValid();
  if ( !mValid )
  {
    QgsDebugMsgLevel( QStringLiteral( "Invalid provider plugin %1" ).arg( QString( mDataSource.toUtf8() ) ), 2 );
    return false;
  }

  if ( mDataProvider->capabilities() & QgsVectorDataProvider::ReadLayerMetadata )
  {
    setMetadata( mDataProvider->layerMetadata() );
    QgsDebugMsgLevel( QStringLiteral( "Set Data provider QgsLayerMetadata identifier[%1]" ).arg( metadata().identifier() ), 4 );
  }

  // TODO: Check if the provider has the capability to send fullExtentCalculated
  connect( mDataProvider, &QgsVectorDataProvider::fullExtentCalculated, this, [ = ] { updateExtents(); } );

  // get and store the feature type
  mWkbType = mDataProvider->wkbType();

  updateFields();

  if ( mProviderKey == QLatin1String( "postgres" ) )
  {
    QgsDebugMsgLevel( QStringLiteral( "Beautifying layer name %1" ).arg( name() ), 3 );

    // adjust the display name for postgres layers
    QRegExp reg( R"lit("[^"]+"\."([^"] + )"( \([^)]+\))?)lit" );
    if ( reg.indexIn( name() ) >= 0 )
    {
      QStringList stuff = reg.capturedTexts();
      QString lName = stuff[1];

      const QMap<QString, QgsMapLayer *> &layers = QgsProject::instance()->mapLayers();

      QMap<QString, QgsMapLayer *>::const_iterator it;
      for ( it = layers.constBegin(); it != layers.constEnd() && ( *it )->name() != lName; ++it )
        ;

      if ( it != layers.constEnd() && stuff.size() > 2 )
      {
        lName += '.' + stuff[2].mid( 2, stuff[2].length() - 3 );
      }

      if ( !lName.isEmpty() )
        setName( lName );
    }
    QgsDebugMsgLevel( QStringLiteral( "Beautified layer name %1" ).arg( name() ), 3 );
  }
  else if ( mProviderKey == QLatin1String( "osm" ) )
  {
    // make sure that the "observer" has been removed from URI to avoid crashes
    mDataSource = mDataProvider->dataSourceUri();
  }
  else if ( provider == QLatin1String( "ogr" ) )
  {
    // make sure that the /vsigzip or /vsizip is added to uri, if applicable
    mDataSource = mDataProvider->dataSourceUri();
    if ( mDataSource.right( 10 ) == QLatin1String( "|layerid=0" ) )
      mDataSource.chop( 10 );
  }
  else if ( provider == QStringLiteral( "memory" ) )
  {
    // required so that source differs between memory layers
    mDataSource = mDataSource + QStringLiteral( "&uid=%1" ).arg( QUuid::createUuid().toString() );
  }

  connect( mDataProvider, &QgsVectorDataProvider::dataChanged, this, &QgsVectorLayer::dataChanged );
  connect( mDataProvider, &QgsVectorDataProvider::dataChanged, this, &QgsVectorLayer::removeSelection );

  return true;
} // QgsVectorLayer:: setDataProvider




/* virtual */
bool QgsVectorLayer::writeXml( QDomNode &layer_node,
                               QDomDocument &document,
                               const QgsReadWriteContext &context ) const
{
  // first get the layer element so that we can append the type attribute

  QDomElement mapLayerNode = layer_node.toElement();

  if ( mapLayerNode.isNull() || ( "maplayer" != mapLayerNode.nodeName() ) )
  {
    QgsDebugMsgLevel( QStringLiteral( "can't find <maplayer>" ), 2 );
    return false;
  }

  mapLayerNode.setAttribute( QStringLiteral( "type" ), QStringLiteral( "vector" ) );

  // set the geometry type
  mapLayerNode.setAttribute( QStringLiteral( "geometry" ), QgsWkbTypes::geometryDisplayString( geometryType() ) );
  mapLayerNode.setAttribute( QStringLiteral( "wkbType" ), qgsEnumValueToKey( wkbType() ) );

  // add provider node
  if ( mDataProvider )
  {
    QDomElement provider  = document.createElement( QStringLiteral( "provider" ) );
    provider.setAttribute( QStringLiteral( "encoding" ), mDataProvider->encoding() );
    QDomText providerText = document.createTextNode( providerType() );
    provider.appendChild( providerText );
    layer_node.appendChild( provider );
  }

  //save joins
  mJoinBuffer->writeXml( layer_node, document );

  // dependencies
  QDomElement dependenciesElement = document.createElement( QStringLiteral( "layerDependencies" ) );
  const auto constDependencies = dependencies();
  for ( const QgsMapLayerDependency &dep : constDependencies )
  {
    if ( dep.type() != QgsMapLayerDependency::PresenceDependency )
      continue;
    QDomElement depElem = document.createElement( QStringLiteral( "layer" ) );
    depElem.setAttribute( QStringLiteral( "id" ), dep.layerId() );
    dependenciesElement.appendChild( depElem );
  }
  layer_node.appendChild( dependenciesElement );

  // change dependencies
  QDomElement dataDependenciesElement = document.createElement( QStringLiteral( "dataDependencies" ) );
  for ( const QgsMapLayerDependency &dep : constDependencies )
  {
    if ( dep.type() != QgsMapLayerDependency::DataDependency )
      continue;
    QDomElement depElem = document.createElement( QStringLiteral( "layer" ) );
    depElem.setAttribute( QStringLiteral( "id" ), dep.layerId() );
    dataDependenciesElement.appendChild( depElem );
  }
  layer_node.appendChild( dataDependenciesElement );

  // legend
  if ( legend() )
  {
    QDomElement legendElement = legend()->writeXml( document, context );
    if ( !legendElement.isNull() )
      layer_node.appendChild( legendElement );
  }

  // save expression fields
  mExpressionFieldBuffer->writeXml( layer_node, document );

  writeStyleManager( layer_node, document );

  // auxiliary layer
  QDomElement asElem = document.createElement( QStringLiteral( "auxiliaryLayer" ) );
  if ( mAuxiliaryLayer )
  {
    const QString pkField = mAuxiliaryLayer->joinInfo().targetFieldName();
    asElem.setAttribute( QStringLiteral( "key" ), pkField );
  }
  layer_node.appendChild( asElem );

  // renderer specific settings
  QString errorMsg;
  return writeSymbology( layer_node, document, errorMsg, context );
}

QString QgsVectorLayer::encodedSource( const QString &source, const QgsReadWriteContext &context ) const
{
  QString src( source );

  // TODO: what about postgres, mysql and others, they should not go through writePath()
  if ( providerType() == QLatin1String( "spatialite" ) )
  {
    QgsDataSourceUri uri( src );
    QString database = context.pathResolver().writePath( uri.database() );
    uri.setConnection( uri.host(), uri.port(), database, uri.username(), uri.password() );
    src = uri.uri();
  }
  else if ( providerType() == QLatin1String( "ogr" ) )
  {
    QStringList theURIParts = src.split( '|' );
    theURIParts[0] = context.pathResolver().writePath( theURIParts[0] );
    src = theURIParts.join( QStringLiteral( "|" ) );
  }
  else if ( providerType() == QLatin1String( "gpx" ) )
  {
    QStringList theURIParts = src.split( '?' );
    theURIParts[0] = context.pathResolver().writePath( theURIParts[0] );
    src = theURIParts.join( QStringLiteral( "?" ) );
  }
  else if ( providerType() == QLatin1String( "delimitedtext" ) )
  {
    QUrl urlSource = QUrl::fromEncoded( src.toLatin1() );
    QUrl urlDest = QUrl::fromLocalFile( context.pathResolver().writePath( urlSource.toLocalFile() ) );
    urlDest.setQueryItems( urlSource.queryItems() );
    src = QString::fromLatin1( urlDest.toEncoded() );
  }
  else if ( providerType() == QLatin1String( "memory" ) )
  {
    // Refetch the source from the provider, because adding fields actually changes the source for this provider.
    src = dataProvider()->dataSourceUri();
  }
  else
  {
    src = context.pathResolver().writePath( src );
  }

  return src;
}

QString QgsVectorLayer::decodedSource( const QString &source, const QString &provider, const QgsReadWriteContext &context ) const
{
  QString src( source );

  if ( provider == QLatin1String( "spatialite" ) )
  {
    QgsDataSourceUri uri( src );
    uri.setDatabase( context.pathResolver().readPath( uri.database() ) );
    src = uri.uri();
  }
  else if ( provider == QLatin1String( "ogr" ) )
  {
    QStringList theURIParts = src.split( '|' );
    theURIParts[0] = context.pathResolver().readPath( theURIParts[0] );
    src = theURIParts.join( QStringLiteral( "|" ) );
  }
  else if ( provider == QLatin1String( "gpx" ) )
  {
    QStringList theURIParts = src.split( '?' );
    theURIParts[0] = context.pathResolver().readPath( theURIParts[0] );
    src = theURIParts.join( QStringLiteral( "?" ) );
  }
  else if ( provider == QLatin1String( "delimitedtext" ) )
  {
    QUrl urlSource = QUrl::fromEncoded( src.toLatin1() );

    if ( !src.startsWith( QLatin1String( "file:" ) ) )
    {
      QUrl file = QUrl::fromLocalFile( src.left( src.indexOf( '?' ) ) );
      urlSource.setScheme( QStringLiteral( "file" ) );
      urlSource.setPath( file.path() );
    }

    QUrl urlDest = QUrl::fromLocalFile( context.pathResolver().readPath( urlSource.toLocalFile() ) );
    urlDest.setQueryItems( urlSource.queryItems() );
    src = QString::fromLatin1( urlDest.toEncoded() );
  }
  else
  {
    src = context.pathResolver().readPath( src );
  }

  return src;
}



void QgsVectorLayer::resolveReferences( QgsProject *project )
{
  QgsMapLayer::resolveReferences( project );
  mJoinBuffer->resolveReferences( project );
}


bool QgsVectorLayer::readSymbology( const QDomNode &layerNode, QString &errorMessage,
                                    QgsReadWriteContext &context, QgsMapLayer::StyleCategories categories )
{
  QgsReadWriteContextCategoryPopper p = context.enterCategory( tr( "Symbology" ) );

  if ( categories.testFlag( Fields ) )
  {
    if ( !mExpressionFieldBuffer )
      mExpressionFieldBuffer = new QgsExpressionFieldBuffer();
    mExpressionFieldBuffer->readXml( layerNode );

    updateFields();
  }

  QDomElement layerElement = layerNode.toElement();

  readCommonStyle( layerElement, context, categories );

  readStyle( layerNode, errorMessage, context, categories );

  if ( categories.testFlag( MapTips ) )
    mMapTipTemplate = layerNode.namedItem( QStringLiteral( "mapTip" ) ).toElement().text();

  if ( categories.testFlag( LayerConfiguration ) )
    mDisplayExpression = layerNode.namedItem( QStringLiteral( "previewExpression" ) ).toElement().text();

  // Try to migrate pre QGIS 3.0 display field property
  QString displayField = layerNode.namedItem( QStringLiteral( "displayfield" ) ).toElement().text();
  if ( mFields.lookupField( displayField ) < 0 )
  {
    // if it's not a field, it's a maptip
    if ( mMapTipTemplate.isEmpty() && categories.testFlag( MapTips ) )
      mMapTipTemplate = displayField;
  }
  else
  {
    if ( mDisplayExpression.isEmpty() && categories.testFlag( LayerConfiguration ) )
      mDisplayExpression = QgsExpression::quotedColumnRef( displayField );
  }

  // process the attribute actions
  if ( categories.testFlag( Actions ) )
    mActions->readXml( layerNode );

  if ( categories.testFlag( Fields ) )
  {
    mAttributeAliasMap.clear();
    QDomNode aliasesNode = layerNode.namedItem( QStringLiteral( "aliases" ) );
    if ( !aliasesNode.isNull() )
    {
      QDomElement aliasElem;

      QDomNodeList aliasNodeList = aliasesNode.toElement().elementsByTagName( QStringLiteral( "alias" ) );
      for ( int i = 0; i < aliasNodeList.size(); ++i )
      {
        aliasElem = aliasNodeList.at( i ).toElement();

        QString field;
        if ( aliasElem.hasAttribute( QStringLiteral( "field" ) ) )
        {
          field = aliasElem.attribute( QStringLiteral( "field" ) );
        }
        else
        {
          int index = aliasElem.attribute( QStringLiteral( "index" ) ).toInt();

          if ( index >= 0 && index < fields().count() )
            field = fields().at( index ).name();
        }

        QString alias;

        if ( !aliasElem.attribute( QStringLiteral( "name" ) ).isEmpty() )
        {
          //if it has alias
          alias = context.projectTranslator()->translate( QStringLiteral( "project:layers:%1:fieldaliases" ).arg( layerNode.namedItem( QStringLiteral( "id" ) ).toElement().text() ), aliasElem.attribute( QStringLiteral( "name" ) ) );
          QgsDebugMsgLevel( "context" + QStringLiteral( "project:layers:%1:fieldaliases" ).arg( layerNode.namedItem( QStringLiteral( "id" ) ).toElement().text() ) + " source " + aliasElem.attribute( QStringLiteral( "name" ) ), 3 );
        }
        else
        {
          //if it has no alias, it should be the fields translation
          alias = context.projectTranslator()->translate( QStringLiteral( "project:layers:%1:fieldaliases" ).arg( layerNode.namedItem( QStringLiteral( "id" ) ).toElement().text() ), field );
          QgsDebugMsgLevel( "context" + QStringLiteral( "project:layers:%1:fieldaliases" ).arg( layerNode.namedItem( QStringLiteral( "id" ) ).toElement().text() ) + " source " + field, 3 );
          //if it gets the exact field value, there has been no translation (or not even translation loaded) - so no alias should be generated;
          if ( alias == aliasElem.attribute( QStringLiteral( "field" ) ) )
            alias.clear();
        }

        QgsDebugMsgLevel( "field " + field + " origalias " + aliasElem.attribute( QStringLiteral( "name" ) ) + " trans " + alias, 3 );
        mAttributeAliasMap.insert( field, alias );
      }
    }

    // default expressions
    mDefaultExpressionMap.clear();
    QDomNode defaultsNode = layerNode.namedItem( QStringLiteral( "defaults" ) );
    if ( !defaultsNode.isNull() )
    {
      QDomNodeList defaultNodeList = defaultsNode.toElement().elementsByTagName( QStringLiteral( "default" ) );
      for ( int i = 0; i < defaultNodeList.size(); ++i )
      {
        QDomElement defaultElem = defaultNodeList.at( i ).toElement();

        QString field = defaultElem.attribute( QStringLiteral( "field" ), QString() );
        QString expression = defaultElem.attribute( QStringLiteral( "expression" ), QString() );
        bool applyOnUpdate = defaultElem.attribute( QStringLiteral( "applyOnUpdate" ), QStringLiteral( "0" ) ) == QLatin1String( "1" );
        if ( field.isEmpty() || expression.isEmpty() )
          continue;

        mDefaultExpressionMap.insert( field, QgsDefaultValue( expression, applyOnUpdate ) );
      }
    }

    // constraints
    mFieldConstraints.clear();
    mFieldConstraintStrength.clear();
    QDomNode constraintsNode = layerNode.namedItem( QStringLiteral( "constraints" ) );
    if ( !constraintsNode.isNull() )
    {
      QDomNodeList constraintNodeList = constraintsNode.toElement().elementsByTagName( QStringLiteral( "constraint" ) );
      for ( int i = 0; i < constraintNodeList.size(); ++i )
      {
        QDomElement constraintElem = constraintNodeList.at( i ).toElement();

        QString field = constraintElem.attribute( QStringLiteral( "field" ), QString() );
        int constraints = constraintElem.attribute( QStringLiteral( "constraints" ), QStringLiteral( "0" ) ).toInt();
        if ( field.isEmpty() || constraints == 0 )
          continue;

        mFieldConstraints.insert( field, static_cast< QgsFieldConstraints::Constraints >( constraints ) );

        int uniqueStrength = constraintElem.attribute( QStringLiteral( "unique_strength" ), QStringLiteral( "1" ) ).toInt();
        int notNullStrength = constraintElem.attribute( QStringLiteral( "notnull_strength" ), QStringLiteral( "1" ) ).toInt();
        int expStrength = constraintElem.attribute( QStringLiteral( "exp_strength" ), QStringLiteral( "1" ) ).toInt();

        mFieldConstraintStrength.insert( qMakePair( field, QgsFieldConstraints::ConstraintUnique ), static_cast< QgsFieldConstraints::ConstraintStrength >( uniqueStrength ) );
        mFieldConstraintStrength.insert( qMakePair( field, QgsFieldConstraints::ConstraintNotNull ), static_cast< QgsFieldConstraints::ConstraintStrength >( notNullStrength ) );
        mFieldConstraintStrength.insert( qMakePair( field, QgsFieldConstraints::ConstraintExpression ), static_cast< QgsFieldConstraints::ConstraintStrength >( expStrength ) );
      }
    }
    mFieldConstraintExpressions.clear();
    QDomNode constraintExpressionsNode = layerNode.namedItem( QStringLiteral( "constraintExpressions" ) );
    if ( !constraintExpressionsNode.isNull() )
    {
      QDomNodeList constraintNodeList = constraintExpressionsNode.toElement().elementsByTagName( QStringLiteral( "constraint" ) );
      for ( int i = 0; i < constraintNodeList.size(); ++i )
      {
        QDomElement constraintElem = constraintNodeList.at( i ).toElement();

        QString field = constraintElem.attribute( QStringLiteral( "field" ), QString() );
        QString exp = constraintElem.attribute( QStringLiteral( "exp" ), QString() );
        QString desc = constraintElem.attribute( QStringLiteral( "desc" ), QString() );
        if ( field.isEmpty() || exp.isEmpty() )
          continue;

        mFieldConstraintExpressions.insert( field, qMakePair( exp, desc ) );
      }
    }

    updateFields();

    //Attributes excluded from WMS and WFS
    mExcludeAttributesWMS.clear();
    QDomNode excludeWMSNode = layerNode.namedItem( QStringLiteral( "excludeAttributesWMS" ) );
    if ( !excludeWMSNode.isNull() )
    {
      QDomNodeList attributeNodeList = excludeWMSNode.toElement().elementsByTagName( QStringLiteral( "attribute" ) );
      for ( int i = 0; i < attributeNodeList.size(); ++i )
      {
        mExcludeAttributesWMS.insert( attributeNodeList.at( i ).toElement().text() );
      }
    }

    mExcludeAttributesWFS.clear();
    QDomNode excludeWFSNode = layerNode.namedItem( QStringLiteral( "excludeAttributesWFS" ) );
    if ( !excludeWFSNode.isNull() )
    {
      QDomNodeList attributeNodeList = excludeWFSNode.toElement().elementsByTagName( QStringLiteral( "attribute" ) );
      for ( int i = 0; i < attributeNodeList.size(); ++i )
      {
        mExcludeAttributesWFS.insert( attributeNodeList.at( i ).toElement().text() );
      }
    }

    // Load editor widget configuration
    QDomElement widgetsElem = layerNode.namedItem( QStringLiteral( "fieldConfiguration" ) ).toElement();
    QDomNodeList fieldConfigurationElementList = widgetsElem.elementsByTagName( QStringLiteral( "field" ) );
    for ( int i = 0; i < fieldConfigurationElementList.size(); ++i )
    {
      const QDomElement fieldConfigElement = fieldConfigurationElementList.at( i ).toElement();
      const QDomElement fieldWidgetElement = fieldConfigElement.elementsByTagName( QStringLiteral( "editWidget" ) ).at( 0 ).toElement();

      QString fieldName = fieldConfigElement.attribute( QStringLiteral( "name" ) );

      const QString widgetType = fieldWidgetElement.attribute( QStringLiteral( "type" ) );
      const QDomElement cfgElem = fieldConfigElement.elementsByTagName( QStringLiteral( "config" ) ).at( 0 ).toElement();
      const QDomElement optionsElem = cfgElem.childNodes().at( 0 ).toElement();
      QVariantMap optionsMap = QgsXmlUtils::readVariant( optionsElem ).toMap();
      if ( widgetType == QStringLiteral( "ValueRelation" ) )
      {
        optionsMap[ QStringLiteral( "Value" ) ] = context.projectTranslator()->translate( QStringLiteral( "project:layers:%1:fields:%2:valuerelationvalue" ).arg( layerNode.namedItem( QStringLiteral( "id" ) ).toElement().text(), fieldName ), optionsMap[ QStringLiteral( "Value" ) ].toString() );
      }
      QgsEditorWidgetSetup setup = QgsEditorWidgetSetup( widgetType, optionsMap );
      mFieldWidgetSetups[fieldName] = setup;
    }
  }

  if ( categories.testFlag( GeometryOptions ) )
    mGeometryOptions->readXml( layerNode.namedItem( QStringLiteral( "geometryOptions" ) ) );

  if ( categories.testFlag( Forms ) )
    mEditFormConfig.readXml( layerNode, context );

  if ( categories.testFlag( AttributeTable ) )
  {
    mAttributeTableConfig.readXml( layerNode );
    mConditionalStyles->readXml( layerNode, context );
  }

  if ( categories.testFlag( CustomProperties ) )
    readCustomProperties( layerNode, QStringLiteral( "variable" ) );

  QDomElement mapLayerNode = layerNode.toElement();
  if ( categories.testFlag( LayerConfiguration )
       && mapLayerNode.attribute( QStringLiteral( "readOnly" ), QStringLiteral( "0" ) ).toInt() == 1 )
    mReadOnly = true;

  updateFields();

  return true;
}

bool QgsVectorLayer::readStyle( const QDomNode &node, QString &errorMessage,
                                QgsReadWriteContext &context, QgsMapLayer::StyleCategories categories )
{
  bool result = true;
  emit readCustomSymbology( node.toElement(), errorMessage );

  if ( isSpatial() )
  {
    // try renderer v2 first
    if ( categories.testFlag( Symbology ) )
    {
      QDomElement rendererElement = node.firstChildElement( RENDERER_TAG_NAME );
      if ( !rendererElement.isNull() )
      {
        QgsFeatureRenderer *r = QgsFeatureRenderer::load( rendererElement, context );
        if ( r )
        {
          setRenderer( r );
        }
        else
        {
          result = false;
        }
      }
      // make sure layer has a renderer - if none exists, fallback to a default renderer
      if ( !renderer() )
      {
        setRenderer( QgsFeatureRenderer::defaultRenderer( geometryType() ) );
      }
    }

    // read labeling definition
    if ( categories.testFlag( Labeling ) )
    {
      QDomElement labelingElement = node.firstChildElement( QStringLiteral( "labeling" ) );
      QgsAbstractVectorLayerLabeling *labeling = nullptr;
      if ( labelingElement.isNull() ||
           ( labelingElement.attribute( QStringLiteral( "type" ) ) == QLatin1String( "simple" ) && labelingElement.firstChildElement( QStringLiteral( "settings" ) ).isNull() ) )
      {
        // make sure we have custom properties for labeling for 2.x projects
        // (custom properties should be already loaded when reading the whole layer from XML,
        // but when reading style, custom properties are not read)
        readCustomProperties( node, QStringLiteral( "labeling" ) );

        // support for pre-QGIS 3 labeling configurations written in custom properties
        labeling = readLabelingFromCustomProperties();
      }
      else
      {
        labeling = QgsAbstractVectorLayerLabeling::create( labelingElement, context );
      }
      setLabeling( labeling );

      if ( node.toElement().hasAttribute( QStringLiteral( "labelsEnabled" ) ) )
        mLabelsEnabled = node.toElement().attribute( QStringLiteral( "labelsEnabled" ) ).toInt();
      else
        mLabelsEnabled = true;
    }

    if ( categories.testFlag( Symbology ) )
    {
      // get and set the blend mode if it exists
      QDomNode blendModeNode = node.namedItem( QStringLiteral( "blendMode" ) );
      if ( !blendModeNode.isNull() )
      {
        QDomElement e = blendModeNode.toElement();
        setBlendMode( QgsPainting::getCompositionMode( static_cast< QgsPainting::BlendMode >( e.text().toInt() ) ) );
      }

      // get and set the feature blend mode if it exists
      QDomNode featureBlendModeNode = node.namedItem( QStringLiteral( "featureBlendMode" ) );
      if ( !featureBlendModeNode.isNull() )
      {
        QDomElement e = featureBlendModeNode.toElement();
        setFeatureBlendMode( QgsPainting::getCompositionMode( static_cast< QgsPainting::BlendMode >( e.text().toInt() ) ) );
      }
    }

    // get and set the layer transparency if it exists
    if ( categories.testFlag( Rendering ) )
    {
      QDomNode layerTransparencyNode = node.namedItem( QStringLiteral( "layerTransparency" ) );
      if ( !layerTransparencyNode.isNull() )
      {
        QDomElement e = layerTransparencyNode.toElement();
        setOpacity( 1.0 - e.text().toInt() / 100.0 );
      }
      QDomNode layerOpacityNode = node.namedItem( QStringLiteral( "layerOpacity" ) );
      if ( !layerOpacityNode.isNull() )
      {
        QDomElement e = layerOpacityNode.toElement();
        setOpacity( e.text().toDouble() );
      }
    }

    if ( categories.testFlag( Rendering ) )
    {
      QDomElement e = node.toElement();

      // get the simplification drawing settings
      mSimplifyMethod.setSimplifyHints( static_cast< QgsVectorSimplifyMethod::SimplifyHints >( e.attribute( QStringLiteral( "simplifyDrawingHints" ), QStringLiteral( "1" ) ).toInt() ) );
      mSimplifyMethod.setSimplifyAlgorithm( static_cast< QgsVectorSimplifyMethod::SimplifyAlgorithm >( e.attribute( QStringLiteral( "simplifyAlgorithm" ), QStringLiteral( "0" ) ).toInt() ) );
      mSimplifyMethod.setThreshold( e.attribute( QStringLiteral( "simplifyDrawingTol" ), QStringLiteral( "1" ) ).toFloat() );
      mSimplifyMethod.setForceLocalOptimization( e.attribute( QStringLiteral( "simplifyLocal" ), QStringLiteral( "1" ) ).toInt() );
      mSimplifyMethod.setMaximumScale( e.attribute( QStringLiteral( "simplifyMaxScale" ), QStringLiteral( "1" ) ).toFloat() );
    }

    //diagram renderer and diagram layer settings
    if ( categories.testFlag( Diagrams ) )
    {
      delete mDiagramRenderer;
      mDiagramRenderer = nullptr;
      QDomElement singleCatDiagramElem = node.firstChildElement( QStringLiteral( "SingleCategoryDiagramRenderer" ) );
      if ( !singleCatDiagramElem.isNull() )
      {
        mDiagramRenderer = new QgsSingleCategoryDiagramRenderer();
        mDiagramRenderer->readXml( singleCatDiagramElem, context );
      }
      QDomElement linearDiagramElem = node.firstChildElement( QStringLiteral( "LinearlyInterpolatedDiagramRenderer" ) );
      if ( !linearDiagramElem.isNull() )
      {
        if ( linearDiagramElem.hasAttribute( QStringLiteral( "classificationAttribute" ) ) )
        {
          // fix project from before QGIS 3.0
          int idx = linearDiagramElem.attribute( QStringLiteral( "classificationAttribute" ) ).toInt();
          if ( idx >= 0 && idx < mFields.count() )
            linearDiagramElem.setAttribute( QStringLiteral( "classificationField" ), mFields.at( idx ).name() );
        }

        mDiagramRenderer = new QgsLinearlyInterpolatedDiagramRenderer();
        mDiagramRenderer->readXml( linearDiagramElem, context );
      }

      if ( mDiagramRenderer )
      {
        QDomElement diagramSettingsElem = node.firstChildElement( QStringLiteral( "DiagramLayerSettings" ) );
        if ( !diagramSettingsElem.isNull() )
        {
          bool oldXPos = diagramSettingsElem.hasAttribute( QStringLiteral( "xPosColumn" ) );
          bool oldYPos = diagramSettingsElem.hasAttribute( QStringLiteral( "yPosColumn" ) );
          bool oldShow = diagramSettingsElem.hasAttribute( QStringLiteral( "showColumn" ) );
          if ( oldXPos || oldYPos || oldShow )
          {
            // fix project from before QGIS 3.0
            QgsPropertyCollection ddp;
            if ( oldXPos )
            {
              int xPosColumn = diagramSettingsElem.attribute( QStringLiteral( "xPosColumn" ) ).toInt();
              if ( xPosColumn >= 0 && xPosColumn < mFields.count() )
                ddp.setProperty( QgsDiagramLayerSettings::PositionX, QgsProperty::fromField( mFields.at( xPosColumn ).name(), true ) );
            }
            if ( oldYPos )
            {
              int yPosColumn = diagramSettingsElem.attribute( QStringLiteral( "yPosColumn" ) ).toInt();
              if ( yPosColumn >= 0 && yPosColumn < mFields.count() )
                ddp.setProperty( QgsDiagramLayerSettings::PositionY, QgsProperty::fromField( mFields.at( yPosColumn ).name(), true ) );
            }
            if ( oldShow )
            {
              int showColumn = diagramSettingsElem.attribute( QStringLiteral( "showColumn" ) ).toInt();
              if ( showColumn >= 0 && showColumn < mFields.count() )
                ddp.setProperty( QgsDiagramLayerSettings::Show, QgsProperty::fromField( mFields.at( showColumn ).name(), true ) );
            }
            QDomElement propertiesElem = diagramSettingsElem.ownerDocument().createElement( QStringLiteral( "properties" ) );
            QgsPropertiesDefinition defs = QgsPropertiesDefinition
            {
              { QgsDiagramLayerSettings::PositionX, QgsPropertyDefinition( "positionX", QObject::tr( "Position (X)" ), QgsPropertyDefinition::Double ) },
              { QgsDiagramLayerSettings::PositionY, QgsPropertyDefinition( "positionY", QObject::tr( "Position (Y)" ), QgsPropertyDefinition::Double ) },
              { QgsDiagramLayerSettings::Show, QgsPropertyDefinition( "show", QObject::tr( "Show diagram" ), QgsPropertyDefinition::Boolean ) },
            };
            ddp.writeXml( propertiesElem, defs );
            diagramSettingsElem.appendChild( propertiesElem );
          }

          delete mDiagramLayerSettings;
          mDiagramLayerSettings = new QgsDiagramLayerSettings();
          mDiagramLayerSettings->readXml( diagramSettingsElem );
        }
      }
    }
    // end diagram
  }
  return result;
}


bool QgsVectorLayer::writeSymbology( QDomNode &node, QDomDocument &doc, QString &errorMessage,
                                     const QgsReadWriteContext &context, QgsMapLayer::StyleCategories categories ) const
{
  QDomElement layerElement = node.toElement();
  writeCommonStyle( layerElement, doc, context, categories );

  ( void )writeStyle( node, doc, errorMessage, context, categories );

  if ( categories.testFlag( GeometryOptions ) )
    mGeometryOptions->writeXml( node );

  if ( categories.testFlag( Fields ) )
  {
    QDomElement fieldConfigurationElement = doc.createElement( QStringLiteral( "fieldConfiguration" ) );
    node.appendChild( fieldConfigurationElement );

    int index = 0;
    for ( const QgsField &field : mFields )
    {
      QDomElement fieldElement = doc.createElement( QStringLiteral( "field" ) );
      fieldElement.setAttribute( QStringLiteral( "name" ), field.name() );

      fieldConfigurationElement.appendChild( fieldElement );

      QgsEditorWidgetSetup widgetSetup = field.editorWidgetSetup();

      // TODO : wrap this part in an if to only save if it was user-modified
      QDomElement editWidgetElement = doc.createElement( QStringLiteral( "editWidget" ) );
      fieldElement.appendChild( editWidgetElement );
      editWidgetElement.setAttribute( QStringLiteral( "type" ), field.editorWidgetSetup().type() );
      QDomElement editWidgetConfigElement = doc.createElement( QStringLiteral( "config" ) );

      editWidgetConfigElement.appendChild( QgsXmlUtils::writeVariant( widgetSetup.config(), doc ) );
      editWidgetElement.appendChild( editWidgetConfigElement );
      // END TODO : wrap this part in an if to only save if it was user-modified

      ++index;
    }

    //attribute aliases
    QDomElement aliasElem = doc.createElement( QStringLiteral( "aliases" ) );
    for ( const QgsField &field : mFields )
    {
      QDomElement aliasEntryElem = doc.createElement( QStringLiteral( "alias" ) );
      aliasEntryElem.setAttribute( QStringLiteral( "field" ), field.name() );
      aliasEntryElem.setAttribute( QStringLiteral( "index" ), mFields.indexFromName( field.name() ) );
      aliasEntryElem.setAttribute( QStringLiteral( "name" ), field.alias() );
      aliasElem.appendChild( aliasEntryElem );
    }
    node.appendChild( aliasElem );

    //exclude attributes WMS
    QDomElement excludeWMSElem = doc.createElement( QStringLiteral( "excludeAttributesWMS" ) );
    QSet<QString>::const_iterator attWMSIt = mExcludeAttributesWMS.constBegin();
    for ( ; attWMSIt != mExcludeAttributesWMS.constEnd(); ++attWMSIt )
    {
      QDomElement attrElem = doc.createElement( QStringLiteral( "attribute" ) );
      QDomText attrText = doc.createTextNode( *attWMSIt );
      attrElem.appendChild( attrText );
      excludeWMSElem.appendChild( attrElem );
    }
    node.appendChild( excludeWMSElem );

    //exclude attributes WFS
    QDomElement excludeWFSElem = doc.createElement( QStringLiteral( "excludeAttributesWFS" ) );
    QSet<QString>::const_iterator attWFSIt = mExcludeAttributesWFS.constBegin();
    for ( ; attWFSIt != mExcludeAttributesWFS.constEnd(); ++attWFSIt )
    {
      QDomElement attrElem = doc.createElement( QStringLiteral( "attribute" ) );
      QDomText attrText = doc.createTextNode( *attWFSIt );
      attrElem.appendChild( attrText );
      excludeWFSElem.appendChild( attrElem );
    }
    node.appendChild( excludeWFSElem );

    //default expressions
    QDomElement defaultsElem = doc.createElement( QStringLiteral( "defaults" ) );
    for ( const QgsField &field : mFields )
    {
      QDomElement defaultElem = doc.createElement( QStringLiteral( "default" ) );
      defaultElem.setAttribute( QStringLiteral( "field" ), field.name() );
      defaultElem.setAttribute( QStringLiteral( "expression" ), field.defaultValueDefinition().expression() );
      defaultElem.setAttribute( QStringLiteral( "applyOnUpdate" ), field.defaultValueDefinition().applyOnUpdate() ? QStringLiteral( "1" ) : QStringLiteral( "0" ) );
      defaultsElem.appendChild( defaultElem );
    }
    node.appendChild( defaultsElem );

    // constraints
    QDomElement constraintsElem = doc.createElement( QStringLiteral( "constraints" ) );
    for ( const QgsField &field : mFields )
    {
      QDomElement constraintElem = doc.createElement( QStringLiteral( "constraint" ) );
      constraintElem.setAttribute( QStringLiteral( "field" ), field.name() );
      constraintElem.setAttribute( QStringLiteral( "constraints" ), field.constraints().constraints() );
      constraintElem.setAttribute( QStringLiteral( "unique_strength" ), field.constraints().constraintStrength( QgsFieldConstraints::ConstraintUnique ) );
      constraintElem.setAttribute( QStringLiteral( "notnull_strength" ), field.constraints().constraintStrength( QgsFieldConstraints::ConstraintNotNull ) );
      constraintElem.setAttribute( QStringLiteral( "exp_strength" ), field.constraints().constraintStrength( QgsFieldConstraints::ConstraintExpression ) );
      constraintsElem.appendChild( constraintElem );
    }
    node.appendChild( constraintsElem );

    // constraint expressions
    QDomElement constraintExpressionsElem = doc.createElement( QStringLiteral( "constraintExpressions" ) );
    for ( const QgsField &field : mFields )
    {
      QDomElement constraintExpressionElem = doc.createElement( QStringLiteral( "constraint" ) );
      constraintExpressionElem.setAttribute( QStringLiteral( "field" ), field.name() );
      constraintExpressionElem.setAttribute( QStringLiteral( "exp" ), field.constraints().constraintExpression() );
      constraintExpressionElem.setAttribute( QStringLiteral( "desc" ), field.constraints().constraintDescription() );
      constraintExpressionsElem.appendChild( constraintExpressionElem );
    }
    node.appendChild( constraintExpressionsElem );

    // save expression fields
    if ( !mExpressionFieldBuffer )
    {
      // can happen when saving style on a invalid layer
      QgsExpressionFieldBuffer dummy;
      dummy.writeXml( node, doc );
    }
    else
    {
      mExpressionFieldBuffer->writeXml( node, doc );
    }
  }

  // add attribute actions
  if ( categories.testFlag( Actions ) )
    mActions->writeXml( node );

  if ( categories.testFlag( AttributeTable ) )
  {
    mAttributeTableConfig.writeXml( node );
    mConditionalStyles->writeXml( node, doc, context );
  }

  if ( categories.testFlag( Forms ) )
    mEditFormConfig.writeXml( node, context );

  // save readonly state
  if ( categories.testFlag( LayerConfiguration ) )
    node.toElement().setAttribute( QStringLiteral( "readOnly" ), mReadOnly );

  // save preview expression
  if ( categories.testFlag( LayerConfiguration ) )
  {
    QDomElement prevExpElem = doc.createElement( QStringLiteral( "previewExpression" ) );
    QDomText prevExpText = doc.createTextNode( mDisplayExpression );
    prevExpElem.appendChild( prevExpText );
    node.appendChild( prevExpElem );
  }

  // save map tip
  if ( categories.testFlag( MapTips ) )
  {
    QDomElement mapTipElem = doc.createElement( QStringLiteral( "mapTip" ) );
    QDomText mapTipText = doc.createTextNode( mMapTipTemplate );
    mapTipElem.appendChild( mapTipText );
    node.toElement().appendChild( mapTipElem );
  }

  return true;
}

bool QgsVectorLayer::writeStyle( QDomNode &node, QDomDocument &doc, QString &errorMessage,
                                 const QgsReadWriteContext &context, QgsMapLayer::StyleCategories categories ) const
{
  QDomElement mapLayerNode = node.toElement();

  emit writeCustomSymbology( mapLayerNode, doc, errorMessage );

  if ( isSpatial() )
  {
    if ( categories.testFlag( Symbology ) )
    {
      if ( mRenderer )
      {
        QDomElement rendererElement = mRenderer->save( doc, context );
        node.appendChild( rendererElement );
      }
    }

    if ( categories.testFlag( Labeling ) )
    {
      if ( mLabeling )
      {
        QDomElement labelingElement = mLabeling->save( doc, context );
        node.appendChild( labelingElement );
      }
      mapLayerNode.setAttribute( QStringLiteral( "labelsEnabled" ), mLabelsEnabled ? QStringLiteral( "1" ) : QStringLiteral( "0" ) );
    }

    // save the simplification drawing settings
    if ( categories.testFlag( Rendering ) )
    {
      mapLayerNode.setAttribute( QStringLiteral( "simplifyDrawingHints" ), QString::number( mSimplifyMethod.simplifyHints() ) );
      mapLayerNode.setAttribute( QStringLiteral( "simplifyAlgorithm" ), QString::number( mSimplifyMethod.simplifyAlgorithm() ) );
      mapLayerNode.setAttribute( QStringLiteral( "simplifyDrawingTol" ), QString::number( mSimplifyMethod.threshold() ) );
      mapLayerNode.setAttribute( QStringLiteral( "simplifyLocal" ), mSimplifyMethod.forceLocalOptimization() ? 1 : 0 );
      mapLayerNode.setAttribute( QStringLiteral( "simplifyMaxScale" ), QString::number( mSimplifyMethod.maximumScale() ) );
    }

    //save customproperties
    if ( categories.testFlag( CustomProperties ) )
    {
      writeCustomProperties( node, doc );
    }

    if ( categories.testFlag( Symbology ) )
    {
      // add the blend mode field
      QDomElement blendModeElem  = doc.createElement( QStringLiteral( "blendMode" ) );
      QDomText blendModeText = doc.createTextNode( QString::number( QgsPainting::getBlendModeEnum( blendMode() ) ) );
      blendModeElem.appendChild( blendModeText );
      node.appendChild( blendModeElem );

      // add the feature blend mode field
      QDomElement featureBlendModeElem  = doc.createElement( QStringLiteral( "featureBlendMode" ) );
      QDomText featureBlendModeText = doc.createTextNode( QString::number( QgsPainting::getBlendModeEnum( featureBlendMode() ) ) );
      featureBlendModeElem.appendChild( featureBlendModeText );
      node.appendChild( featureBlendModeElem );
    }

    // add the layer opacity
    if ( categories.testFlag( Rendering ) )
    {
      QDomElement layerOpacityElem  = doc.createElement( QStringLiteral( "layerOpacity" ) );
      QDomText layerOpacityText = doc.createTextNode( QString::number( opacity() ) );
      layerOpacityElem.appendChild( layerOpacityText );
      node.appendChild( layerOpacityElem );
    }

    if ( categories.testFlag( Diagrams ) && mDiagramRenderer )
    {
      mDiagramRenderer->writeXml( mapLayerNode, doc, context );
      if ( mDiagramLayerSettings )
        mDiagramLayerSettings->writeXml( mapLayerNode, doc );
    }
  }
  return true;
}

bool QgsVectorLayer::readSld( const QDomNode &node, QString &errorMessage )
{
  // get the Name element
  QDomElement nameElem = node.firstChildElement( QStringLiteral( "Name" ) );
  if ( nameElem.isNull() )
  {
    errorMessage = QStringLiteral( "Warning: Name element not found within NamedLayer while it's required." );
  }

  if ( isSpatial() )
  {
    QgsFeatureRenderer *r = QgsFeatureRenderer::loadSld( node, geometryType(), errorMessage );
    if ( !r )
      return false;

    setRenderer( r );

    // labeling
    readSldLabeling( node );
  }
  return true;
}

bool QgsVectorLayer::writeSld( QDomNode &node, QDomDocument &doc, QString &errorMessage, const QgsStringMap &props ) const
{
  Q_UNUSED( errorMessage );

  QgsStringMap localProps = QgsStringMap( props );
  if ( hasScaleBasedVisibility() )
  {
    QgsSymbolLayerUtils::mergeScaleDependencies( maximumScale(), minimumScale(), localProps );
  }

  if ( isSpatial() )
  {
    // store the Name element
    QDomElement nameNode = doc.createElement( QStringLiteral( "se:Name" ) );
    nameNode.appendChild( doc.createTextNode( name() ) );
    node.appendChild( nameNode );

    QDomElement userStyleElem = doc.createElement( QStringLiteral( "UserStyle" ) );
    node.appendChild( userStyleElem );

    QDomElement nameElem = doc.createElement( QStringLiteral( "se:Name" ) );
    nameElem.appendChild( doc.createTextNode( name() ) );

    userStyleElem.appendChild( nameElem );

    QDomElement featureTypeStyleElem = doc.createElement( QStringLiteral( "se:FeatureTypeStyle" ) );
    userStyleElem.appendChild( featureTypeStyleElem );

    mRenderer->toSld( doc, featureTypeStyleElem, localProps );
    if ( labelsEnabled() )
    {
      mLabeling->toSld( featureTypeStyleElem, localProps );
    }
  }
  return true;
}


bool QgsVectorLayer::changeGeometry( QgsFeatureId fid, QgsGeometry &geom, bool skipDefaultValue )
{
  if ( !mEditBuffer || !mDataProvider )
  {
    return false;
  }

  if ( mGeometryOptions->isActive() )
    mGeometryOptions->apply( geom );

  updateExtents();

  bool result = mEditBuffer->changeGeometry( fid, geom );

  if ( result )
  {
    updateExtents();
    if ( !skipDefaultValue && !mDefaultValueOnUpdateFields.isEmpty() )
      updateDefaultValues( fid );
  }
  return result;
}


bool QgsVectorLayer::changeAttributeValue( QgsFeatureId fid, int field, const QVariant &newValue, const QVariant &oldValue, bool skipDefaultValues )
{
  bool result = false;

  switch ( fields().fieldOrigin( field ) )
  {
    case QgsFields::OriginJoin:
      result = mJoinBuffer->changeAttributeValue( fid, field, newValue, oldValue );
      break;

    case QgsFields::OriginProvider:
    case QgsFields::OriginEdit:
    case QgsFields::OriginExpression:
    {
      if ( mEditBuffer && mDataProvider )
        result = mEditBuffer->changeAttributeValue( fid, field, newValue, oldValue );
      break;
    }

    case QgsFields::OriginUnknown:
      break;
  }

  if ( result && !skipDefaultValues && !mDefaultValueOnUpdateFields.isEmpty() )
    updateDefaultValues( fid );

  return result;
}

bool QgsVectorLayer::changeAttributeValues( QgsFeatureId fid, const QgsAttributeMap &newValues, const QgsAttributeMap &oldValues, bool skipDefaultValues )
{
  bool result = true;

  QgsAttributeMap newValuesJoin;
  QgsAttributeMap oldValuesJoin;

  QgsAttributeMap newValuesNotJoin;
  QgsAttributeMap oldValuesNotJoin;

  for ( auto it = newValues.constBegin(); it != newValues.constEnd(); ++it )
  {
    const int field = it.key();
    const QVariant newValue = it.value();
    QVariant oldValue;

    if ( oldValues.contains( field ) )
      oldValue = oldValues[field];

    switch ( fields().fieldOrigin( field ) )
    {
      case QgsFields::OriginJoin:
        newValuesJoin[field] = newValue;
        oldValuesJoin[field] = oldValue;
        break;

      case QgsFields::OriginProvider:
      case QgsFields::OriginEdit:
      case QgsFields::OriginExpression:
      {
        newValuesNotJoin[field] = newValue;
        oldValuesNotJoin[field] = oldValue;
        break;
      }

      case QgsFields::OriginUnknown:
        break;
    }
  }

  if ( ! newValuesJoin.isEmpty() && mJoinBuffer )
  {
    result = mJoinBuffer->changeAttributeValues( fid, newValuesJoin, oldValuesJoin );
  }

  if ( ! newValuesNotJoin.isEmpty() && mEditBuffer && mDataProvider )
  {
    result &= mEditBuffer->changeAttributeValues( fid, newValues, oldValues );
  }

  if ( result && !skipDefaultValues && !mDefaultValueOnUpdateFields.isEmpty() )
  {
    updateDefaultValues( fid );
  }

  return result;
}

bool QgsVectorLayer::addAttribute( const QgsField &field )
{
  if ( !mEditBuffer || !mDataProvider )
    return false;

  return mEditBuffer->addAttribute( field );
}

void QgsVectorLayer::removeFieldAlias( int attIndex )
{
  if ( attIndex < 0 || attIndex >= fields().count() )
    return;

  QString name = fields().at( attIndex ).name();
  mFields[ attIndex ].setAlias( QString() );
  if ( mAttributeAliasMap.contains( name ) )
  {
    mAttributeAliasMap.remove( name );
    updateFields();
    mEditFormConfig.setFields( mFields );
    emit layerModified();
  }
}

bool QgsVectorLayer::renameAttribute( int index, const QString &newName )
{
  if ( index < 0 || index >= fields().count() )
    return false;

  switch ( mFields.fieldOrigin( index ) )
  {
    case QgsFields::OriginExpression:
    {
      if ( mExpressionFieldBuffer )
      {
        int oi = mFields.fieldOriginIndex( index );
        mExpressionFieldBuffer->renameExpression( oi, newName );
        updateFields();
        return true;
      }
      else
      {
        return false;
      }
    }

    case QgsFields::OriginProvider:
    case QgsFields::OriginEdit:

      if ( !mEditBuffer || !mDataProvider )
        return false;

      return mEditBuffer->renameAttribute( index, newName );

    case QgsFields::OriginJoin:
    case QgsFields::OriginUnknown:
      return false;

  }

  return false; // avoid warning
}

void QgsVectorLayer::setFieldAlias( int attIndex, const QString &aliasString )
{
  if ( attIndex < 0 || attIndex >= fields().count() )
    return;

  QString name = fields().at( attIndex ).name();

  mAttributeAliasMap.insert( name, aliasString );
  mFields[ attIndex ].setAlias( aliasString );
  mEditFormConfig.setFields( mFields );
  emit layerModified(); // TODO[MD]: should have a different signal?
}

QString QgsVectorLayer::attributeAlias( int index ) const
{
  if ( index < 0 || index >= fields().count() )
    return QString();

  return fields().at( index ).alias();
}

QString QgsVectorLayer::attributeDisplayName( int index ) const
{
  if ( index >= 0 && index < mFields.count() )
    return mFields.at( index ).displayName();
  else
    return QString();
}

QgsStringMap QgsVectorLayer::attributeAliases() const
{
  return mAttributeAliasMap;
}

bool QgsVectorLayer::deleteAttribute( int index )
{
  if ( index < 0 || index >= fields().count() )
    return false;

  if ( mFields.fieldOrigin( index ) == QgsFields::OriginExpression )
  {
    removeExpressionField( index );
    return true;
  }

  if ( !mEditBuffer || !mDataProvider )
    return false;

  return mEditBuffer->deleteAttribute( index );
}

bool QgsVectorLayer::deleteAttributes( const QList<int> &attrs )
{
  bool deleted = false;

  // Remove multiple occurrences of same attribute
  QList<int> attrList = attrs.toSet().toList();

  std::sort( attrList.begin(), attrList.end(), std::greater<int>() );

  for ( int attr : qgis::as_const( attrList ) )
  {
    if ( deleteAttribute( attr ) )
    {
      deleted = true;
    }
  }

  return deleted;
}

bool QgsVectorLayer::deleteFeature( QgsFeatureId fid )
{
  if ( !mEditBuffer )
    return false;

  if ( mJoinBuffer->containsJoins() )
    mJoinBuffer->deleteFeature( fid );

  bool res = mEditBuffer->deleteFeature( fid );
  if ( res )
  {
    mSelectedFeatureIds.remove( fid ); // remove it from selection
    updateExtents();
  }

  return res;
}

bool QgsVectorLayer::deleteFeatures( const QgsFeatureIds &fids )
{
  if ( !mEditBuffer )
  {
    QgsDebugMsgLevel( QStringLiteral( "Cannot delete features (mEditBuffer==NULL)" ), 1 );
    return false;
  }

  if ( mJoinBuffer->containsJoins() )
    mJoinBuffer->deleteFeatures( fids );

  bool res = mEditBuffer->deleteFeatures( fids );

  if ( res )
  {
    mSelectedFeatureIds.subtract( fids ); // remove it from selection
    updateExtents();
  }

  return res;
}

QgsFields QgsVectorLayer::fields() const
{
  return mFields;
}

QgsAttributeList QgsVectorLayer::primaryKeyAttributes() const
{
  QgsAttributeList pkAttributesList;

  QgsAttributeList providerIndexes = mDataProvider->pkAttributeIndexes();
  for ( int i = 0; i < mFields.count(); ++i )
  {
    if ( mFields.fieldOrigin( i ) == QgsFields::OriginProvider &&
         providerIndexes.contains( mFields.fieldOriginIndex( i ) ) )
      pkAttributesList << i;
  }

  return pkAttributesList;
}

long QgsVectorLayer::featureCount() const
{
  if ( ! mDataProvider )
    return -1;
  return mDataProvider->featureCount() +
         ( mEditBuffer ? mEditBuffer->mAddedFeatures.size() - mEditBuffer->mDeletedFeatureIds.size() : 0 );
}

QgsFeatureSource::FeatureAvailability QgsVectorLayer::hasFeatures() const
{
  const QgsFeatureIds deletedFeatures( mEditBuffer ? mEditBuffer->deletedFeatureIds() : QgsFeatureIds() );
  const QgsFeatureMap addedFeatures( mEditBuffer ? mEditBuffer->addedFeatures() : QgsFeatureMap() );

  if ( mEditBuffer && !deletedFeatures.empty() )
  {
    if ( addedFeatures.size() > deletedFeatures.size() )
      return QgsFeatureSource::FeatureAvailability::FeaturesAvailable;
    else
      return QgsFeatureSource::FeatureAvailability::FeaturesMaybeAvailable;
  }

  if ( ( !mEditBuffer || addedFeatures.empty() ) && mDataProvider->empty() )
    return QgsFeatureSource::FeatureAvailability::NoFeaturesAvailable;
  else
    return QgsFeatureSource::FeatureAvailability::FeaturesAvailable;
}

bool QgsVectorLayer::commitChanges()
{
  mCommitErrors.clear();

  if ( !mDataProvider )
  {
    mCommitErrors << tr( "ERROR: no provider" );
    return false;
  }

  if ( !mEditBuffer )
  {
    mCommitErrors << tr( "ERROR: layer not editable" );
    return false;
  }

  emit beforeCommitChanges();

  if ( !mAllowCommit )
    return false;

  bool success = mEditBuffer->commitChanges( mCommitErrors );

  if ( success )
  {
    delete mEditBuffer;
    mEditBuffer = nullptr;
    undoStack()->clear();
    emit editingStopped();
  }
  else
  {
    QgsMessageLog::logMessage( tr( "Commit errors:\n  %1" ).arg( mCommitErrors.join( QStringLiteral( "\n  " ) ) ) );
  }

  updateFields();
  mDataProvider->updateExtents();

  mDataProvider->leaveUpdateMode();

  emit repaintRequested();

  return success;
}

QStringList QgsVectorLayer::commitErrors() const
{
  return mCommitErrors;
}

bool QgsVectorLayer::rollBack( bool deleteBuffer )
{
  if ( !mEditBuffer )
  {
    return false;
  }

  bool rollbackExtent = !mEditBuffer->mDeletedFeatureIds.isEmpty() ||
                        !mEditBuffer->mAddedFeatures.isEmpty() ||
                        !mEditBuffer->mChangedGeometries.isEmpty();

  emit beforeRollBack();

  mEditBuffer->rollBack();

  emit afterRollBack();

  if ( isModified() )
  {
    // new undo stack roll back method
    // old method of calling every undo could cause many canvas refreshes
    undoStack()->setIndex( 0 );
  }

  updateFields();

  if ( deleteBuffer )
  {
    delete mEditBuffer;
    mEditBuffer = nullptr;
    undoStack()->clear();
  }
  emit editingStopped();

  if ( rollbackExtent )
    updateExtents();

  mDataProvider->leaveUpdateMode();

  emit repaintRequested();
  return true;
}

int QgsVectorLayer::selectedFeatureCount() const
{
  return mSelectedFeatureIds.size();
}

const QgsFeatureIds &QgsVectorLayer::selectedFeatureIds() const
{
  return mSelectedFeatureIds;
}

QgsFeatureList QgsVectorLayer::selectedFeatures() const
{
  QgsFeatureList features;
  features.reserve( mSelectedFeatureIds.count() );
  QgsFeature f;

  if ( mSelectedFeatureIds.count() <= 8 )
  {
    // for small amount of selected features, fetch them directly
    // because request with FilterFids would go iterate over the whole layer
    const auto constMSelectedFeatureIds = mSelectedFeatureIds;
    for ( QgsFeatureId fid : constMSelectedFeatureIds )
    {
      getFeatures( QgsFeatureRequest( fid ) ).nextFeature( f );
      features << f;
    }
  }
  else
  {
    QgsFeatureIterator it = getSelectedFeatures();

    while ( it.nextFeature( f ) )
    {
      features.push_back( f );
    }
  }

  return features;
}

QgsFeatureIterator QgsVectorLayer::getSelectedFeatures( QgsFeatureRequest request ) const
{
  if ( mSelectedFeatureIds.isEmpty() )
    return QgsFeatureIterator();

  if ( geometryType() == QgsWkbTypes::NullGeometry )
    request.setFlags( QgsFeatureRequest::NoGeometry );

  if ( mSelectedFeatureIds.count() == 1 )
    request.setFilterFid( *mSelectedFeatureIds.constBegin() );
  else
    request.setFilterFids( mSelectedFeatureIds );

  return getFeatures( request );
}

bool QgsVectorLayer::addFeatures( QgsFeatureList &features, Flags )
{
  if ( !mEditBuffer || !mDataProvider )
    return false;

  if ( mGeometryOptions->isActive() )
  {
    for ( auto feature = features.begin(); feature != features.end(); ++feature )
    {
      QgsGeometry geom = feature->geometry();
      mGeometryOptions->apply( geom );
      feature->setGeometry( geom );
    }
  }

  bool res = mEditBuffer->addFeatures( features );
  updateExtents();

  if ( res && mJoinBuffer->containsJoins() )
    res = mJoinBuffer->addFeatures( features );

  return res;
}

void QgsVectorLayer::setCoordinateSystem()
{
  // if layer is not spatial, it has not CRS!
  setCrs( isSpatial() ? mDataProvider->crs() : QgsCoordinateReferenceSystem() );
}

QString QgsVectorLayer::displayField() const
{
  QgsExpression exp( mDisplayExpression );
  if ( exp.isField() )
  {
    return static_cast<const QgsExpressionNodeColumnRef *>( exp.rootNode() )->name();
  }

  return QString();
}

void QgsVectorLayer::setDisplayExpression( const QString &displayExpression )
{
  if ( mDisplayExpression == displayExpression )
    return;

  mDisplayExpression = displayExpression;
  emit displayExpressionChanged();
}

QString QgsVectorLayer::displayExpression() const
{
  if ( !mDisplayExpression.isEmpty() || mFields.isEmpty() )
  {
    return mDisplayExpression;
  }
  else
  {
    QString idxName;

    const auto constMFields = mFields;
    for ( const QgsField &field : constMFields )
    {
      QString fldName = field.name();

      // Check the fields and keep the first one that matches.
      // We assume that the user has organized the data with the
      // more "interesting" field names first. As such, name should
      // be selected before oldname, othername, etc.
      if ( fldName.indexOf( QLatin1String( "name" ), 0, Qt::CaseInsensitive ) > -1 )
      {
        idxName = fldName;
        break;
      }
      if ( fldName.indexOf( QLatin1String( "descrip" ), 0, Qt::CaseInsensitive ) > -1 )
      {
        idxName = fldName;
        break;
      }
      if ( fldName.indexOf( QLatin1String( "id" ), 0, Qt::CaseInsensitive ) > -1 )
      {
        idxName = fldName;
        break;
      }
    }

    if ( !idxName.isNull() )
    {
      return QgsExpression::quotedColumnRef( idxName );
    }
    else
    {
      return QgsExpression::quotedColumnRef( mFields.at( 0 ).name() );
    }
  }
}

bool QgsVectorLayer::isEditable() const
{
  return ( mEditBuffer && mDataProvider );
}

bool QgsVectorLayer::isSpatial() const
{
  QgsWkbTypes::GeometryType t = geometryType();
  return t != QgsWkbTypes::NullGeometry && t != QgsWkbTypes::UnknownGeometry;
}

bool QgsVectorLayer::isReadOnly() const
{
  return mReadOnly;
}

bool QgsVectorLayer::setReadOnly( bool readonly )
{
  // exit if the layer is in editing mode
  if ( readonly && mEditBuffer )
    return false;

  mReadOnly = readonly;
  emit readOnlyChanged();
  return true;
}

bool QgsVectorLayer::isModified() const
{
  emit beforeModifiedCheck();
  return mEditBuffer && mEditBuffer->isModified();
}

bool QgsVectorLayer::isAuxiliaryField( int index, int &srcIndex ) const
{
  bool auxiliaryField = false;
  srcIndex = -1;

  if ( !auxiliaryLayer() )
    return auxiliaryField;

  if ( index >= 0 && fields().fieldOrigin( index ) == QgsFields::OriginJoin )
  {
    const QgsVectorLayerJoinInfo *info = mJoinBuffer->joinForFieldIndex( index, fields(), srcIndex );

    if ( info && info->joinLayerId() == auxiliaryLayer()->id() )
      auxiliaryField = true;
  }

  return auxiliaryField;
}

void QgsVectorLayer::setRenderer( QgsFeatureRenderer *r )
{
  if ( !isSpatial() )
    return;

  if ( r != mRenderer )
  {
    delete mRenderer;
    mRenderer = r;
    mSymbolFeatureCounted = false;
    mSymbolFeatureCountMap.clear();

    emit rendererChanged();
    emit styleChanged();
  }
}

void QgsVectorLayer::beginEditCommand( const QString &text )
{
  if ( !mDataProvider )
  {
    return;
  }
  if ( mDataProvider->transaction() )
  {
    QString ignoredError;
    mDataProvider->transaction()->createSavepoint( ignoredError );
  }
  undoStack()->beginMacro( text );
  mEditCommandActive = true;
  emit editCommandStarted( text );
}

void QgsVectorLayer::endEditCommand()
{
  if ( !mDataProvider )
  {
    return;
  }
  undoStack()->endMacro();
  mEditCommandActive = false;
  if ( !mDeletedFids.isEmpty() )
  {
    emit featuresDeleted( mDeletedFids );
    mDeletedFids.clear();
  }
  emit editCommandEnded();
}

void QgsVectorLayer::destroyEditCommand()
{
  if ( !mDataProvider )
  {
    return;
  }
  undoStack()->endMacro();
  undoStack()->undo();

  // it's not directly possible to pop the last command off the stack (the destroyed one)
  // and delete, so we add a dummy obsolete command to force this to occur.
  // Pushing the new command deletes the destroyed one, and since the new
  // command is obsolete it's automatically deleted by the undo stack.
  std::unique_ptr< QUndoCommand > command = qgis::make_unique< QUndoCommand >();
  command->setObsolete( true );
  undoStack()->push( command.release() );

  mEditCommandActive = false;
  mDeletedFids.clear();
  emit editCommandDestroyed();
}

bool QgsVectorLayer::addJoin( const QgsVectorLayerJoinInfo &joinInfo )
{
  return mJoinBuffer->addJoin( joinInfo );
}


bool QgsVectorLayer::removeJoin( const QString &joinLayerId )
{
  return mJoinBuffer->removeJoin( joinLayerId );
}

const QList< QgsVectorLayerJoinInfo > QgsVectorLayer::vectorJoins() const
{
  return mJoinBuffer->vectorJoins();
}

int QgsVectorLayer::addExpressionField( const QString &exp, const QgsField &fld )
{
  emit beforeAddingExpressionField( fld.name() );
  mExpressionFieldBuffer->addExpression( exp, fld );
  updateFields();
  int idx = mFields.indexFromName( fld.name() );
  emit attributeAdded( idx );
  return idx;
}

void QgsVectorLayer::removeExpressionField( int index )
{
  emit beforeRemovingExpressionField( index );
  int oi = mFields.fieldOriginIndex( index );
  mExpressionFieldBuffer->removeExpression( oi );
  updateFields();
  emit attributeDeleted( index );
}

QString QgsVectorLayer::expressionField( int index ) const
{
  int oi = mFields.fieldOriginIndex( index );
  if ( oi < 0 || oi >= mExpressionFieldBuffer->expressions().size() )
    return QString();

  return mExpressionFieldBuffer->expressions().at( oi ).cachedExpression.expression();
}

void QgsVectorLayer::updateExpressionField( int index, const QString &exp )
{
  int oi = mFields.fieldOriginIndex( index );
  mExpressionFieldBuffer->updateExpression( oi, exp );
}

void QgsVectorLayer::updateFields()
{
  if ( !mDataProvider )
    return;

  QgsFields oldFields = mFields;

  mFields = mDataProvider->fields();

  // added / removed fields
  if ( mEditBuffer )
    mEditBuffer->updateFields( mFields );

  // joined fields
  if ( mJoinBuffer->containsJoins() )
    mJoinBuffer->updateFields( mFields );

  if ( mExpressionFieldBuffer )
    mExpressionFieldBuffer->updateFields( mFields );

  // set aliases and default values
  QMap< QString, QString >::const_iterator aliasIt = mAttributeAliasMap.constBegin();
  for ( ; aliasIt != mAttributeAliasMap.constEnd(); ++aliasIt )
  {
    int index = mFields.lookupField( aliasIt.key() );
    if ( index < 0 )
      continue;

    mFields[ index ].setAlias( aliasIt.value() );
  }

  // Update default values
  mDefaultValueOnUpdateFields.clear();
  QMap< QString, QgsDefaultValue >::const_iterator defaultIt = mDefaultExpressionMap.constBegin();
  for ( ; defaultIt != mDefaultExpressionMap.constEnd(); ++defaultIt )
  {
    int index = mFields.lookupField( defaultIt.key() );
    if ( index < 0 )
      continue;

    mFields[ index ].setDefaultValueDefinition( defaultIt.value() );
    if ( defaultIt.value().applyOnUpdate() )
      mDefaultValueOnUpdateFields.insert( index );
  }

  QMap< QString, QgsFieldConstraints::Constraints >::const_iterator constraintIt = mFieldConstraints.constBegin();
  for ( ; constraintIt != mFieldConstraints.constEnd(); ++constraintIt )
  {
    int index = mFields.lookupField( constraintIt.key() );
    if ( index < 0 )
      continue;

    QgsFieldConstraints constraints = mFields.at( index ).constraints();

    // always keep provider constraints intact
    if ( !( constraints.constraints() & QgsFieldConstraints::ConstraintNotNull ) && ( constraintIt.value() & QgsFieldConstraints::ConstraintNotNull ) )
      constraints.setConstraint( QgsFieldConstraints::ConstraintNotNull, QgsFieldConstraints::ConstraintOriginLayer );
    if ( !( constraints.constraints() & QgsFieldConstraints::ConstraintUnique ) && ( constraintIt.value() & QgsFieldConstraints::ConstraintUnique ) )
      constraints.setConstraint( QgsFieldConstraints::ConstraintUnique, QgsFieldConstraints::ConstraintOriginLayer );
    if ( !( constraints.constraints() & QgsFieldConstraints::ConstraintExpression ) && ( constraintIt.value() & QgsFieldConstraints::ConstraintExpression ) )
      constraints.setConstraint( QgsFieldConstraints::ConstraintExpression, QgsFieldConstraints::ConstraintOriginLayer );
    mFields[ index ].setConstraints( constraints );
  }

  QMap< QString, QPair< QString, QString > >::const_iterator constraintExpIt = mFieldConstraintExpressions.constBegin();
  for ( ; constraintExpIt != mFieldConstraintExpressions.constEnd(); ++constraintExpIt )
  {
    int index = mFields.lookupField( constraintExpIt.key() );
    if ( index < 0 )
      continue;

    QgsFieldConstraints constraints = mFields.at( index ).constraints();

    // always keep provider constraints intact
    if ( constraints.constraintOrigin( QgsFieldConstraints::ConstraintExpression ) == QgsFieldConstraints::ConstraintOriginProvider )
      continue;

    constraints.setConstraintExpression( constraintExpIt.value().first, constraintExpIt.value().second );
    mFields[ index ].setConstraints( constraints );
  }

  QMap< QPair< QString, QgsFieldConstraints::Constraint >, QgsFieldConstraints::ConstraintStrength >::const_iterator constraintStrengthIt = mFieldConstraintStrength.constBegin();
  for ( ; constraintStrengthIt != mFieldConstraintStrength.constEnd(); ++constraintStrengthIt )
  {
    int index = mFields.lookupField( constraintStrengthIt.key().first );
    if ( index < 0 )
      continue;

    QgsFieldConstraints constraints = mFields.at( index ).constraints();

    // always keep provider constraints intact
    if ( constraints.constraintOrigin( QgsFieldConstraints::ConstraintExpression ) == QgsFieldConstraints::ConstraintOriginProvider )
      continue;

    constraints.setConstraintStrength( constraintStrengthIt.key().second, constraintStrengthIt.value() );
    mFields[ index ].setConstraints( constraints );
  }

  auto fieldWidgetIterator = mFieldWidgetSetups.constBegin();
  for ( ; fieldWidgetIterator != mFieldWidgetSetups.constEnd(); ++ fieldWidgetIterator )
  {
    int index = mFields.indexOf( fieldWidgetIterator.key() );
    if ( index < 0 )
      continue;

    mFields[index].setEditorWidgetSetup( fieldWidgetIterator.value() );
  }

  if ( oldFields != mFields )
  {
    emit updatedFields();
    mEditFormConfig.setFields( mFields );
  }
}


QVariant QgsVectorLayer::defaultValue( int index, const QgsFeature &feature, QgsExpressionContext *context ) const
{
  if ( index < 0 || index >= mFields.count() )
    return QVariant();

  QString expression = mFields.at( index ).defaultValueDefinition().expression();
  if ( expression.isEmpty() )
    return mDataProvider->defaultValue( index );

  QgsExpressionContext *evalContext = context;
  std::unique_ptr< QgsExpressionContext > tempContext;
  if ( !evalContext )
  {
    // no context passed, so we create a default one
    tempContext.reset( new QgsExpressionContext( QgsExpressionContextUtils::globalProjectLayerScopes( this ) ) );
    evalContext = tempContext.get();
  }

  if ( feature.isValid() )
  {
    QgsExpressionContextScope *featScope = new QgsExpressionContextScope();
    featScope->setFeature( feature );
    featScope->setFields( feature.fields() );
    evalContext->appendScope( featScope );
  }

  QVariant val;
  QgsExpression exp( expression );
  exp.prepare( evalContext );
  if ( exp.hasEvalError() )
  {
    QgsLogger::warning( "Error evaluating default value: " + exp.evalErrorString() );
  }
  else
  {
    val = exp.evaluate( evalContext );
  }

  if ( feature.isValid() )
  {
    delete evalContext->popScope();
  }

  return val;
}

void QgsVectorLayer::setDefaultValueDefinition( int index, const QgsDefaultValue &definition )
{
  if ( index < 0 || index >= mFields.count() )
    return;

  if ( definition.isValid() )
  {
    mDefaultExpressionMap.insert( mFields.at( index ).name(), definition );
  }
  else
  {
    mDefaultExpressionMap.remove( mFields.at( index ).name() );
  }
  updateFields();
}

QgsDefaultValue QgsVectorLayer::defaultValueDefinition( int index ) const
{
  if ( index < 0 || index >= mFields.count() )
    return QgsDefaultValue();
  else
    return mFields.at( index ).defaultValueDefinition();
}

QSet<QVariant> QgsVectorLayer::uniqueValues( int index, int limit ) const
{
  QSet<QVariant> uniqueValues;
  if ( !mDataProvider )
  {
    return uniqueValues;
  }

  QgsFields::FieldOrigin origin = mFields.fieldOrigin( index );
  switch ( origin )
  {
    case QgsFields::OriginUnknown:
      return uniqueValues;

    case QgsFields::OriginProvider: //a provider field
    {
      uniqueValues = mDataProvider->uniqueValues( index, limit );

      if ( mEditBuffer )
      {
        QSet<QString> vals;
        const auto constUniqueValues = uniqueValues;
        for ( const QVariant &v : constUniqueValues )
        {
          vals << v.toString();
        }

        QgsFeatureMap added = mEditBuffer->addedFeatures();
        QMapIterator< QgsFeatureId, QgsFeature > addedIt( added );
        while ( addedIt.hasNext() && ( limit < 0 || uniqueValues.count() < limit ) )
        {
          addedIt.next();
          QVariant v = addedIt.value().attribute( index );
          if ( v.isValid() )
          {
            QString vs = v.toString();
            if ( !vals.contains( vs ) )
            {
              vals << vs;
              uniqueValues << v;
            }
          }
        }

        QMapIterator< QgsFeatureId, QgsAttributeMap > it( mEditBuffer->changedAttributeValues() );
        while ( it.hasNext() && ( limit < 0 || uniqueValues.count() < limit ) )
        {
          it.next();
          QVariant v = it.value().value( index );
          if ( v.isValid() )
          {
            QString vs = v.toString();
            if ( !vals.contains( vs ) )
            {
              vals << vs;
              uniqueValues << v;
            }
          }
        }
      }

      return uniqueValues;
    }

    case QgsFields::OriginEdit:
      // the layer is editable, but in certain cases it can still be avoided going through all features
      if ( mEditBuffer->mDeletedFeatureIds.isEmpty() &&
           mEditBuffer->mAddedFeatures.isEmpty() &&
           !mEditBuffer->mDeletedAttributeIds.contains( index ) &&
           mEditBuffer->mChangedAttributeValues.isEmpty() )
      {
        uniqueValues = mDataProvider->uniqueValues( index, limit );
        return uniqueValues;
      }
      FALLTHROUGH
    //we need to go through each feature
    case QgsFields::OriginJoin:
    case QgsFields::OriginExpression:
    {
      QgsAttributeList attList;
      attList << index;

      QgsFeatureIterator fit = getFeatures( QgsFeatureRequest()
                                            .setFlags( QgsFeatureRequest::NoGeometry )
                                            .setSubsetOfAttributes( attList ) );

      QgsFeature f;
      QVariant currentValue;
      QHash<QString, QVariant> val;
      while ( fit.nextFeature( f ) )
      {
        currentValue = f.attribute( index );
        val.insert( currentValue.toString(), currentValue );
        if ( limit >= 0 && val.size() >= limit )
        {
          break;
        }
      }

      return val.values().toSet();
    }
  }

  Q_ASSERT_X( false, "QgsVectorLayer::uniqueValues()", "Unknown source of the field!" );
  return uniqueValues;
}

QStringList QgsVectorLayer::uniqueStringsMatching( int index, const QString &substring, int limit, QgsFeedback *feedback ) const
{
  QStringList results;
  if ( !mDataProvider )
  {
    return results;
  }

  QgsFields::FieldOrigin origin = mFields.fieldOrigin( index );
  switch ( origin )
  {
    case QgsFields::OriginUnknown:
      return results;

    case QgsFields::OriginProvider: //a provider field
    {
      results = mDataProvider->uniqueStringsMatching( index, substring, limit, feedback );

      if ( mEditBuffer )
      {
        QgsFeatureMap added = mEditBuffer->addedFeatures();
        QMapIterator< QgsFeatureId, QgsFeature > addedIt( added );
        while ( addedIt.hasNext() && ( limit < 0 || results.count() < limit ) && ( !feedback || !feedback->isCanceled() ) )
        {
          addedIt.next();
          QVariant v = addedIt.value().attribute( index );
          if ( v.isValid() )
          {
            QString vs = v.toString();
            if ( vs.contains( substring, Qt::CaseInsensitive ) && !results.contains( vs ) )
            {
              results << vs;
            }
          }
        }

        QMapIterator< QgsFeatureId, QgsAttributeMap > it( mEditBuffer->changedAttributeValues() );
        while ( it.hasNext() && ( limit < 0 || results.count() < limit ) && ( !feedback || !feedback->isCanceled() ) )
        {
          it.next();
          QVariant v = it.value().value( index );
          if ( v.isValid() )
          {
            QString vs = v.toString();
            if ( vs.contains( substring, Qt::CaseInsensitive ) && !results.contains( vs ) )
            {
              results << vs;
            }
          }
        }
      }

      return results;
    }

    case QgsFields::OriginEdit:
      // the layer is editable, but in certain cases it can still be avoided going through all features
      if ( mEditBuffer->mDeletedFeatureIds.isEmpty() &&
           mEditBuffer->mAddedFeatures.isEmpty() &&
           !mEditBuffer->mDeletedAttributeIds.contains( index ) &&
           mEditBuffer->mChangedAttributeValues.isEmpty() )
      {
        return mDataProvider->uniqueStringsMatching( index, substring, limit, feedback );
      }
      FALLTHROUGH
    //we need to go through each feature
    case QgsFields::OriginJoin:
    case QgsFields::OriginExpression:
    {
      QgsAttributeList attList;
      attList << index;

      QgsFeatureRequest request;
      request.setSubsetOfAttributes( attList );
      request.setFlags( QgsFeatureRequest::NoGeometry );
      QString fieldName = mFields.at( index ).name();
      request.setFilterExpression( QStringLiteral( "\"%1\" ILIKE '%%2%'" ).arg( fieldName, substring ) );
      QgsFeatureIterator fit = getFeatures( request );

      QgsFeature f;
      QString currentValue;
      while ( fit.nextFeature( f ) )
      {
        currentValue = f.attribute( index ).toString();
        if ( !results.contains( currentValue ) )
          results << currentValue;

        if ( ( limit >= 0 && results.size() >= limit ) || ( feedback && feedback->isCanceled() ) )
        {
          break;
        }
      }

      return results;
    }
  }

  Q_ASSERT_X( false, "QgsVectorLayer::uniqueStringsMatching()", "Unknown source of the field!" );
  return results;
}

QVariant QgsVectorLayer::minimumValue( int index ) const
{
  if ( !mDataProvider )
  {
    return QVariant();
  }

  QgsFields::FieldOrigin origin = mFields.fieldOrigin( index );

  switch ( origin )
  {
    case QgsFields::OriginUnknown:
      return QVariant();

    case QgsFields::OriginProvider: //a provider field
    {
      QVariant min = mDataProvider->minimumValue( index );
      if ( mEditBuffer )
      {
        QgsFeatureMap added = mEditBuffer->addedFeatures();
        QMapIterator< QgsFeatureId, QgsFeature > addedIt( added );
        while ( addedIt.hasNext() )
        {
          addedIt.next();
          QVariant v = addedIt.value().attribute( index );
          if ( v.isValid() && qgsVariantLessThan( v, min ) )
          {
            min = v;
          }
        }

        QMapIterator< QgsFeatureId, QgsAttributeMap > it( mEditBuffer->changedAttributeValues() );
        while ( it.hasNext() )
        {
          it.next();
          QVariant v = it.value().value( index );
          if ( v.isValid() && qgsVariantLessThan( v, min ) )
          {
            min = v;
          }
        }
      }
      return min;
    }

    case QgsFields::OriginEdit:
    {
      // the layer is editable, but in certain cases it can still be avoided going through all features
      if ( mEditBuffer->mDeletedFeatureIds.isEmpty() &&
           mEditBuffer->mAddedFeatures.isEmpty() && !
           mEditBuffer->mDeletedAttributeIds.contains( index ) &&
           mEditBuffer->mChangedAttributeValues.isEmpty() )
      {
        return mDataProvider->minimumValue( index );
      }
    }
    FALLTHROUGH
    // no choice but to go through all features
    case QgsFields::OriginExpression:
    case QgsFields::OriginJoin:
    {
      // we need to go through each feature
      QgsAttributeList attList;
      attList << index;

      QgsFeatureIterator fit = getFeatures( QgsFeatureRequest()
                                            .setFlags( QgsFeatureRequest::NoGeometry )
                                            .setSubsetOfAttributes( attList ) );

      QgsFeature f;
      double minimumValue = std::numeric_limits<double>::max();
      double currentValue = 0;
      while ( fit.nextFeature( f ) )
      {
        currentValue = f.attribute( index ).toDouble();
        if ( currentValue < minimumValue )
        {
          minimumValue = currentValue;
        }
      }
      return QVariant( minimumValue );
    }
  }

  Q_ASSERT_X( false, "QgsVectorLayer::minimumValue()", "Unknown source of the field!" );
  return QVariant();
}

QVariant QgsVectorLayer::maximumValue( int index ) const
{
  if ( !mDataProvider )
  {
    return QVariant();
  }

  QgsFields::FieldOrigin origin = mFields.fieldOrigin( index );
  switch ( origin )
  {
    case QgsFields::OriginUnknown:
      return QVariant();

    case QgsFields::OriginProvider: //a provider field
    {
      QVariant min = mDataProvider->maximumValue( index );
      if ( mEditBuffer )
      {
        QgsFeatureMap added = mEditBuffer->addedFeatures();
        QMapIterator< QgsFeatureId, QgsFeature > addedIt( added );
        while ( addedIt.hasNext() )
        {
          addedIt.next();
          QVariant v = addedIt.value().attribute( index );
          if ( v.isValid() && qgsVariantGreaterThan( v, min ) )
          {
            min = v;
          }
        }

        QMapIterator< QgsFeatureId, QgsAttributeMap > it( mEditBuffer->changedAttributeValues() );
        while ( it.hasNext() )
        {
          it.next();
          QVariant v = it.value().value( index );
          if ( v.isValid() && qgsVariantGreaterThan( v, min ) )
          {
            min = v;
          }
        }
      }
      return min;
    }

    case QgsFields::OriginEdit:
      // the layer is editable, but in certain cases it can still be avoided going through all features
      if ( mEditBuffer->mDeletedFeatureIds.isEmpty() &&
           mEditBuffer->mAddedFeatures.isEmpty() &&
           !mEditBuffer->mDeletedAttributeIds.contains( index ) &&
           mEditBuffer->mChangedAttributeValues.isEmpty() )
      {
        return mDataProvider->maximumValue( index );
      }

      FALLTHROUGH
    //no choice but to go through each feature
    case QgsFields::OriginJoin:
    case QgsFields::OriginExpression:
    {
      QgsAttributeList attList;
      attList << index;

      QgsFeatureIterator fit = getFeatures( QgsFeatureRequest()
                                            .setFlags( QgsFeatureRequest::NoGeometry )
                                            .setSubsetOfAttributes( attList ) );

      QgsFeature f;
      double maximumValue = -std::numeric_limits<double>::max();
      double currentValue = 0;
      while ( fit.nextFeature( f ) )
      {
        currentValue = f.attribute( index ).toDouble();
        if ( currentValue > maximumValue )
        {
          maximumValue = currentValue;
        }
      }
      return QVariant( maximumValue );
    }
  }

  Q_ASSERT_X( false, "QgsVectorLayer::maximumValue()", "Unknown source of the field!" );
  return QVariant();
}

QVariant QgsVectorLayer::aggregate( QgsAggregateCalculator::Aggregate aggregate, const QString &fieldOrExpression,
                                    const QgsAggregateCalculator::AggregateParameters &parameters, QgsExpressionContext *context, bool *ok ) const
{
  if ( ok )
    *ok = false;

  if ( !mDataProvider )
  {
    return QVariant();
  }

  // test if we are calculating based on a field
  int attrIndex = mFields.lookupField( fieldOrExpression );
  if ( attrIndex >= 0 )
  {
    // aggregate is based on a field - if it's a provider field, we could possibly hand over the calculation
    // to the provider itself
    QgsFields::FieldOrigin origin = mFields.fieldOrigin( attrIndex );
    if ( origin == QgsFields::OriginProvider )
    {
      bool providerOk = false;
      QVariant val = mDataProvider->aggregate( aggregate, attrIndex, parameters, context, providerOk );
      if ( providerOk )
      {
        // provider handled calculation
        if ( ok )
          *ok = true;
        return val;
      }
    }
  }

  // fallback to using aggregate calculator to determine aggregate
  QgsAggregateCalculator c( this );
  c.setParameters( parameters );
  return c.calculate( aggregate, fieldOrExpression, context, ok );
}

void QgsVectorLayer::setFeatureBlendMode( QPainter::CompositionMode featureBlendMode )
{
  if ( mFeatureBlendMode == featureBlendMode )
    return;

  mFeatureBlendMode = featureBlendMode;
  emit featureBlendModeChanged( featureBlendMode );
  emit styleChanged();
}

QPainter::CompositionMode QgsVectorLayer::featureBlendMode() const
{
  return mFeatureBlendMode;
}

void QgsVectorLayer::setOpacity( double opacity )
{
  if ( qgsDoubleNear( mLayerOpacity, opacity ) )
    return;
  mLayerOpacity = opacity;
  emit opacityChanged( opacity );
  emit styleChanged();
}

double QgsVectorLayer::opacity() const
{
  return mLayerOpacity;
}



void QgsVectorLayer::readSldLabeling( const QDomNode &node )
{
  setLabeling( nullptr ); // start with no labeling
  setLabelsEnabled( false );

  QDomElement element = node.toElement();
  if ( element.isNull() )
    return;

  QDomElement userStyleElem = element.firstChildElement( QStringLiteral( "UserStyle" ) );
  if ( userStyleElem.isNull() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Info: UserStyle element not found." ), 4 );
    return;
  }

  QDomElement featureTypeStyleElem = userStyleElem.firstChildElement( QStringLiteral( "FeatureTypeStyle" ) );
  if ( featureTypeStyleElem.isNull() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Info: FeatureTypeStyle element not found." ), 4 );
    return;
  }

  // use last rule
  QDomElement ruleElem = featureTypeStyleElem.lastChildElement( QStringLiteral( "Rule" ) );
  if ( ruleElem.isNull() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Info: Rule element not found." ), 4 );
    return;
  }

  // use last text symbolizer
  QDomElement textSymbolizerElem = ruleElem.lastChildElement( QStringLiteral( "TextSymbolizer" ) );
  if ( textSymbolizerElem.isNull() )
  {
    QgsDebugMsgLevel( QStringLiteral( "Info: TextSymbolizer element not found." ), 4 );
    return;
  }

  QgsPalLayerSettings settings;

  // Label
  QDomElement labelElem = textSymbolizerElem.firstChildElement( QStringLiteral( "Label" ) );
  if ( !labelElem.isNull() )
  {
    QDomElement propertyNameElem = labelElem.firstChildElement( QStringLiteral( "PropertyName" ) );
    if ( !propertyNameElem.isNull() )
    {
      // set labeling defaults

      // label attribute
      QString labelAttribute = propertyNameElem.text();
      settings.fieldName = labelAttribute;
      settings.isExpression = false;

      int fieldIndex = mFields.lookupField( labelAttribute );
      if ( fieldIndex == -1 )
      {
        // label attribute is not in columns, check if it is an expression
        QgsExpression exp( labelAttribute );
        if ( !exp.hasEvalError() )
        {
          settings.isExpression = true;
        }
        else
        {
          QgsDebugMsgLevel( QStringLiteral( "SLD label attribute error: %1" ).arg( exp.evalErrorString() ), 3 );
        }
      }
    }
    else
    {
      QgsDebugMsgLevel( QStringLiteral( "Info: PropertyName element not found." ), 4 );
      return;
    }
  }
  else
  {
    QgsDebugMsgLevel( QStringLiteral( "Info: Label element not found." ), 4 );
    return;
  }

  QString fontFamily = QStringLiteral( "Sans-Serif" );
  int fontPointSize = 10;
  int fontWeight = -1;
  bool fontItalic = false;
  bool fontUnderline = false;

  // Font
  QDomElement fontElem = textSymbolizerElem.firstChildElement( QStringLiteral( "Font" ) );
  if ( !fontElem.isNull() )
  {
    QString cssName;
    QString elemText;
    QDomElement cssElem = fontElem.firstChildElement( QStringLiteral( "CssParameter" ) );
    while ( !cssElem.isNull() )
    {
      cssName = cssElem.attribute( QStringLiteral( "name" ), QStringLiteral( "not_found" ) );
      if ( cssName != QLatin1String( "not_found" ) )
      {
        elemText = cssElem.text();
        if ( cssName == QLatin1String( "font-family" ) )
        {
          fontFamily = elemText;
        }
        else if ( cssName == QLatin1String( "font-style" ) )
        {
          fontItalic = ( elemText == QLatin1String( "italic" ) ) || ( elemText == QLatin1String( "Italic" ) );
        }
        else if ( cssName == QLatin1String( "font-size" ) )
        {
          bool ok;
          int fontSize = elemText.toInt( &ok );
          if ( ok )
          {
            fontPointSize = fontSize;
          }
        }
        else if ( cssName == QLatin1String( "font-weight" ) )
        {
          if ( ( elemText == QLatin1String( "bold" ) ) || ( elemText == QLatin1String( "Bold" ) ) )
            fontWeight = QFont::Bold;
        }
        else if ( cssName == QLatin1String( "font-underline" ) )
        {
          fontUnderline = ( elemText == QLatin1String( "underline" ) ) || ( elemText == QLatin1String( "Underline" ) );
        }
      }

      cssElem = cssElem.nextSiblingElement( QStringLiteral( "CssParameter" ) );
    }
  }

  QgsTextFormat format;
  QFont font( fontFamily, fontPointSize, fontWeight, fontItalic );
  font.setUnderline( fontUnderline );
  format.setFont( font );
  format.setSize( fontPointSize );

  // Fill
  QColor textColor = QgsOgcUtils::colorFromOgcFill( textSymbolizerElem.firstChildElement( QStringLiteral( "Fill" ) ) );
  if ( textColor.isValid() )
  {
    format.setColor( textColor );
  }

  QgsTextBufferSettings bufferSettings;

  // Halo
  QDomElement haloElem = textSymbolizerElem.firstChildElement( QStringLiteral( "Halo" ) );
  if ( !haloElem.isNull() )
  {
    bufferSettings.setEnabled( true );
    bufferSettings.setSize( 1 );

    QDomElement radiusElem = haloElem.firstChildElement( QStringLiteral( "Radius" ) );
    if ( !radiusElem.isNull() )
    {
      bool ok;
      double bufferSize = radiusElem.text().toDouble( &ok );
      if ( ok )
      {
        bufferSettings.setSize( bufferSize );
      }
    }

    QColor bufferColor = QgsOgcUtils::colorFromOgcFill( haloElem.firstChildElement( QStringLiteral( "Fill" ) ) );
    if ( bufferColor.isValid() )
    {
      bufferSettings.setColor( bufferColor );
    }
  }

  // LabelPlacement
  QDomElement labelPlacementElem = textSymbolizerElem.firstChildElement( QStringLiteral( "LabelPlacement" ) );
  if ( !labelPlacementElem.isNull() )
  {
    // PointPlacement
    QDomElement pointPlacementElem = labelPlacementElem.firstChildElement( QStringLiteral( "PointPlacement" ) );
    if ( !pointPlacementElem.isNull() )
    {
      settings.placement = QgsPalLayerSettings::OverPoint;

      QDomElement displacementElem = pointPlacementElem.firstChildElement( QStringLiteral( "Displacement" ) );
      if ( !displacementElem.isNull() )
      {
        QDomElement displacementXElem = displacementElem.firstChildElement( QStringLiteral( "DisplacementX" ) );
        if ( !displacementXElem.isNull() )
        {
          bool ok;
          double xOffset = displacementXElem.text().toDouble( &ok );
          if ( ok )
          {
            settings.xOffset = xOffset;
          }
        }
        QDomElement displacementYElem = displacementElem.firstChildElement( QStringLiteral( "DisplacementY" ) );
        if ( !displacementYElem.isNull() )
        {
          bool ok;
          double yOffset = displacementYElem.text().toDouble( &ok );
          if ( ok )
          {
            settings.yOffset = yOffset;
          }
        }
      }

      QDomElement rotationElem = pointPlacementElem.firstChildElement( QStringLiteral( "Rotation" ) );
      if ( !rotationElem.isNull() )
      {
        bool ok;
        double rotation = rotationElem.text().toDouble( &ok );
        if ( ok )
        {
          settings.angleOffset = 360 - rotation;
        }
      }
    }
  }

  format.setBuffer( bufferSettings );
  settings.setFormat( format );
  setLabeling( new QgsVectorLayerSimpleLabeling( settings ) );
  setLabelsEnabled( true );
}

QgsEditFormConfig QgsVectorLayer::editFormConfig() const
{
  return mEditFormConfig;
}

void QgsVectorLayer::setEditFormConfig( const QgsEditFormConfig &editFormConfig )
{
  if ( mEditFormConfig == editFormConfig )
    return;

  mEditFormConfig = editFormConfig;
  mEditFormConfig.onRelationsLoaded();
  emit editFormConfigChanged();
}

QString QgsVectorLayer::mapTipTemplate() const
{
  return mMapTipTemplate;
}

void QgsVectorLayer::setMapTipTemplate( const QString &mapTip )
{
  if ( mMapTipTemplate == mapTip )
    return;

  mMapTipTemplate = mapTip;
  emit mapTipTemplateChanged();
}

QgsAttributeTableConfig QgsVectorLayer::attributeTableConfig() const
{
  QgsAttributeTableConfig config = mAttributeTableConfig;

  if ( config.isEmpty() )
    config.update( fields() );

  return config;
}

void QgsVectorLayer::setAttributeTableConfig( const QgsAttributeTableConfig &attributeTableConfig )
{
  if ( mAttributeTableConfig != attributeTableConfig )
  {
    mAttributeTableConfig = attributeTableConfig;
    emit configChanged();
  }
}

QgsExpressionContext QgsVectorLayer::createExpressionContext() const
{
  return QgsExpressionContext( QgsExpressionContextUtils::globalProjectLayerScopes( this ) );
}

QgsExpressionContextScope *QgsVectorLayer::createExpressionContextScope() const
{
  return QgsExpressionContextUtils::layerScope( this );
}

void QgsVectorLayer::setDiagramLayerSettings( const QgsDiagramLayerSettings &s )
{
  if ( !mDiagramLayerSettings )
    mDiagramLayerSettings = new QgsDiagramLayerSettings();
  *mDiagramLayerSettings = s;
}

QString QgsVectorLayer::htmlMetadata() const
{
  QgsLayerMetadataFormatter htmlFormatter( metadata() );
  QString myMetadata = QStringLiteral( "<html>\n<body>\n" );

  // Begin Provider section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "Information from provider" ) + QStringLiteral( "</h1>\n<hr>\n" );
  myMetadata += QLatin1String( "<table class=\"list-view\">\n" );

  // name
  myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Name" ) + QStringLiteral( "</td><td>" ) + name() + QStringLiteral( "</td></tr>\n" );

  // local path
  QVariantMap uriComponents = QgsProviderRegistry::instance()->decodeUri( mProviderKey, publicSource() );
  QString path;
  if ( uriComponents.contains( QStringLiteral( "path" ) ) )
  {
    path = uriComponents[QStringLiteral( "path" )].toString();
    if ( QFile::exists( path ) )
      myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Path" ) + QStringLiteral( "</td><td>%1" ).arg( QStringLiteral( "<a href=\"%1\">%2</a>" ).arg( QUrl::fromLocalFile( path ).toString(), QDir::toNativeSeparators( path ) ) ) + QStringLiteral( "</td></tr>\n" );
  }
  if ( uriComponents.contains( QStringLiteral( "url" ) ) )
  {
    const QString url = uriComponents[QStringLiteral( "url" )].toString();
    myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "URL" ) + QStringLiteral( "</td><td>%1" ).arg( QStringLiteral( "<a href=\"%1\">%2</a>" ).arg( QUrl( url ).toString(), url ) ) + QStringLiteral( "</td></tr>\n" );
  }

  // data source
  if ( publicSource() != path )
    myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Source" ) + QStringLiteral( "</td><td>%1" ).arg( publicSource() ) + QStringLiteral( "</td></tr>\n" );

  // storage type
  myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Storage" ) + QStringLiteral( "</td><td>" ) + storageType() + QStringLiteral( "</td></tr>\n" );

  // comment
  myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Comment" ) + QStringLiteral( "</td><td>" ) + dataComment() + QStringLiteral( "</td></tr>\n" );

  // encoding
  myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Encoding" ) + QStringLiteral( "</td><td>" ) + dataProvider()->encoding() + QStringLiteral( "</td></tr>\n" );

  if ( isSpatial() )
  {
    // geom type
    QgsWkbTypes::GeometryType type = geometryType();
    if ( type < 0 || type > QgsWkbTypes::NullGeometry )
    {
      QgsDebugMsgLevel( QStringLiteral( "Invalid vector type" ), 2 );
    }
    else
    {
      QString typeString( QStringLiteral( "%1 (%2)" ).arg( QgsWkbTypes::geometryDisplayString( geometryType() ),
                          QgsWkbTypes::displayString( wkbType() ) ) );
      myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Geometry" ) + QStringLiteral( "</td><td>" ) + typeString + QStringLiteral( "</td></tr>\n" );
    }

    // EPSG
    myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "CRS" ) + QStringLiteral( "</td><td>" );
    if ( crs().isValid() )
    {
      myMetadata += crs().authid() + QStringLiteral( " - " );
      myMetadata += crs().description() + QStringLiteral( " - " );
      if ( crs().isGeographic() )
        myMetadata += tr( "Geographic" );
      else
        myMetadata += tr( "Projected" );
    }
    myMetadata += QLatin1String( "</td></tr>\n" );

    // Extent
    myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Extent" ) + QStringLiteral( "</td><td>" ) + extent().toString() + QStringLiteral( "</td></tr>\n" );

    // unit
    myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Unit" ) + QStringLiteral( "</td><td>" ) + QgsUnitTypes::toString( crs().mapUnits() ) + QStringLiteral( "</td></tr>\n" );

  }

  // feature count
  QLocale locale = QLocale();
  locale.setNumberOptions( locale.numberOptions() &= ~QLocale::NumberOption::OmitGroupSeparator );
  myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" )
                + tr( "Feature count" ) + QStringLiteral( "</td><td>" )
                + ( featureCount() == -1 ? tr( "unknown" ) : locale.toString( static_cast<qlonglong>( featureCount() ) ) )
                + QStringLiteral( "</td></tr>\n" );

  // End Provider section
  myMetadata += QLatin1String( "</table>\n<br><br>" );

  // identification section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "Identification" ) + QStringLiteral( "</h1>\n<hr>\n" );
  myMetadata += htmlFormatter.identificationSectionHtml( );
  myMetadata += QLatin1String( "<br><br>\n" );

  // extent section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "Extent" ) + QStringLiteral( "</h1>\n<hr>\n" );
  myMetadata += htmlFormatter.extentSectionHtml( isSpatial() );
  myMetadata += QLatin1String( "<br><br>\n" );

  // Start the Access section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "Access" ) + QStringLiteral( "</h1>\n<hr>\n" );
  myMetadata += htmlFormatter.accessSectionHtml( );
  myMetadata += QLatin1String( "<br><br>\n" );

  // Fields section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "Fields" ) + QStringLiteral( "</h1>\n<hr>\n<table class=\"list-view\">\n" );

  // primary key
  QgsAttributeList pkAttrList = primaryKeyAttributes();
  if ( !pkAttrList.isEmpty() )
  {
    myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Primary key attributes" ) + QStringLiteral( "</td><td>" );
    const auto constPkAttrList = pkAttrList;
    for ( int idx : constPkAttrList )
    {
      myMetadata += fields().at( idx ).name() + ' ';
    }
    myMetadata += QLatin1String( "</td></tr>\n" );
  }

  const QgsFields myFields = fields();

  // count fields
  myMetadata += QStringLiteral( "<tr><td class=\"highlight\">" ) + tr( "Count" ) + QStringLiteral( "</td><td>" ) + QString::number( myFields.size() ) + QStringLiteral( "</td></tr>\n" );

  myMetadata += QLatin1String( "</table>\n<br><table width=\"100%\" class=\"tabular-view\">\n" );
  myMetadata += QLatin1String( "<tr><th>" ) + tr( "Field" ) + QLatin1String( "</th><th>" ) + tr( "Type" ) + QLatin1String( "</th><th>" ) + tr( "Length" ) + QLatin1String( "</th><th>" ) + tr( "Precision" ) + QLatin1String( "</th><th>" ) + tr( "Comment" ) + QLatin1String( "</th></tr>\n" );

  for ( int i = 0; i < myFields.size(); ++i )
  {
    QgsField myField = myFields.at( i );
    QString rowClass;
    if ( i % 2 )
      rowClass = QStringLiteral( "class=\"odd-row\"" );
    myMetadata += QLatin1String( "<tr " ) + rowClass + QLatin1String( "><td>" ) + myField.name() + QLatin1String( "</td><td>" ) + myField.typeName() + QLatin1String( "</td><td>" ) + QString::number( myField.length() ) + QLatin1String( "</td><td>" ) + QString::number( myField.precision() ) + QLatin1String( "</td><td>" ) + myField.comment() + QLatin1String( "</td></tr>\n" );
  }

  //close field list
  myMetadata += QLatin1String( "</table>\n<br><br>" );

  // Start the contacts section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "Contacts" ) + QStringLiteral( "</h1>\n<hr>\n" );
  myMetadata += htmlFormatter.contactsSectionHtml( );
  myMetadata += QLatin1String( "<br><br>\n" );

  // Start the links section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "Links" ) + QStringLiteral( "</h1>\n<hr>\n" );
  myMetadata += htmlFormatter.linksSectionHtml( );
  myMetadata += QLatin1String( "<br><br>\n" );

  // Start the history section
  myMetadata += QStringLiteral( "<h1>" ) + tr( "History" ) + QStringLiteral( "</h1>\n<hr>\n" );
  myMetadata += htmlFormatter.historySectionHtml( );
  myMetadata += QLatin1String( "<br><br>\n" );

  myMetadata += QStringLiteral( "\n</body>\n</html>\n" );
  return myMetadata;
}

void QgsVectorLayer::invalidateSymbolCountedFlag()
{
  mSymbolFeatureCounted = false;
}

void QgsVectorLayer::onFeatureCounterCompleted()
{
  onSymbolsCounted();
  mFeatureCounter = nullptr;
}

void QgsVectorLayer::onFeatureCounterTerminated()
{
  mFeatureCounter = nullptr;
}

void QgsVectorLayer::onJoinedFieldsChanged()
{
  // some of the fields of joined layers have changed -> we need to update this layer's fields too
  updateFields();
}

void QgsVectorLayer::onFeatureDeleted( QgsFeatureId fid )
{
  if ( mEditCommandActive )
    mDeletedFids << fid;
  else
    emit featuresDeleted( QgsFeatureIds() << fid );

  emit featureDeleted( fid );
}

void QgsVectorLayer::onRelationsLoaded()
{
  mEditFormConfig.onRelationsLoaded();
}

void QgsVectorLayer::onSymbolsCounted()
{
  if ( mFeatureCounter )
  {
    mSymbolFeatureCountMap = mFeatureCounter->symbolFeatureCountMap();
    mSymbolFeatureCounted = true;
    emit symbolFeatureCountMapChanged();
  }
}

QList<QgsRelation> QgsVectorLayer::referencingRelations( int idx ) const
{
  return QgsProject::instance()->relationManager()->referencingRelations( this, idx );
}

int QgsVectorLayer::listStylesInDatabase( QStringList &ids, QStringList &names, QStringList &descriptions, QString &msgError )
{
  std::unique_ptr<QLibrary> myLib( QgsProviderRegistry::instance()->createProviderLibrary( mProviderKey ) );
  if ( !myLib )
  {
    msgError = QObject::tr( "Unable to load %1 provider" ).arg( mProviderKey );
    return -1;
  }
  listStyles_t *listStylesExternalMethod = reinterpret_cast< listStyles_t * >( cast_to_fptr( myLib->resolve( "listStyles" ) ) );

  if ( !listStylesExternalMethod )
  {
    msgError = QObject::tr( "Provider %1 has no %2 method" ).arg( mProviderKey, QStringLiteral( "listStyles" ) );
    return -1;
  }

  return listStylesExternalMethod( mDataSource, ids, names, descriptions, msgError );
}

QString QgsVectorLayer::getStyleFromDatabase( const QString &styleId, QString &msgError )
{
  std::unique_ptr<QLibrary> myLib( QgsProviderRegistry::instance()->createProviderLibrary( mProviderKey ) );
  if ( !myLib )
  {
    msgError = QObject::tr( "Unable to load %1 provider" ).arg( mProviderKey );
    return QString();
  }
  getStyleById_t *getStyleByIdMethod = reinterpret_cast< getStyleById_t * >( cast_to_fptr( myLib->resolve( "getStyleById" ) ) );

  if ( !getStyleByIdMethod )
  {
    msgError = QObject::tr( "Provider %1 has no %2 method" ).arg( mProviderKey, QStringLiteral( "getStyleById" ) );
    return QString();
  }

  return getStyleByIdMethod( mDataSource, styleId, msgError );
}

bool QgsVectorLayer::deleteStyleFromDatabase( const QString &styleId, QString &msgError )
{
  std::unique_ptr<QLibrary> myLib( QgsProviderRegistry::instance()->createProviderLibrary( mProviderKey ) );
  if ( !myLib )
  {
    msgError = QObject::tr( "Unable to load %1 provider" ).arg( mProviderKey );
    return false;
  }
  deleteStyleById_t *deleteStyleByIdMethod = reinterpret_cast< deleteStyleById_t * >( cast_to_fptr( myLib->resolve( "deleteStyleById" ) ) );
  if ( !deleteStyleByIdMethod )
  {
    msgError = QObject::tr( "Provider %1 has no %2 method" ).arg( mProviderKey, QStringLiteral( "deleteStyleById" ) );
    return false;
  }
  return deleteStyleByIdMethod( mDataSource, styleId, msgError );
}


void QgsVectorLayer::saveStyleToDatabase( const QString &name, const QString &description,
    bool useAsDefault, const QString &uiFileContent, QString &msgError )
{

  QString sldStyle, qmlStyle;
  std::unique_ptr<QLibrary> myLib( QgsProviderRegistry::instance()->createProviderLibrary( mProviderKey ) );
  if ( !myLib )
  {
    msgError = QObject::tr( "Unable to load %1 provider" ).arg( mProviderKey );
    return;
  }
  saveStyle_t *saveStyleExternalMethod = reinterpret_cast< saveStyle_t * >( cast_to_fptr( myLib->resolve( "saveStyle" ) ) );

  if ( !saveStyleExternalMethod )
  {
    msgError = QObject::tr( "Provider %1 has no %2 method" ).arg( mProviderKey, QStringLiteral( "saveStyle" ) );
    return;
  }

  QDomDocument qmlDocument, sldDocument;
  QgsReadWriteContext context;
  exportNamedStyle( qmlDocument, msgError, context );
  if ( !msgError.isNull() )
  {
    return;
  }
  qmlStyle = qmlDocument.toString();

  this->exportSldStyle( sldDocument, msgError );
  if ( !msgError.isNull() )
  {
    return;
  }
  sldStyle = sldDocument.toString();

  saveStyleExternalMethod( mDataSource, qmlStyle, sldStyle, name,
                           description, uiFileContent, useAsDefault, msgError );
}



QString QgsVectorLayer::loadNamedStyle( const QString &theURI, bool &resultFlag, QgsMapLayer::StyleCategories categories )
{
  return loadNamedStyle( theURI, resultFlag, false, categories );
}

bool QgsVectorLayer::loadAuxiliaryLayer( const QgsAuxiliaryStorage &storage, const QString &key )
{
  bool rc = false;

  QString joinKey = mAuxiliaryLayerKey;
  if ( !key.isEmpty() )
    joinKey = key;

  if ( storage.isValid() && !joinKey.isEmpty() )
  {
    QgsAuxiliaryLayer *alayer = nullptr;

    int idx = fields().lookupField( joinKey );

    if ( idx >= 0 )
    {
      alayer = storage.createAuxiliaryLayer( fields().field( idx ), this );

      if ( alayer )
      {
        setAuxiliaryLayer( alayer );
        rc = true;
      }
    }
  }

  return rc;
}

void QgsVectorLayer::setAuxiliaryLayer( QgsAuxiliaryLayer *alayer )
{
  mAuxiliaryLayerKey.clear();

  if ( mAuxiliaryLayer )
    removeJoin( mAuxiliaryLayer->id() );

  if ( alayer )
  {
    addJoin( alayer->joinInfo() );

    if ( !alayer->isEditable() )
      alayer->startEditing();

    mAuxiliaryLayerKey = alayer->joinInfo().targetFieldName();
  }

  mAuxiliaryLayer.reset( alayer );
  if ( mAuxiliaryLayer )
    mAuxiliaryLayer->setParent( this );
  updateFields();
}

const QgsAuxiliaryLayer *QgsVectorLayer::auxiliaryLayer() const
{
  return mAuxiliaryLayer.get();
}

QgsAuxiliaryLayer *QgsVectorLayer::auxiliaryLayer()
{
  return mAuxiliaryLayer.get();
}

QString QgsVectorLayer::loadNamedStyle( const QString &theURI, bool &resultFlag, bool loadFromLocalDB, QgsMapLayer::StyleCategories categories )
{
  QgsDataSourceUri dsUri( theURI );
  if ( !loadFromLocalDB && mDataProvider && mDataProvider->isSaveAndLoadStyleToDatabaseSupported() )
  {
    std::unique_ptr<QLibrary> myLib( QgsProviderRegistry::instance()->createProviderLibrary( mProviderKey ) );
    if ( myLib )
    {
      loadStyle_t *loadStyleExternalMethod = reinterpret_cast< loadStyle_t * >( cast_to_fptr( myLib->resolve( "loadStyle" ) ) );
      if ( loadStyleExternalMethod )
      {
        QString qml, errorMsg;
        qml = loadStyleExternalMethod( mDataSource, errorMsg );
        if ( !qml.isEmpty() )
        {
          QDomDocument myDocument( QStringLiteral( "qgis" ) );
          myDocument.setContent( qml );
          resultFlag = importNamedStyle( myDocument, errorMsg );
          return QObject::tr( "Loaded from Provider" );
        }
      }
    }
  }

  return QgsMapLayer::loadNamedStyle( theURI, resultFlag, categories );
}

QSet<QgsMapLayerDependency> QgsVectorLayer::dependencies() const
{
  if ( mDataProvider )
    return mDataProvider->dependencies() + mDependencies;
  return mDependencies;
}

bool QgsVectorLayer::setDependencies( const QSet<QgsMapLayerDependency> &oDeps )
{
  QSet<QgsMapLayerDependency> deps;
  const auto constODeps = oDeps;
  for ( const QgsMapLayerDependency &dep : constODeps )
  {
    if ( dep.origin() == QgsMapLayerDependency::FromUser )
      deps << dep;
  }
  if ( hasDependencyCycle( deps ) )
    return false;

  QSet<QgsMapLayerDependency> toAdd = deps - dependencies();

  // disconnect layers that are not present in the list of dependencies anymore
  for ( const QgsMapLayerDependency &dep : qgis::as_const( mDependencies ) )
  {
    QgsVectorLayer *lyr = static_cast<QgsVectorLayer *>( QgsProject::instance()->mapLayer( dep.layerId() ) );
    if ( !lyr )
      continue;
    disconnect( lyr, &QgsVectorLayer::featureAdded, this, &QgsVectorLayer::dataChanged );
    disconnect( lyr, &QgsVectorLayer::featureDeleted, this, &QgsVectorLayer::dataChanged );
    disconnect( lyr, &QgsVectorLayer::geometryChanged, this, &QgsVectorLayer::dataChanged );
    disconnect( lyr, &QgsVectorLayer::dataChanged, this, &QgsVectorLayer::dataChanged );
    disconnect( lyr, &QgsVectorLayer::repaintRequested, this, &QgsVectorLayer::triggerRepaint );
  }

  // assign new dependencies
  if ( mDataProvider )
    mDependencies = mDataProvider->dependencies() + deps;
  else
    mDependencies = deps;
  emit dependenciesChanged();

  // connect to new layers
  for ( const QgsMapLayerDependency &dep : qgis::as_const( mDependencies ) )
  {
    QgsVectorLayer *lyr = static_cast<QgsVectorLayer *>( QgsProject::instance()->mapLayer( dep.layerId() ) );
    if ( !lyr )
      continue;
    connect( lyr, &QgsVectorLayer::featureAdded, this, &QgsVectorLayer::dataChanged );
    connect( lyr, &QgsVectorLayer::featureDeleted, this, &QgsVectorLayer::dataChanged );
    connect( lyr, &QgsVectorLayer::geometryChanged, this, &QgsVectorLayer::dataChanged );
    connect( lyr, &QgsVectorLayer::dataChanged, this, &QgsVectorLayer::dataChanged );
    connect( lyr, &QgsVectorLayer::repaintRequested, this, &QgsVectorLayer::triggerRepaint );
  }

  // if new layers are present, emit a data change
  if ( ! toAdd.isEmpty() )
    emit dataChanged();

  return true;
}

QgsFieldConstraints::Constraints QgsVectorLayer::fieldConstraints( int fieldIndex ) const
{
  if ( fieldIndex < 0 || fieldIndex >= mFields.count() )
    return nullptr;

  QgsFieldConstraints::Constraints constraints = mFields.at( fieldIndex ).constraints().constraints();

  // make sure provider constraints are always present!
  if ( mFields.fieldOrigin( fieldIndex ) == QgsFields::OriginProvider )
  {
    constraints |= mDataProvider->fieldConstraints( mFields.fieldOriginIndex( fieldIndex ) );
  }

  return constraints;
}

QMap< QgsFieldConstraints::Constraint, QgsFieldConstraints::ConstraintStrength> QgsVectorLayer::fieldConstraintsAndStrength( int fieldIndex ) const
{
  QMap< QgsFieldConstraints::Constraint, QgsFieldConstraints::ConstraintStrength > m;

  if ( fieldIndex < 0 || fieldIndex >= mFields.count() )
    return m;

  QString name = mFields.at( fieldIndex ).name();

  QMap< QPair< QString, QgsFieldConstraints::Constraint >, QgsFieldConstraints::ConstraintStrength >::const_iterator conIt = mFieldConstraintStrength.constBegin();
  for ( ; conIt != mFieldConstraintStrength.constEnd(); ++conIt )
  {
    if ( conIt.key().first == name )
    {
      m[ conIt.key().second ] = mFieldConstraintStrength.value( conIt.key() );
    }
  }

  return m;
}

void QgsVectorLayer::setFieldConstraint( int index, QgsFieldConstraints::Constraint constraint, QgsFieldConstraints::ConstraintStrength strength )
{
  if ( index < 0 || index >= mFields.count() )
    return;

  QString name = mFields.at( index ).name();

  // add constraint to existing constraints
  QgsFieldConstraints::Constraints constraints = mFieldConstraints.value( name, nullptr );
  constraints |= constraint;
  mFieldConstraints.insert( name, constraints );

  mFieldConstraintStrength.insert( qMakePair( name, constraint ), strength );

  updateFields();
}

void QgsVectorLayer::removeFieldConstraint( int index, QgsFieldConstraints::Constraint constraint )
{
  if ( index < 0 || index >= mFields.count() )
    return;

  QString name = mFields.at( index ).name();

  // remove constraint from existing constraints
  QgsFieldConstraints::Constraints constraints = mFieldConstraints.value( name, nullptr );
  constraints &= ~constraint;
  mFieldConstraints.insert( name, constraints );

  mFieldConstraintStrength.remove( qMakePair( name, constraint ) );

  updateFields();
}

QString QgsVectorLayer::constraintExpression( int index ) const
{
  if ( index < 0 || index >= mFields.count() )
    return QString();

  return mFields.at( index ).constraints().constraintExpression();
}

QString QgsVectorLayer::constraintDescription( int index ) const
{
  if ( index < 0 || index >= mFields.count() )
    return QString();

  return mFields.at( index ).constraints().constraintDescription();
}

void QgsVectorLayer::setConstraintExpression( int index, const QString &expression, const QString &description )
{
  if ( index < 0 || index >= mFields.count() )
    return;

  if ( expression.isEmpty() )
  {
    mFieldConstraintExpressions.remove( mFields.at( index ).name() );
  }
  else
  {
    mFieldConstraintExpressions.insert( mFields.at( index ).name(), qMakePair( expression, description ) );
  }
  updateFields();
}

void QgsVectorLayer::setEditorWidgetSetup( int index, const QgsEditorWidgetSetup &setup )
{
  if ( index < 0 || index >= mFields.count() )
    return;

  if ( setup.isNull() )
    mFieldWidgetSetups.remove( mFields.at( index ).name() );
  else
    mFieldWidgetSetups.insert( mFields.at( index ).name(), setup );
  updateFields();
}

QgsEditorWidgetSetup QgsVectorLayer::editorWidgetSetup( int index ) const
{

  if ( index < 0 || index >= mFields.count() )
    return QgsEditorWidgetSetup();

  return mFields.at( index ).editorWidgetSetup();
}

QgsAbstractVectorLayerLabeling *QgsVectorLayer::readLabelingFromCustomProperties()
{
  QgsAbstractVectorLayerLabeling *labeling = nullptr;
  if ( customProperty( QStringLiteral( "labeling" ) ).toString() == QLatin1String( "pal" ) )
  {
    if ( customProperty( QStringLiteral( "labeling/enabled" ), QVariant( false ) ).toBool() )
    {
      // try to load from custom properties
      QgsPalLayerSettings settings;
      settings.readFromLayerCustomProperties( this );
      labeling = new QgsVectorLayerSimpleLabeling( settings );
    }

    // also clear old-style labeling config
    removeCustomProperty( QStringLiteral( "labeling" ) );
    const auto constCustomPropertyKeys = customPropertyKeys();
    for ( const QString &key : constCustomPropertyKeys )
    {
      if ( key.startsWith( QLatin1String( "labeling/" ) ) )
        removeCustomProperty( key );
    }
  }

  return labeling;
}

bool QgsVectorLayer::allowCommit() const
{
  return mAllowCommit;
}

void QgsVectorLayer::setAllowCommit( bool allowCommit )
{
  if ( mAllowCommit == allowCommit )
    return;

  mAllowCommit = allowCommit;
  emit allowCommitChanged();
}

QgsGeometryOptions *QgsVectorLayer::geometryOptions() const
{
  return mGeometryOptions.get();
}

void QgsVectorLayer::setReadExtentFromXml( bool readExtentFromXml )
{
  mReadExtentFromXml = readExtentFromXml;
}

bool QgsVectorLayer::readExtentFromXml() const
{
  return mReadExtentFromXml;
}

void QgsVectorLayer::onDirtyTransaction( const QString &sql, const QString &name )
{
  QgsTransaction *tr = dataProvider()->transaction();
  if ( tr && mEditBuffer )
  {
    qobject_cast<QgsVectorLayerEditPassthrough *>( mEditBuffer )->update( tr, sql, name );
  }
}
