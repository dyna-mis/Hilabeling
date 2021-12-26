/***************************************************************************
                         testqgsmapthemecollection.cpp
                         ----------------------
    begin                : April 2018
    copyright            : (C) 2018 by Martin Dobias
    email                : wonder dot sk at gmail dot com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgstest.h"
#include "qgslayertree.h"
#include "qgslayertreemodel.h"
#include "qgsmapthemecollection.h"
#include "qgsrulebasedrenderer.h"
#include "qgsvectorlayer.h"

class TestQgsMapThemeCollection : public QObject
{
    Q_OBJECT

  public:
    TestQgsMapThemeCollection() = default;

  private slots:
    void initTestCase();// will be called before the first testfunction is executed.
    void cleanupTestCase();// will be called after the last testfunction was executed.

    void expandedState();

  private:
    QgsVectorLayer *mPointsLayer = nullptr;
    QgsVectorLayer *mPolysLayer = nullptr;
    QgsVectorLayer *mLinesLayer = nullptr;
    QgsLayerTree *mLayerTree = nullptr;
    QgsLayerTreeModel *mLayerTreeModel = nullptr;
    QgsProject *mProject = nullptr;

    QgsLayerTreeGroup *mNodeGroup1;
    QgsLayerTreeGroup *mNodeGroup2;
    QgsLayerTreeGroup *mNodeGroup3;
    QgsLayerTreeLayer *mNodeLayerPoints;
    QgsLayerTreeLayer *mNodeLayerLines;
    QgsLayerTreeLayer *mNodeLayerPolys;
};

void TestQgsMapThemeCollection::initTestCase()
{
  QgsApplication::init();
  QgsApplication::initQgis();

  QFileInfo pointFileInfo( QStringLiteral( TEST_DATA_DIR ) + "/points.shp" );
  mPointsLayer = new QgsVectorLayer( pointFileInfo.filePath(),
                                     pointFileInfo.completeBaseName(), QStringLiteral( "ogr" ) );

  QFileInfo polyFileInfo( QStringLiteral( TEST_DATA_DIR ) + "/polys.shp" );
  mPolysLayer = new QgsVectorLayer( polyFileInfo.filePath(),
                                    polyFileInfo.completeBaseName(), QStringLiteral( "ogr" ) );

  QFileInfo lineFileInfo( QStringLiteral( TEST_DATA_DIR ) + "/lines.shp" );
  mLinesLayer = new QgsVectorLayer( lineFileInfo.filePath(),
                                    lineFileInfo.completeBaseName(), QStringLiteral( "ogr" ) );

  // rule-based renderer for points layer with two levels of rules
  QgsRuleBasedRenderer::Rule *rule0 = new QgsRuleBasedRenderer::Rule( nullptr );
  QgsRuleBasedRenderer::Rule *ruleA = new QgsRuleBasedRenderer::Rule( nullptr, 0, 0, QString(), "A" );
  QgsRuleBasedRenderer::Rule *ruleB = new QgsRuleBasedRenderer::Rule( nullptr, 0, 0, QString(), "B" );
  QgsRuleBasedRenderer::Rule *ruleC = new QgsRuleBasedRenderer::Rule( nullptr, 0, 0, QString(), "C" );
  QgsRuleBasedRenderer::Rule *ruleD = new QgsRuleBasedRenderer::Rule( nullptr, 0, 0, QString(), "D" );
  rule0->appendChild( ruleA );
  ruleA->appendChild( ruleB );
  rule0->appendChild( ruleC );
  ruleC->appendChild( ruleD );
  mPointsLayer->setRenderer( new QgsRuleBasedRenderer( rule0 ) );

  mLayerTree = new QgsLayerTree;
  mNodeGroup1 = mLayerTree->addGroup( "group1" );
  mNodeGroup2 = mLayerTree->addGroup( "group2" );
  mNodeGroup3 = mNodeGroup2->addGroup( "group3" );
  mNodeLayerPoints = mNodeGroup1->addLayer( mPointsLayer );
  mNodeLayerLines = mNodeGroup2->addLayer( mLinesLayer );
  mNodeLayerPolys = mNodeGroup3->addLayer( mPolysLayer );

  mLayerTreeModel = new QgsLayerTreeModel( mLayerTree );

  // layers need to be in project to be accessible from map themes
  mProject = new QgsProject;
  mProject->addMapLayers( QList<QgsMapLayer *>() << mPointsLayer << mPolysLayer << mLinesLayer );
}

void TestQgsMapThemeCollection::cleanupTestCase()
{
  delete mLayerTreeModel;
  delete mLayerTree;
  delete mProject;

  QgsApplication::exitQgis();
}


static QgsMapThemeCollection::MapThemeLayerRecord _recordForLayer( QgsMapLayer *l, const QgsMapThemeCollection::MapThemeRecord &rec )
{
  const QList<QgsMapThemeCollection::MapThemeLayerRecord> layerRecords = rec.layerRecords();
  for ( const QgsMapThemeCollection::MapThemeLayerRecord &layerRec : layerRecords )
    if ( layerRec.layer() == l )
      return layerRec;
  Q_ASSERT( false );
  return QgsMapThemeCollection::MapThemeLayerRecord();
}


void TestQgsMapThemeCollection::expandedState()
{
  QgsMapThemeCollection themes( mProject );

  QStringList pointLayerRootLegendNodes;
  const QList<QgsLayerTreeModelLegendNode *> legendNodes = mLayerTreeModel->layerLegendNodes( mNodeLayerPoints, true );
  for ( QgsLayerTreeModelLegendNode *legendNode : legendNodes )
  {
    QString key = legendNode->data( QgsLayerTreeModelLegendNode::RuleKeyRole ).toString();
    pointLayerRootLegendNodes << key;
  }
  QCOMPARE( pointLayerRootLegendNodes.count(), 4 );

  // make two themes: 1. all expanded, 2. all collapsed

  mNodeGroup1->setExpanded( false );
  mNodeGroup2->setExpanded( false );
  mNodeGroup3->setExpanded( false );
  mNodeLayerPoints->setExpanded( false );
  mNodeLayerLines->setExpanded( false );
  mNodeLayerPolys->setExpanded( false );
  mNodeLayerPoints->setCustomProperty( "expandedLegendNodes", QStringList() );

  themes.insert( "all-collapsed", QgsMapThemeCollection::createThemeFromCurrentState( mLayerTree, mLayerTreeModel ) );

  mNodeGroup1->setExpanded( true );
  mNodeGroup2->setExpanded( true );
  mNodeGroup3->setExpanded( true );
  mNodeLayerPoints->setExpanded( true );
  mNodeLayerLines->setExpanded( true );
  mNodeLayerPolys->setExpanded( true );
  mNodeLayerPoints->setCustomProperty( "expandedLegendNodes", pointLayerRootLegendNodes );

  themes.insert( "all-expanded", QgsMapThemeCollection::createThemeFromCurrentState( mLayerTree, mLayerTreeModel ) );

  // check theme data

  QgsMapThemeCollection::MapThemeRecord recCollapsed = themes.mapThemeState( "all-collapsed" );
  QVERIFY( recCollapsed.hasExpandedStateInfo() );
  QCOMPARE( recCollapsed.expandedGroupNodes().count(), 0 );
  QgsMapThemeCollection::MapThemeLayerRecord recCollapsedPoints = _recordForLayer( mPointsLayer, recCollapsed );
  QCOMPARE( recCollapsedPoints.expandedLayerNode, false );
  QVERIFY( recCollapsedPoints.expandedLegendItems.isEmpty() );

  QgsMapThemeCollection::MapThemeRecord recExpanded = themes.mapThemeState( "all-expanded" );
  QVERIFY( recExpanded.hasExpandedStateInfo() );
  QCOMPARE( recExpanded.expandedGroupNodes().count(), 3 );
  QgsMapThemeCollection::MapThemeLayerRecord recExpandedPoints = _recordForLayer( mPointsLayer, recExpanded );
  QCOMPARE( recExpandedPoints.expandedLayerNode, true );
  QCOMPARE( recExpandedPoints.expandedLegendItems.count(), 4 );

  // switch themes

  themes.applyTheme( "all-collapsed", mLayerTree, mLayerTreeModel );
  QCOMPARE( mNodeGroup1->isExpanded(), false );
  QCOMPARE( mNodeGroup3->isExpanded(), false );
  QCOMPARE( mNodeLayerPolys->isExpanded(), false );
  QVERIFY( mNodeLayerPoints->customProperty( "expandedLegendNodes" ).toStringList().isEmpty() );

  themes.applyTheme( "all-expanded", mLayerTree, mLayerTreeModel );
  QCOMPARE( mNodeGroup1->isExpanded(), true );
  QCOMPARE( mNodeGroup3->isExpanded(), true );
  QCOMPARE( mNodeLayerPolys->isExpanded(), true );
  QVERIFY( mNodeLayerPoints->customProperty( "expandedLegendNodes" ).toStringList().count() == 4 );

  // test read+write

  QDomDocument doc;
  QDomElement elemRoot = doc.createElement( "qgis" );
  doc.appendChild( elemRoot );
  themes.writeXml( doc );

  QgsMapThemeCollection themes2( mProject );
  themes2.readXml( doc );

  QgsMapThemeCollection::MapThemeRecord recCollapsed2 = themes2.mapThemeState( "all-collapsed" );
  QVERIFY( recCollapsed2.hasExpandedStateInfo() );
  QCOMPARE( recCollapsed2.expandedGroupNodes().count(), 0 );
  QgsMapThemeCollection::MapThemeLayerRecord recCollapsedPoints2 = _recordForLayer( mPointsLayer, recCollapsed2 );
  QCOMPARE( recCollapsedPoints2.expandedLayerNode, false );
  QVERIFY( recCollapsedPoints2.expandedLegendItems.isEmpty() );

  QgsMapThemeCollection::MapThemeRecord recExpanded2 = themes2.mapThemeState( "all-expanded" );
  QVERIFY( recExpanded2.hasExpandedStateInfo() );
  QCOMPARE( recExpanded2.expandedGroupNodes().count(), 3 );
  QgsMapThemeCollection::MapThemeLayerRecord recExpandedPoints2 = _recordForLayer( mPointsLayer, recExpanded2 );
  QCOMPARE( recExpandedPoints2.expandedLayerNode, true );
  QCOMPARE( recExpandedPoints2.expandedLegendItems.count(), 4 );
}

QGSTEST_MAIN( TestQgsMapThemeCollection )
#include "testqgsmapthemecollection.moc"
