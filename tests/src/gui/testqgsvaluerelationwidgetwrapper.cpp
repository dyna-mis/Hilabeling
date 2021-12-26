/***************************************************************************
    testqgsvaluerelationwidgetwrapper.cpp
     --------------------------------------
    Date                 : 21 07 2017
    Copyright            : (C) 2017 Paul Blottiere
    Email                : paul dot blottiere at oslandia dot com
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#include "qgstest.h"
#include <QScrollBar>

#include <editorwidgets/core/qgseditorwidgetregistry.h>
#include <qgsapplication.h>
#include <qgsproject.h>
#include <qgsvectorlayer.h>
#include "qgseditorwidgetwrapper.h"
#include <editorwidgets/qgsvaluerelationwidgetwrapper.h>
#include <QTableWidget>
#include <QComboBox>
#include "qgsgui.h"
#include <gdal_version.h>

class TestQgsValueRelationWidgetWrapper : public QObject
{
    Q_OBJECT
  public:
    TestQgsValueRelationWidgetWrapper() = default;

  private:
    QTemporaryDir tempDir;

  private slots:
    void initTestCase(); // will be called before the first testfunction is executed.
    void cleanupTestCase(); // will be called after the last testfunction was executed.
    void init(); // will be called before each testfunction is executed.
    void cleanup(); // will be called after every testfunction.

    void testScrollBarUnlocked();
    void testDrillDown();
    void testDrillDownMulti();
    //! Checks that a related value of 0 is not interpreted as a NULL
    void testZeroIndexInRelatedTable();
    void testWithJsonInPostgres();
    void testWithJsonInGPKG();
};

void TestQgsValueRelationWidgetWrapper::initTestCase()
{
  QgsApplication::init();
  QgsApplication::initQgis();
  QgsGui::editorWidgetRegistry()->initEditors();
}

void TestQgsValueRelationWidgetWrapper::cleanupTestCase()
{
  QgsApplication::exitQgis();
}

void TestQgsValueRelationWidgetWrapper::init()
{
}

void TestQgsValueRelationWidgetWrapper::cleanup()
{

}

void TestQgsValueRelationWidgetWrapper::testScrollBarUnlocked()
{
  // create a vector layer
  QgsVectorLayer vl1( QStringLiteral( "LineString?crs=epsg:3111&field=pk:int&field=fk|:int" ), QStringLiteral( "vl1" ), QStringLiteral( "memory" ) );
  QgsProject::instance()->addMapLayer( &vl1, false, false );

  // build a value relation widget wrapper
  QgsValueRelationWidgetWrapper w( &vl1, 0, nullptr, nullptr );

  QVariantMap config;
  config.insert( QStringLiteral( "AllowMulti" ), true );
  w.setConfig( config );
  w.widget();
  w.setEnabled( true );

  // add an item virtually
  QTableWidgetItem item;
  item.setText( QStringLiteral( "MyText" ) );
  w.mTableWidget->setItem( 0, 0, &item );

  QCOMPARE( w.mTableWidget->item( 0, 0 )->text(), QString( "MyText" ) );

  // when the widget wrapper is enabled, the container should be enabled
  // as well as items
  w.setEnabled( true );

  QCOMPARE( w.widget()->isEnabled(), true );

  bool itemEnabled = w.mTableWidget->item( 0, 0 )->flags() & Qt::ItemIsEnabled;
  QCOMPARE( itemEnabled, true );

  // when the widget wrapper is disabled, the container should still be enabled
  // to keep the scrollbar available but items should be disabled to avoid
  // edition
  w.setEnabled( false );

  itemEnabled = w.mTableWidget->item( 0, 0 )->flags() & Qt::ItemIsEnabled;
  QCOMPARE( itemEnabled, false );

  QCOMPARE( w.widget()->isEnabled(), true );

  // recheck after re-enabled
  w.setEnabled( true );

  QCOMPARE( w.widget()->isEnabled(), true );
  itemEnabled = w.mTableWidget->item( 0, 0 )->flags() & Qt::ItemIsEnabled;
  QCOMPARE( itemEnabled, true );
}

void TestQgsValueRelationWidgetWrapper::testDrillDown()
{
  // create a vector layer
  QgsVectorLayer vl1( QStringLiteral( "Polygon?crs=epsg:4326&field=pk:int&field=province:int&field=municipality:string" ), QStringLiteral( "vl1" ), QStringLiteral( "memory" ) );
  QgsVectorLayer vl2( QStringLiteral( "Point?crs=epsg:4326&field=pk:int&field=fk_province:int&field=fk_municipality:int" ), QStringLiteral( "vl2" ), QStringLiteral( "memory" ) );
  QgsProject::instance()->addMapLayer( &vl1, false, false );
  QgsProject::instance()->addMapLayer( &vl2, false, false );

  // insert some features
  QgsFeature f1( vl1.fields() );
  f1.setAttribute( QStringLiteral( "pk" ), 1 );
  f1.setAttribute( QStringLiteral( "province" ), 123 );
  f1.setAttribute( QStringLiteral( "municipality" ), QStringLiteral( "Some Place By The River" ) );
  f1.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POLYGON(( 0 0, 0 1, 1 1, 1 0, 0 0 ))" ) ) );
  QVERIFY( f1.isValid() );
  QgsFeature f2( vl1.fields() );
  f2.setAttribute( QStringLiteral( "pk" ), 2 );
  f2.setAttribute( QStringLiteral( "province" ), 245 );
  f2.setAttribute( QStringLiteral( "municipality" ), QStringLiteral( "Dreamland By The Clouds" ) );
  f2.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POLYGON(( 1 0, 1 1, 2 1, 2 0, 1 0 ))" ) ) );
  QVERIFY( f2.isValid() );
  QVERIFY( vl1.dataProvider()->addFeatures( QgsFeatureList() << f1 << f2 ) );

  QgsFeature f3( vl2.fields() );
  f3.setAttribute( QStringLiteral( "fk_province" ), 123 );
  f3.setAttribute( QStringLiteral( "fk_municipality" ), 1 );
  f3.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POINT( 0.5 0.5)" ) ) );
  QVERIFY( f3.isValid() );
  QVERIFY( f3.geometry().isGeosValid() );
  QVERIFY( vl2.dataProvider()->addFeature( f3 ) );

  // build a value relation widget wrapper for municipality
  QgsValueRelationWidgetWrapper w_municipality( &vl2, vl2.fields().indexOf( QStringLiteral( "fk_municipality" ) ), nullptr, nullptr );
  QVariantMap cfg_municipality;
  cfg_municipality.insert( QStringLiteral( "Layer" ), vl1.id() );
  cfg_municipality.insert( QStringLiteral( "Key" ),  QStringLiteral( "pk" ) );
  cfg_municipality.insert( QStringLiteral( "Value" ), QStringLiteral( "municipality" ) );
  cfg_municipality.insert( QStringLiteral( "AllowMulti" ), false );
  cfg_municipality.insert( QStringLiteral( "NofColumns" ), 1 );
  cfg_municipality.insert( QStringLiteral( "AllowNull" ), false );
  cfg_municipality.insert( QStringLiteral( "OrderByValue" ), true );
  cfg_municipality.insert( QStringLiteral( "FilterExpression" ), QStringLiteral( "\"province\" =  current_value('fk_province')" ) );
  cfg_municipality.insert( QStringLiteral( "UseCompleter" ), false );
  w_municipality.setConfig( cfg_municipality );
  w_municipality.widget();
  w_municipality.setEnabled( true );

  QCOMPARE( w_municipality.mCache.size(), 2 );
  QCOMPARE( w_municipality.mComboBox->count(), 2 );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mCache.size(), 1 );

  // Check first is selected
  QCOMPARE( w_municipality.mComboBox->count(), 1 );
  QCOMPARE( w_municipality.mComboBox->itemText( 0 ), QStringLiteral( "Some Place By The River" ) );
  QCOMPARE( w_municipality.value().toString(), QStringLiteral( "1" ) );

  // Filter by geometry
  cfg_municipality[ QStringLiteral( "FilterExpression" ) ] = QStringLiteral( "contains(buffer(@current_geometry, 1 ), $geometry)" );
  w_municipality.setConfig( cfg_municipality );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mComboBox->count(), 1 );
  QCOMPARE( w_municipality.mComboBox->itemText( 0 ), QStringLiteral( "Some Place By The River" ) );

  // Move the point to 1.5 0.5
  f3.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POINT( 1.5 0.5)" ) ) );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mComboBox->count(), 1 );
  QCOMPARE( w_municipality.mComboBox->itemText( 0 ), QStringLiteral( "Dreamland By The Clouds" ) );

  // Enlarge the buffer
  cfg_municipality[ QStringLiteral( "FilterExpression" ) ] = QStringLiteral( "contains(buffer(@current_geometry, 3 ), $geometry)" );
  w_municipality.setConfig( cfg_municipality );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mComboBox->count(), 2 );
  QCOMPARE( w_municipality.mComboBox->itemText( 0 ), QStringLiteral( "Dreamland By The Clouds" ) );
  QCOMPARE( w_municipality.mComboBox->itemText( 1 ), QStringLiteral( "Some Place By The River" ) );

  // Check with allow null
  cfg_municipality[QStringLiteral( "AllowNull" )] = true;
  w_municipality.setConfig( cfg_municipality );
  w_municipality.setFeature( QgsFeature() );

  // Check null is selected
  QCOMPARE( w_municipality.mComboBox->count(), 3 );
  QCOMPARE( w_municipality.mComboBox->itemText( 0 ), QStringLiteral( "(no selection)" ) );
  QVERIFY( w_municipality.value().isNull() );
  QCOMPARE( w_municipality.value().toString(), QString() );

  // Check order by value false
  cfg_municipality[QStringLiteral( "AllowNull" )] = false;
  cfg_municipality[QStringLiteral( "OrderByValue" )] = false;
  w_municipality.setConfig( cfg_municipality );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mComboBox->itemText( 1 ), QStringLiteral( "Dreamland By The Clouds" ) );
  QCOMPARE( w_municipality.mComboBox->itemText( 0 ), QStringLiteral( "Some Place By The River" ) );

}


void TestQgsValueRelationWidgetWrapper::testDrillDownMulti()
{
  // create a vector layer
  QgsVectorLayer vl1( QStringLiteral( "Polygon?crs=epsg:4326&field=pk:int&field=province:int&field=municipality:string" ), QStringLiteral( "vl1" ), QStringLiteral( "memory" ) );
  QgsVectorLayer vl2( QStringLiteral( "Point?crs=epsg:4326&field=pk:int&field=fk_province:int&field=fk_municipality:int" ), QStringLiteral( "vl2" ), QStringLiteral( "memory" ) );
  QgsProject::instance()->addMapLayer( &vl1, false, false );
  QgsProject::instance()->addMapLayer( &vl2, false, false );

  // insert some features
  QgsFeature f1( vl1.fields() );
  f1.setAttribute( QStringLiteral( "pk" ), 1 );
  f1.setAttribute( QStringLiteral( "province" ), 123 );
  f1.setAttribute( QStringLiteral( "municipality" ), QStringLiteral( "Some Place By The River" ) );
  f1.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POLYGON(( 0 0, 0 1, 1 1, 1 0, 0 0 ))" ) ) );
  QVERIFY( f1.isValid() );
  QgsFeature f2( vl1.fields() );
  f2.setAttribute( QStringLiteral( "pk" ), 2 );
  f2.setAttribute( QStringLiteral( "province" ), 245 );
  f2.setAttribute( QStringLiteral( "municipality" ), QStringLiteral( "Dreamland By The Clouds" ) );
  f2.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POLYGON(( 1 0, 1 1, 2 1, 2 0, 1 0 ))" ) ) );
  QVERIFY( f2.isValid() );
  QVERIFY( vl1.dataProvider()->addFeatures( QgsFeatureList() << f1 << f2 ) );

  QgsFeature f3( vl2.fields() );
  f3.setAttribute( QStringLiteral( "fk_province" ), 123 );
  f3.setAttribute( QStringLiteral( "fk_municipality" ), QStringLiteral( "{1}" ) );
  f3.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POINT( 0.5 0.5)" ) ) );
  QVERIFY( f3.isValid() );
  QVERIFY( f3.geometry().isGeosValid() );
  QVERIFY( vl2.dataProvider()->addFeature( f3 ) );

  // build a value relation widget wrapper for municipality
  QgsValueRelationWidgetWrapper w_municipality( &vl2, vl2.fields().indexOf( QStringLiteral( "fk_municipality" ) ), nullptr, nullptr );
  QVariantMap cfg_municipality;
  cfg_municipality.insert( QStringLiteral( "Layer" ), vl1.id() );
  cfg_municipality.insert( QStringLiteral( "Key" ),  QStringLiteral( "pk" ) );
  cfg_municipality.insert( QStringLiteral( "Value" ), QStringLiteral( "municipality" ) );
  cfg_municipality.insert( QStringLiteral( "AllowMulti" ), true );
  cfg_municipality.insert( QStringLiteral( "NofColumns" ), 1 );
  cfg_municipality.insert( QStringLiteral( "AllowNull" ), false );
  cfg_municipality.insert( QStringLiteral( "OrderByValue" ), true );
  cfg_municipality.insert( QStringLiteral( "FilterExpression" ), QStringLiteral( "\"province\" =  current_value('fk_province')" ) );
  cfg_municipality.insert( QStringLiteral( "UseCompleter" ), false );
  w_municipality.setConfig( cfg_municipality );
  w_municipality.widget();
  w_municipality.setEnabled( true );

  QCOMPARE( w_municipality.mCache.size(), 2 );
  QCOMPARE( w_municipality.mTableWidget->rowCount(), 2 );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mCache.size(), 1 );

  QCOMPARE( w_municipality.mTableWidget->rowCount(), 1 );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Some Place By The River" ) );
  QCOMPARE( w_municipality.value().toString(), QStringLiteral( "{1}" ) );

  // Filter by geometry
  cfg_municipality[ QStringLiteral( "FilterExpression" ) ] = QStringLiteral( "contains(buffer(@current_geometry, 1 ), $geometry)" );
  w_municipality.setConfig( cfg_municipality );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mTableWidget->rowCount(), 1 );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Some Place By The River" ) );

  // Move the point to 1.5 0.5
  f3.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POINT( 1.5 0.5)" ) ) );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mTableWidget->rowCount(), 1 );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Dreamland By The Clouds" ) );

  // Enlarge the buffer
  cfg_municipality[ QStringLiteral( "FilterExpression" ) ] = QStringLiteral( "contains(buffer(@current_geometry, 3 ), $geometry)" );
  w_municipality.setConfig( cfg_municipality );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mTableWidget->rowCount(), 2 );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Dreamland By The Clouds" ) );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "2" ) );
  QCOMPARE( w_municipality.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Some Place By The River" ) );
  QCOMPARE( w_municipality.mTableWidget->item( 1, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "1" ) );
  QCOMPARE( w_municipality.value().toString(), QStringLiteral( "{1}" ) );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_municipality.mTableWidget->item( 1, 0 )->checkState(), Qt::Checked );
  w_municipality.setValue( QStringLiteral( "{1,2}" ) );
  QCOMPARE( w_municipality.value().toString(), QStringLiteral( "{2,1}" ) );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_municipality.mTableWidget->item( 1, 0 )->checkState(), Qt::Checked );

  // Check values are checked
  f3.setAttribute( QStringLiteral( "fk_municipality" ), QStringLiteral( "{1,2}" ) );
  w_municipality.setFeature( f3 );
  QCOMPARE( w_municipality.mTableWidget->rowCount(), 2 );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Dreamland By The Clouds" ) );
  QCOMPARE( w_municipality.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Some Place By The River" ) );
  QCOMPARE( w_municipality.mTableWidget->item( 0, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_municipality.mTableWidget->item( 1, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_municipality.value().toString(), QStringLiteral( "{2,1}" ) );

}

void TestQgsValueRelationWidgetWrapper::testZeroIndexInRelatedTable()
{
  // findData fails to tell a 0 from a NULL
  // See: "Value relation, value 0 = NULL" - https://issues.qgis.org/issues/19981

  // create a vector layer
  QgsVectorLayer vl1( QStringLiteral( "Polygon?crs=epsg:4326&field=pk:int&field=province:int&field=municipality:string" ), QStringLiteral( "vl1" ), QStringLiteral( "memory" ) );
  QgsVectorLayer vl2( QStringLiteral( "Point?crs=epsg:4326&field=pk:int&field=fk_province:int&field=fk_municipality:int" ), QStringLiteral( "vl2" ), QStringLiteral( "memory" ) );
  QgsProject::instance()->addMapLayer( &vl1, false, false );
  QgsProject::instance()->addMapLayer( &vl2, false, false );

  // insert some features
  QgsFeature f1( vl1.fields() );
  f1.setAttribute( QStringLiteral( "pk" ), 0 );  // !!! Notice: pk 0
  f1.setAttribute( QStringLiteral( "province" ), 123 );
  f1.setAttribute( QStringLiteral( "municipality" ), QStringLiteral( "Some Place By The River" ) );
  f1.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POLYGON(( 0 0, 0 1, 1 1, 1 0, 0 0 ))" ) ) );
  QVERIFY( f1.isValid() );
  QgsFeature f2( vl1.fields() );
  f2.setAttribute( QStringLiteral( "pk" ), 2 );
  f2.setAttribute( QStringLiteral( "province" ), 245 );
  f2.setAttribute( QStringLiteral( "municipality" ), QStringLiteral( "Dreamland By The Clouds" ) );
  f2.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POLYGON(( 1 0, 1 1, 2 1, 2 0, 1 0 ))" ) ) );
  QVERIFY( f2.isValid() );
  QVERIFY( vl1.dataProvider()->addFeatures( QgsFeatureList() << f1 << f2 ) );

  QgsFeature f3( vl2.fields() );
  f3.setAttribute( QStringLiteral( "fk_province" ), 123 );
  f3.setAttribute( QStringLiteral( "fk_municipality" ), QStringLiteral( "0" ) );
  f3.setGeometry( QgsGeometry::fromWkt( QStringLiteral( "POINT( 0.5 0.5)" ) ) );
  QVERIFY( f3.isValid() );
  QVERIFY( f3.geometry().isGeosValid() );
  QVERIFY( vl2.dataProvider()->addFeature( f3 ) );

  // build a value relation widget wrapper for municipality
  QgsValueRelationWidgetWrapper w_municipality( &vl2, vl2.fields().indexOf( QStringLiteral( "fk_municipality" ) ), nullptr, nullptr );
  QVariantMap cfg_municipality;
  cfg_municipality.insert( QStringLiteral( "Layer" ), vl1.id() );
  cfg_municipality.insert( QStringLiteral( "Key" ),  QStringLiteral( "pk" ) );
  cfg_municipality.insert( QStringLiteral( "Value" ), QStringLiteral( "municipality" ) );
  cfg_municipality.insert( QStringLiteral( "AllowMulti" ), false );
  cfg_municipality.insert( QStringLiteral( "NofColumns" ), 1 );
  cfg_municipality.insert( QStringLiteral( "AllowNull" ), true );
  cfg_municipality.insert( QStringLiteral( "OrderByValue" ), false );
  cfg_municipality.insert( QStringLiteral( "UseCompleter" ), false );
  w_municipality.setConfig( cfg_municipality );
  w_municipality.widget();
  w_municipality.setEnabled( true );

  w_municipality.setValue( 0 );
  QCOMPARE( w_municipality.mComboBox->currentIndex(), 1 );
  QCOMPARE( w_municipality.mComboBox->currentText(), QStringLiteral( "Some Place By The River" ) );
}

void TestQgsValueRelationWidgetWrapper::testWithJsonInPostgres()
{
  //this is only reading

  // create pg layers
  QString dbConn = getenv( "QGIS_PGTEST_DB" );
  if ( dbConn.isEmpty() )
  {
    dbConn = "service=\"qgis_test\"";
  }
  QgsVectorLayer *vl_json = new QgsVectorLayer( QStringLiteral( "%1 sslmode=disable key=\"pk\" table=\"qgis_test\".\"json\" sql=" ).arg( dbConn ), QStringLiteral( "json" ), QStringLiteral( "postgres" ) );
  QgsVectorLayer *vl_authors = new QgsVectorLayer( QStringLiteral( "%1 sslmode=disable key='pk' table=\"qgis_test\".\"authors\" sql=" ).arg( dbConn ), QStringLiteral( "authors" ), QStringLiteral( "postgres" ) );
  QVERIFY( vl_json->isValid() );
  QVERIFY( vl_authors->isValid() );

  QgsProject::instance()->addMapLayer( vl_json, false, false );
  QgsProject::instance()->addMapLayer( vl_authors, false, false );

  QCOMPARE( vl_json->fields().at( 1 ).type(), QVariant::Map );

  // build a value relation widget wrapper for json field
  QgsValueRelationWidgetWrapper w_favoriteauthors( vl_json, vl_json->fields().indexOf( QStringLiteral( "jvalue" ) ), nullptr, nullptr );
  QVariantMap cfg_favoriteautors;
  cfg_favoriteautors.insert( QStringLiteral( "Layer" ), vl_authors->id() );
  cfg_favoriteautors.insert( QStringLiteral( "Key" ),  QStringLiteral( "pk" ) );
  cfg_favoriteautors.insert( QStringLiteral( "Value" ), QStringLiteral( "name" ) );
  cfg_favoriteautors.insert( QStringLiteral( "AllowMulti" ), true );
  cfg_favoriteautors.insert( QStringLiteral( "NofColumns" ), 1 );
  cfg_favoriteautors.insert( QStringLiteral( "AllowNull" ), false );
  cfg_favoriteautors.insert( QStringLiteral( "OrderByValue" ), false );
  cfg_favoriteautors.insert( QStringLiteral( "UseCompleter" ), false );
  w_favoriteauthors.setConfig( cfg_favoriteautors );
  w_favoriteauthors.widget();
  w_favoriteauthors.setEnabled( true );

  //check if set up nice
  QCOMPARE( w_favoriteauthors.mTableWidget->rowCount(), 7 );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "1" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "2" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "3" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "4" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "5" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "6" ) );

  //check if first feature checked correctly (should be 1,2,3 and the rest is not)
  w_favoriteauthors.setFeature( vl_json->getFeature( 1 ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->checkState(), Qt::Unchecked );

  // build a value relation widget wrapper for jsonb field
  QgsValueRelationWidgetWrapper w_favoriteauthors_b( vl_json, vl_json->fields().indexOf( QStringLiteral( "jbvalue" ) ), nullptr, nullptr );
  QVariantMap cfg_favoriteautors_b;
  cfg_favoriteautors_b.insert( QStringLiteral( "Layer" ), vl_authors->id() );
  cfg_favoriteautors_b.insert( QStringLiteral( "Key" ),  QStringLiteral( "pk" ) );
  cfg_favoriteautors_b.insert( QStringLiteral( "Value" ), QStringLiteral( "name" ) );
  cfg_favoriteautors_b.insert( QStringLiteral( "AllowMulti" ), true );
  cfg_favoriteautors_b.insert( QStringLiteral( "NofColumns" ), 1 );
  cfg_favoriteautors_b.insert( QStringLiteral( "AllowNull" ), false );
  cfg_favoriteautors_b.insert( QStringLiteral( "OrderByValue" ), false );
  cfg_favoriteautors_b.insert( QStringLiteral( "UseCompleter" ), false );
  w_favoriteauthors_b.setConfig( cfg_favoriteautors_b );
  w_favoriteauthors_b.widget();
  w_favoriteauthors_b.setEnabled( true );

  //check if set up nice
  QCOMPARE( w_favoriteauthors_b.mTableWidget->rowCount(), 7 );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 0, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "1" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 1, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "2" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 2, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "3" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 3, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "4" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 4, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "5" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 5, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "6" ) );

  //check if second feature checked correctly (should be 4,5,6 and the rest is not)
  w_favoriteauthors_b.setFeature( vl_json->getFeature( 1 ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 0, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 1, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 2, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 3, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 4, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors_b.mTableWidget->item( 5, 0 )->checkState(), Qt::Checked );
}


void TestQgsValueRelationWidgetWrapper::testWithJsonInGPKG()
{
#if GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(2,4,0)
  // create ogr gpkg layers
  QString myFileName( TEST_DATA_DIR ); //defined in CmakeLists.txt
  QString myTempDirName = tempDir.path();
  QFile::copy( myFileName + "/provider/test_json.gpkg", myTempDirName + "/test_json.gpkg" );
  QString myTempFileName = myTempDirName + "/test_json.gpkg";
  QFileInfo myMapFileInfo( myTempFileName );
  QgsVectorLayer *vl_json = new QgsVectorLayer( myMapFileInfo.filePath() + "|layername=foo", "test", QStringLiteral( "ogr" ) );
  QgsVectorLayer *vl_authors = new QgsVectorLayer( myMapFileInfo.filePath() + "|layername=author", "test", QStringLiteral( "ogr" ) );
  QVERIFY( vl_json->isValid() );
  QVERIFY( vl_authors->isValid() );

  QgsProject::instance()->addMapLayer( vl_json, false, false );
  QgsProject::instance()->addMapLayer( vl_authors, false, false );
  vl_json->startEditing();

  // build a value relation widget wrapper for authors
  QgsValueRelationWidgetWrapper w_favoriteauthors( vl_json, vl_json->fields().indexOf( QStringLiteral( "json_content" ) ), nullptr, nullptr );
  QVariantMap cfg_favoriteautors;
  cfg_favoriteautors.insert( QStringLiteral( "Layer" ), vl_authors->id() );
  cfg_favoriteautors.insert( QStringLiteral( "Key" ),  QStringLiteral( "fid" ) );
  cfg_favoriteautors.insert( QStringLiteral( "Value" ), QStringLiteral( "NAME" ) );
  cfg_favoriteautors.insert( QStringLiteral( "AllowMulti" ), true );
  cfg_favoriteautors.insert( QStringLiteral( "NofColumns" ), 1 );
  cfg_favoriteautors.insert( QStringLiteral( "AllowNull" ), false );
  cfg_favoriteautors.insert( QStringLiteral( "OrderByValue" ), false );
  cfg_favoriteautors.insert( QStringLiteral( "UseCompleter" ), false );
  w_favoriteauthors.setConfig( cfg_favoriteautors );
  w_favoriteauthors.widget();
  w_favoriteauthors.setEnabled( true );

  //check if set up nice
  QCOMPARE( w_favoriteauthors.mTableWidget->rowCount(), 6 );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "1" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "2" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "3" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "4" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "5" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->data( Qt::UserRole ).toString(), QStringLiteral( "6" ) );

  w_favoriteauthors.setFeature( vl_json->getFeature( 1 ) );

  //check if first feature checked correctly (none)
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->checkState(), Qt::Unchecked );

  //check other authors
  w_favoriteauthors.mTableWidget->item( 0, 0 )->setCheckState( Qt::Checked );
  w_favoriteauthors.mTableWidget->item( 2, 0 )->setCheckState( Qt::Checked );
  w_favoriteauthors.mTableWidget->item( 4, 0 )->setCheckState( Qt::Checked );

  //check if first feature checked correctly 0, 2, 4 (means the fids 1, 3, 5
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->checkState(), Qt::Unchecked );

  //we do jump over the part with QgsAttributeForm::saveEdits
  vl_json->changeAttributeValue( 1, 4, w_favoriteauthors.value() );

  w_favoriteauthors.setFeature( vl_json->getFeature( 2 ) );
  //check if second feature checked correctly (none)
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->checkState(), Qt::Unchecked );

  w_favoriteauthors.setFeature( vl_json->getFeature( 1 ) );
  //check if first feature checked correctly 0, 2, 4
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->text(), QStringLiteral( "Erich Gamma" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 0, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->text(), QStringLiteral( "Richard Helm" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 1, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->text(), QStringLiteral( "Ralph Johnson" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 2, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->text(), QStringLiteral( "John Vlissides" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 3, 0 )->checkState(), Qt::Unchecked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->text(), QStringLiteral( "Douglas Adams" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 4, 0 )->checkState(), Qt::Checked );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->text(), QStringLiteral( "Ken Follett" ) );
  QCOMPARE( w_favoriteauthors.mTableWidget->item( 5, 0 )->checkState(), Qt::Unchecked );

  // check if stored correctly
  vl_json->commitChanges();
  QVariantList expected_vl;
  expected_vl << "1" << "3" << "5";

  QgsFeature f = vl_json->getFeature( 1 );
  QVariant attribute = f.attribute( QStringLiteral( "json_content" ) );
  QList<QVariant> value = attribute.toList();
  QCOMPARE( value, expected_vl );
#endif
}


QGSTEST_MAIN( TestQgsValueRelationWidgetWrapper )
#include "testqgsvaluerelationwidgetwrapper.moc"
