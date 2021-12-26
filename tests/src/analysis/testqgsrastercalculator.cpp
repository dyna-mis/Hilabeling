/***************************************************************************
  testqgsrastercalculator.cpp
  --------------------------------------
Date                 : Jun-2015
Copyright            : (C) 2015 by Nyall Dawson
Email                : nyall dot dawson at gmail dot com
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "qgstest.h"

#ifdef HAVE_OPENCL
#include "qgsopenclutils.h"
#endif

#include "qgsrastercalculator.h"
#include "qgsrastercalcnode.h"
#include "qgsrasterdataprovider.h"
#include "qgsrasterlayer.h"
#include "qgsrastermatrix.h"
#include "qgsapplication.h"
#include "qgsproject.h"

Q_DECLARE_METATYPE( QgsRasterCalcNode::Operator )

class TestQgsRasterCalculator : public QObject
{
    Q_OBJECT

  public:
    TestQgsRasterCalculator() = default;

  private slots:
    void initTestCase();// will be called before the first testfunction is executed.
    void cleanupTestCase();// will be called after the last testfunction was executed.
    void init() ;// will be called before each testfunction is executed.
    void cleanup() ;// will be called after every testfunction.

    void dualOp_data();
    void dualOp(); //test operators which operate on a left&right node

    void singleOp_data();
    void singleOp(); //test operators which operate on a single value

    void singleOpMatrices(); // test single op using matrix
    void dualOpNumberMatrix(); // test dual op run on number and matrix
    void dualOpMatrixNumber(); // test dual op run on matrix and number
    void dualOpMatrixMatrix(); // test dual op run on matrix and matrix

    void rasterRefOp();
    void dualOpRasterRaster(); //test dual op on raster ref and raster ref

    void calcWithLayers();
    void calcWithReprojectedLayers();

    void errors();
    void toString();
    void findNodes();

    void testRasterEntries();
    void calcFormulasWithReprojectedLayers();

  private:

    QgsRasterLayer *mpLandsatRasterLayer = nullptr;
    QgsRasterLayer *mpLandsatRasterLayer4326 = nullptr;
};


void  TestQgsRasterCalculator::initTestCase()
{
  //
  // Runs once before any tests are run
  //
  // Set up the QgsSettings environment
  QCoreApplication::setOrganizationName( QStringLiteral( "QGIS" ) );
  QCoreApplication::setOrganizationDomain( QStringLiteral( "qgis.org" ) );
  QCoreApplication::setApplicationName( QStringLiteral( "QGIS-TEST" ) );

  QgsApplication::init();
  QgsApplication::initQgis();


  QString testDataDir = QStringLiteral( TEST_DATA_DIR ) + '/'; //defined in CmakeLists.txt

  QString landsatFileName = testDataDir + "landsat.tif";
  QFileInfo landsatRasterFileInfo( landsatFileName );
  mpLandsatRasterLayer = new QgsRasterLayer( landsatRasterFileInfo.filePath(),
      landsatRasterFileInfo.completeBaseName() );


  QString landsat4326FileName = testDataDir + "landsat_4326.tif";
  QFileInfo landsat4326RasterFileInfo( landsat4326FileName );
  mpLandsatRasterLayer4326 = new QgsRasterLayer( landsat4326RasterFileInfo.filePath(),
      landsat4326RasterFileInfo.completeBaseName() );

  QgsProject::instance()->addMapLayers(
    QList<QgsMapLayer *>() << mpLandsatRasterLayer << mpLandsatRasterLayer4326 );
}

void  TestQgsRasterCalculator::cleanupTestCase()
{
  QgsApplication::exitQgis();
}

void  TestQgsRasterCalculator::init()
{
#ifdef HAVE_OPENCL
  QgsOpenClUtils::setEnabled( false );
  // Reset to default in case some tests mess it up
  QgsOpenClUtils::setSourcePath( QDir( QgsApplication::pkgDataPath() ).absoluteFilePath( QStringLiteral( "resources/opencl_programs" ) ) );
#endif
}

void  TestQgsRasterCalculator::cleanup()
{

}

void TestQgsRasterCalculator::dualOp_data()
{
  QTest::addColumn< QgsRasterCalcNode::Operator >( "op" );
  QTest::addColumn<double>( "left" );
  QTest::addColumn<double>( "right" );
  QTest::addColumn<double>( "expected" );

  QTest::newRow( "opPlus" ) << QgsRasterCalcNode::opPLUS << 5.5 << 2.2 << 7.7;
  QTest::newRow( "opMINUS" ) << QgsRasterCalcNode::opMINUS << 5.0 << 2.5 << 2.5;
  QTest::newRow( "opMUL" ) << QgsRasterCalcNode::opMUL << 2.5 << 4.0 << 10.0;
  QTest::newRow( "opDIV" ) << QgsRasterCalcNode::opDIV << 2.5 << 2.0 << 1.25;
  QTest::newRow( "opDIV by 0" ) << QgsRasterCalcNode::opDIV << 2.5 << 0.0 << -9999.0;
  QTest::newRow( "opPOW" ) << QgsRasterCalcNode::opPOW << 3.0 << 2.0 << 9.0;
  QTest::newRow( "opPOW negative" ) << QgsRasterCalcNode::opPOW << 4.0 << -2.0 << 0.0625;
  QTest::newRow( "opPOW sqrt" ) << QgsRasterCalcNode::opPOW << 4.0 << 0.5 << 2.0;
  QTest::newRow( "opPOW complex" ) << QgsRasterCalcNode::opPOW << -2.0 << 0.5 << -9999.0;
  QTest::newRow( "opEQ true" ) << QgsRasterCalcNode::opEQ << 1.0 << 1.0 << 1.0;
  QTest::newRow( "opEQ false" ) << QgsRasterCalcNode::opEQ << 0.5 << 1.0 << 0.0;
  QTest::newRow( "opNE equal" ) << QgsRasterCalcNode::opNE << 1.0 << 1.0 << 0.0;
  QTest::newRow( "opNE not equal" ) << QgsRasterCalcNode::opNE << 0.5 << 1.0 << 1.0;
  QTest::newRow( "opGT >" ) << QgsRasterCalcNode::opGT << 1.0 << 0.5 << 1.0;
  QTest::newRow( "opGT =" ) << QgsRasterCalcNode::opGT << 0.5 << 0.5 << 0.0;
  QTest::newRow( "opGT <" ) << QgsRasterCalcNode::opGT << 0.5 << 1.0 << 0.0;
  QTest::newRow( "opLT >" ) << QgsRasterCalcNode::opLT << 1.0 << 0.5 << 0.0;
  QTest::newRow( "opLT =" ) << QgsRasterCalcNode::opLT << 0.5 << 0.5 << 0.0;
  QTest::newRow( "opLT <" ) << QgsRasterCalcNode::opLT << 0.5 << 1.0 << 1.0;
  QTest::newRow( "opGE >" ) << QgsRasterCalcNode::opGE << 1.0 << 0.5 << 1.0;
  QTest::newRow( "opGE =" ) << QgsRasterCalcNode::opGE << 0.5 << 0.5 << 1.0;
  QTest::newRow( "opGE <" ) << QgsRasterCalcNode::opGE << 0.5 << 1.0 << 0.0;
  QTest::newRow( "opLE >" ) << QgsRasterCalcNode::opLE << 1.0 << 0.5 << 0.0;
  QTest::newRow( "opLE =" ) << QgsRasterCalcNode::opLE << 0.5 << 0.5 << 1.0;
  QTest::newRow( "opLE <" ) << QgsRasterCalcNode::opLE << 0.5 << 1.0 << 1.0;
  QTest::newRow( "opAND 0/0" ) << QgsRasterCalcNode::opAND << 0.0 << 0.0 << 0.0;
  QTest::newRow( "opAND 0/1" ) << QgsRasterCalcNode::opAND << 0.0 << 1.0 << 0.0;
  QTest::newRow( "opAND 1/0" ) << QgsRasterCalcNode::opAND << 1.0 << 0.0 << 0.0;
  QTest::newRow( "opAND 1/1" ) << QgsRasterCalcNode::opAND << 1.0 << 1.0 << 1.0;
  QTest::newRow( "opOR 0/0" ) << QgsRasterCalcNode::opOR << 0.0 << 0.0 << 0.0;
  QTest::newRow( "opOR 0/1" ) << QgsRasterCalcNode::opOR << 0.0 << 1.0 << 1.0;
  QTest::newRow( "opOR 1/0" ) << QgsRasterCalcNode::opOR << 1.0 << 0.0 << 1.0;
  QTest::newRow( "opOR 1/1" ) << QgsRasterCalcNode::opOR << 1.0 << 1.0 << 1.0;
}

void TestQgsRasterCalculator::dualOp()
{
  QFETCH( QgsRasterCalcNode::Operator, op );
  QFETCH( double, left );
  QFETCH( double, right );
  QFETCH( double, expected );

  QgsRasterCalcNode node( op, new QgsRasterCalcNode( left ), new QgsRasterCalcNode( right ) );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );
  QMap<QString, QgsRasterBlock *> rasterData;

  QVERIFY( node.calculate( rasterData, result ) );

  qDebug() << "Result: " << result.number() << " expected: " << expected;
  QCOMPARE( result.number(), expected );

}

void TestQgsRasterCalculator::singleOp_data()
{
  QTest::addColumn< QgsRasterCalcNode::Operator >( "op" );
  QTest::addColumn<double>( "value" );
  QTest::addColumn<double>( "expected" );

  QTest::newRow( "opSQRT" ) << QgsRasterCalcNode::opSQRT << 16.0 << 4.0;
  QTest::newRow( "opSQRT negative" ) << QgsRasterCalcNode::opSQRT << -16.0 << -9999.0;
  QTest::newRow( "opSIN 0" ) << QgsRasterCalcNode::opSIN << 0.0 << 0.0;
  QTest::newRow( "opSIN pi/2" ) << QgsRasterCalcNode::opSIN << M_PI / 2.0 << 1.0;
  QTest::newRow( "opCOS 0" ) << QgsRasterCalcNode::opCOS << 0.0 << 1.0;
  QTest::newRow( "opCOS pi" ) << QgsRasterCalcNode::opCOS << M_PI << -1.0;
  QTest::newRow( "opTAN 0" ) << QgsRasterCalcNode::opTAN << 0.0 << 0.0;
  QTest::newRow( "opTAN pi" ) << QgsRasterCalcNode::opTAN << M_PI << 0.0;
  QTest::newRow( "opASIN 0" ) << QgsRasterCalcNode::opASIN << 0.0 << 0.0;
  QTest::newRow( "opASIN pi/2" ) << QgsRasterCalcNode::opASIN << 1.0 << M_PI / 2.0;
  QTest::newRow( "opACOS 0" ) << QgsRasterCalcNode::opACOS << 1.0 << 0.0;
  QTest::newRow( "opACOS pi/2" ) << QgsRasterCalcNode::opACOS << -1.0 << M_PI;
  QTest::newRow( "opATAN 0" ) << QgsRasterCalcNode::opATAN << 0.0 << 0.0;
  QTest::newRow( "opATAN 1.0" ) << QgsRasterCalcNode::opATAN << 1.0 << 0.7853981634;
  QTest::newRow( "opSIGN +" ) << QgsRasterCalcNode::opSIGN << 1.0 << -1.0;
  QTest::newRow( "opSIGN -" ) << QgsRasterCalcNode::opSIGN << -1.0 << 1.0;
  QTest::newRow( "opLOG -1" ) << QgsRasterCalcNode::opLOG << -1.0 << -9999.0;
  QTest::newRow( "opLOG 0" ) << QgsRasterCalcNode::opLOG << 0.0 << -9999.0;
  QTest::newRow( "opLOG 1" ) << QgsRasterCalcNode::opLOG << 1.0 << 0.0;
  QTest::newRow( "opLOG10 -1" ) << QgsRasterCalcNode::opLOG10 << -1.0 << -9999.0;
  QTest::newRow( "opLOG10 0" ) << QgsRasterCalcNode::opLOG10 << 0.0 << -9999.0;
  QTest::newRow( "opLOG10 1" ) << QgsRasterCalcNode::opLOG10 << 1.0 << 0.0;
  QTest::newRow( "opLOG10 10" ) << QgsRasterCalcNode::opLOG10 << 10.0 << 1.0;
}

void TestQgsRasterCalculator::singleOp()
{
  QFETCH( QgsRasterCalcNode::Operator, op );
  QFETCH( double, value );
  QFETCH( double, expected );

  QgsRasterCalcNode node( op, new QgsRasterCalcNode( value ), nullptr );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );
  QMap<QString, QgsRasterBlock *> rasterData;

  QVERIFY( node.calculate( rasterData, result ) );

  qDebug() << "Result: " << result.number() << " expected: " << expected;
  QGSCOMPARENEAR( result.number(), expected, 0.0000000001 );

}

void TestQgsRasterCalculator::singleOpMatrices()
{
  // test single op run on matrix
  double *d = new double[6];
  d[0] = 1.0;
  d[1] = 2.0;
  d[2] = 3.0;
  d[3] = 4.0;
  d[4] = 5.0;
  d[5] = -1.0;

  QgsRasterMatrix m( 2, 3, d, -1.0 );

  QgsRasterCalcNode node( QgsRasterCalcNode::opSIGN, new QgsRasterCalcNode( &m ), nullptr );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );
  QMap<QString, QgsRasterBlock *> rasterData;

  QVERIFY( node.calculate( rasterData, result ) );

  QCOMPARE( result.data()[0], -d[0] );
  QCOMPARE( result.data()[1], -d[1] );
  QCOMPARE( result.data()[2], -d[2] );
  QCOMPARE( result.data()[3], -d[3] );
  QCOMPARE( result.data()[4], -d[4] );
  QCOMPARE( result.data()[5], -9999.0 );
}

void TestQgsRasterCalculator::dualOpNumberMatrix()
{
  // test dual op run on number and matrix
  double *d = new double[6];
  d[0] = 1.0;
  d[1] = 2.0;
  d[2] = 3.0;
  d[3] = 4.0;
  d[4] = 5.0;
  d[5] = -1.0;

  QgsRasterMatrix m( 2, 3, d, -1.0 );

  QgsRasterCalcNode node( QgsRasterCalcNode::opPLUS, new QgsRasterCalcNode( 5.0 ), new QgsRasterCalcNode( &m ) );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );
  QMap<QString, QgsRasterBlock *> rasterData;

  QVERIFY( node.calculate( rasterData, result ) );

  QCOMPARE( result.data()[0], 6.0 );
  QCOMPARE( result.data()[1], 7.0 );
  QCOMPARE( result.data()[2], 8.0 );
  QCOMPARE( result.data()[3], 9.0 );
  QCOMPARE( result.data()[4], 10.0 );
  QCOMPARE( result.data()[5], -9999.0 );

  //also check adding no data number
  QgsRasterCalcNode nodeNoData( QgsRasterCalcNode::opPLUS, new QgsRasterCalcNode( -9999 ), new QgsRasterCalcNode( &m ) );
  QVERIFY( nodeNoData.calculate( rasterData, result ) );
  QCOMPARE( result.data()[0], -9999.0 );
}

void TestQgsRasterCalculator::dualOpMatrixNumber()
{
  // test dual op run on matrix and number
  double *d = new double[6];
  d[0] = 1.0;
  d[1] = 2.0;
  d[2] = 3.0;
  d[3] = 4.0;
  d[4] = 5.0;
  d[5] = -1.0;

  QgsRasterMatrix m( 2, 3, d, -1.0 );

  QgsRasterCalcNode node( QgsRasterCalcNode::opPLUS, new QgsRasterCalcNode( &m ), new QgsRasterCalcNode( 5.0 ) );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );
  QMap<QString, QgsRasterBlock *> rasterData;

  QVERIFY( node.calculate( rasterData, result ) );

  QCOMPARE( result.data()[0], 6.0 );
  QCOMPARE( result.data()[1], 7.0 );
  QCOMPARE( result.data()[2], 8.0 );
  QCOMPARE( result.data()[3], 9.0 );
  QCOMPARE( result.data()[4], 10.0 );
  QCOMPARE( result.data()[5], -9999.0 );

  //also check adding no data number
  QgsRasterCalcNode nodeNoData( QgsRasterCalcNode::opPLUS, new QgsRasterCalcNode( &m ), new QgsRasterCalcNode( -9999 ) );
  QVERIFY( nodeNoData.calculate( rasterData, result ) );
  QCOMPARE( result.data()[0], -9999.0 );
}

void TestQgsRasterCalculator::dualOpMatrixMatrix()
{
  // test dual op run on matrix and matrix
  double *d = new double[6];
  d[0] = 1.0;
  d[1] = 2.0;
  d[2] = -2.0;
  d[3] = -1.0; //nodata
  d[4] = 5.0;
  d[5] = -1.0; //nodata
  QgsRasterMatrix m1( 2, 3, d, -1.0 );

  double *d2 = new double[6];
  d2[0] = -1.0;
  d2[1] = -2.0; //nodata
  d2[2] = 13.0;
  d2[3] = -2.0; //nodata
  d2[4] = 15.0;
  d2[5] = -1.0;
  QgsRasterMatrix m2( 2, 3, d2, -2.0 ); //different no data value

  QgsRasterCalcNode node( QgsRasterCalcNode::opPLUS, new QgsRasterCalcNode( &m1 ), new QgsRasterCalcNode( &m2 ) );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );
  QMap<QString, QgsRasterBlock *> rasterData;

  QVERIFY( node.calculate( rasterData, result ) );

  QCOMPARE( result.data()[0], 0.0 );
  QCOMPARE( result.data()[1], -9999.0 );
  QCOMPARE( result.data()[2], 11.0 );
  QCOMPARE( result.data()[3], -9999.0 );
  QCOMPARE( result.data()[4], 20.0 );
  QCOMPARE( result.data()[5], -9999.0 );
}

void TestQgsRasterCalculator::rasterRefOp()
{
  // test single op run on raster ref
  QgsRasterCalcNode node( QgsRasterCalcNode::opSIGN, new QgsRasterCalcNode( QStringLiteral( "raster" ) ), nullptr );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );
  QMap<QString, QgsRasterBlock *> rasterData;

  //first test invalid raster ref
  QVERIFY( !node.calculate( rasterData, result ) );

  //now create raster ref
  QgsRasterBlock m( Qgis::Float32, 2, 3 );
  m.setNoDataValue( -1.0 );
  m.setValue( 0, 0, 1.0 );
  m.setValue( 0, 1, 2.0 );
  m.setValue( 1, 0, 3.0 );
  m.setValue( 1, 1, 4.0 );
  m.setValue( 2, 0, 5.0 );
  m.setValue( 2, 1, -1.0 );
  rasterData.insert( QStringLiteral( "raster" ), &m );

  QVERIFY( node.calculate( rasterData, result ) );
  QCOMPARE( result.data()[0], -1.0 );
  QCOMPARE( result.data()[1], -2.0 );
  QCOMPARE( result.data()[2], -3.0 );
  QCOMPARE( result.data()[3], -4.0 );
  QCOMPARE( result.data()[4], -5.0 );
  QCOMPARE( result.data()[5], -9999.0 );
}

void TestQgsRasterCalculator::dualOpRasterRaster()
{
  // test dual op run on matrix and matrix

  QgsRasterBlock m1( Qgis::Float32, 2, 3 );
  m1.setNoDataValue( -1.0 );
  m1.setValue( 0, 0, 1.0 );
  m1.setValue( 0, 1, 2.0 );
  m1.setValue( 1, 0, -2.0 );
  m1.setValue( 1, 1, -1.0 ); //nodata
  m1.setValue( 2, 0, 5.0 );
  m1.setValue( 2, 1, -1.0 ); //nodata
  QMap<QString, QgsRasterBlock *> rasterData;
  rasterData.insert( QStringLiteral( "raster1" ), &m1 );

  QgsRasterBlock m2( Qgis::Float32, 2, 3 );
  m2.setNoDataValue( -2.0 ); //different no data value
  m2.setValue( 0, 0, -1.0 );
  m2.setValue( 0, 1, -2.0 ); //nodata
  m2.setValue( 1, 0, 13.0 );
  m2.setValue( 1, 1, -2.0 ); //nodata
  m2.setValue( 2, 0, 15.0 );
  m2.setValue( 2, 1, -1.0 );
  rasterData.insert( QStringLiteral( "raster2" ), &m2 );

  QgsRasterCalcNode node( QgsRasterCalcNode::opPLUS, new QgsRasterCalcNode( QStringLiteral( "raster1" ) ), new QgsRasterCalcNode( QStringLiteral( "raster2" ) ) );

  QgsRasterMatrix result;
  result.setNodataValue( -9999 );

  QVERIFY( node.calculate( rasterData, result ) );
  QCOMPARE( result.data()[0], 0.0 );
  QCOMPARE( result.data()[1], -9999.0 );
  QCOMPARE( result.data()[2], 11.0 );
  QCOMPARE( result.data()[3], -9999.0 );
  QCOMPARE( result.data()[4], 20.0 );
  QCOMPARE( result.data()[5], -9999.0 );
}

void TestQgsRasterCalculator::calcWithLayers()
{
  QgsRasterCalculatorEntry entry1;
  entry1.bandNumber = 1;
  entry1.raster = mpLandsatRasterLayer;
  entry1.ref = QStringLiteral( "landsat@1" );

  QgsRasterCalculatorEntry entry2;
  entry2.bandNumber = 2;
  entry2.raster = mpLandsatRasterLayer;
  entry2.ref = QStringLiteral( "landsat@2" );

  QVector<QgsRasterCalculatorEntry> entries;
  entries << entry1 << entry2;

  QgsCoordinateReferenceSystem crs;
  crs.createFromId( 32633, QgsCoordinateReferenceSystem::EpsgCrsId );
  QgsRectangle extent( 783235, 3348110, 783350, 3347960 );

  QTemporaryFile tmpFile;
  tmpFile.open(); // fileName is no avialable until open
  QString tmpName = tmpFile.fileName();
  tmpFile.close();

  QgsRasterCalculator rc( QStringLiteral( "\"landsat@1\" + 2" ),
                          tmpName,
                          QStringLiteral( "GTiff" ),
                          extent, crs, 2, 3, entries,
                          QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation() ), 0 );

  //open output file and check results
  QgsRasterLayer *result = new QgsRasterLayer( tmpName, QStringLiteral( "result" ) );
  QCOMPARE( result->width(), 2 );
  QCOMPARE( result->height(), 3 );
  QgsRasterBlock *block = result->dataProvider()->block( 1, extent, 2, 3 );
  QCOMPARE( block->value( 0, 0 ), 127.0 );
  QCOMPARE( block->value( 0, 1 ), 127.0 );
  QCOMPARE( block->value( 1, 0 ), 126.0 );
  QCOMPARE( block->value( 1, 1 ), 127.0 );
  QCOMPARE( block->value( 2, 0 ), 127.0 );
  QCOMPARE( block->value( 2, 1 ), 126.0 );
  delete result;
  delete block;

  //now try with 2 raster bands
  QgsRasterCalculator rc2( QStringLiteral( "\"landsat@1\" + \"landsat@2\"" ),
                           tmpName,
                           QStringLiteral( "GTiff" ),
                           extent, crs, 2, 3, entries,
                           QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc2.processCalculation() ), 0 );

  //open output file and check results
  result = new QgsRasterLayer( tmpName, QStringLiteral( "result" ) );
  QCOMPARE( result->width(), 2 );
  QCOMPARE( result->height(), 3 );
  block = result->dataProvider()->block( 1, extent, 2, 3 );
  QCOMPARE( block->value( 0, 0 ), 265.0 );
  QCOMPARE( block->value( 0, 1 ), 263.0 );
  QCOMPARE( block->value( 1, 0 ), 263.0 );
  QCOMPARE( block->value( 1, 1 ), 264.0 );
  QCOMPARE( block->value( 2, 0 ), 266.0 );
  QCOMPARE( block->value( 2, 1 ), 261.0 );
  delete result;
  delete block;
}

void TestQgsRasterCalculator::calcWithReprojectedLayers()
{
  QgsRasterCalculatorEntry entry1;
  entry1.bandNumber = 1;
  entry1.raster = mpLandsatRasterLayer;
  entry1.ref = QStringLiteral( "landsat@1" );

  QgsRasterCalculatorEntry entry2;
  entry2.bandNumber = 2;
  entry2.raster = mpLandsatRasterLayer4326;
  entry2.ref = QStringLiteral( "landsat_4326@2" );

  QVector<QgsRasterCalculatorEntry> entries;
  entries << entry1 << entry2;

  QgsCoordinateReferenceSystem crs;
  crs.createFromId( 32633, QgsCoordinateReferenceSystem::EpsgCrsId );
  QgsRectangle extent( 783235, 3348110, 783350, 3347960 );

  QTemporaryFile tmpFile;
  tmpFile.open(); // fileName is not available until open
  QString tmpName = tmpFile.fileName();
  tmpFile.close();

  QgsRasterCalculator rc( QStringLiteral( "\"landsat@1\" + \"landsat_4326@2\"" ),
                          tmpName,
                          QStringLiteral( "GTiff" ),
                          extent, crs, 2, 3, entries,
                          QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation() ), 0 );

  //open output file and check results
  QgsRasterLayer *result = new QgsRasterLayer( tmpName, QStringLiteral( "result" ) );
  QCOMPARE( result->width(), 2 );
  QCOMPARE( result->height(), 3 );
  QgsRasterBlock *block = result->dataProvider()->block( 1, extent, 2, 3 );
  QCOMPARE( block->value( 0, 0 ), 264.0 );
  QCOMPARE( block->value( 0, 1 ), 263.0 );
  QCOMPARE( block->value( 1, 0 ), 264.0 );
  QCOMPARE( block->value( 1, 1 ), 264.0 );
  QCOMPARE( block->value( 2, 0 ), 266.0 );
  QCOMPARE( block->value( 2, 1 ), 261.0 );
  delete result;
  delete block;
}

void TestQgsRasterCalculator::findNodes()
{

  std::unique_ptr< QgsRasterCalcNode > calcNode;

  auto _test =
    [ & ]( QString exp, const QgsRasterCalcNode::Type type ) -> QList<const QgsRasterCalcNode *>
  {
    QString error;
    calcNode.reset( QgsRasterCalcNode::parseRasterCalcString( exp, error ) );
    return calcNode->findNodes( type );
  };

  QCOMPARE( _test( QStringLiteral( "atan(\"raster@1\") * cos( 3  +  2 )" ), QgsRasterCalcNode::Type::tOperator ).length(), 4 );
  QCOMPARE( _test( QStringLiteral( "\"raster@1\"" ), QgsRasterCalcNode::Type::tOperator ).length(), 0 );
  QCOMPARE( _test( QStringLiteral( "\"raster@1\"" ), QgsRasterCalcNode::Type::tRasterRef ).length(), 1 );
  QCOMPARE( _test( QStringLiteral( "\"raster@1\"" ), QgsRasterCalcNode::Type::tMatrix ).length(), 0 );
  QCOMPARE( _test( QStringLiteral( "2 + 3" ), QgsRasterCalcNode::Type::tNumber ).length(), 2 );
  QCOMPARE( _test( QStringLiteral( "2 + 3" ), QgsRasterCalcNode::Type::tOperator ).length(), 1 );

}

void TestQgsRasterCalculator::testRasterEntries()
{
  // Create some test layers
  QList<QgsMapLayer *> layers;
  QgsRasterLayer *rlayer = new QgsRasterLayer( QStringLiteral( TEST_DATA_DIR ) + "/analysis/dem.tif",  QStringLiteral( "dem" ) );
  layers << rlayer;
  // Duplicate name, same source
  rlayer = new QgsRasterLayer( QStringLiteral( TEST_DATA_DIR ) + "/analysis/dem.tif",  QStringLiteral( "dem" ) );
  layers << rlayer;
  // Duplicated name different source
  rlayer = new QgsRasterLayer( QStringLiteral( TEST_DATA_DIR ) + "/analysis/dem_int16.tif",  QStringLiteral( "dem" ) );
  layers << rlayer;
  // Different name and different source
  rlayer = new QgsRasterLayer( QStringLiteral( TEST_DATA_DIR ) + "/analysis/slope.tif",  QStringLiteral( "slope" ) );
  layers << rlayer ;
  // Different name and same source
  rlayer = new QgsRasterLayer( QStringLiteral( TEST_DATA_DIR ) + "/analysis/slope.tif",  QStringLiteral( "slope2" ) );
  layers << rlayer ;
  QgsProject::instance()->addMapLayers( layers );
  QVector<QgsRasterCalculatorEntry> availableRasterBands = QgsRasterCalculatorEntry::rasterEntries();
  QMap<QString, QgsRasterCalculatorEntry> entryMap;
  for ( const auto &rb : qgis::as_const( availableRasterBands ) )
  {
    entryMap[rb.ref] = rb;
  }
  QStringList keys( entryMap.keys() );
  keys.sort();
  QCOMPARE( keys.join( ',' ), QStringLiteral( "dem@1,dem_1@1,landsat@1,landsat@2,landsat@3,landsat@4,"
            "landsat@5,landsat@6,landsat@7,landsat@8,landsat@9,"
            "landsat_4326@1,landsat_4326@2,landsat_4326@3,landsat_4326@4,"
            "landsat_4326@5,landsat_4326@6,landsat_4326@7,landsat_4326@8,landsat_4326@9,slope2@1" ) );
}

void TestQgsRasterCalculator::errors( )
{
  QgsRasterCalculatorEntry entry1;
  entry1.bandNumber = 0; // bad band
  entry1.raster = mpLandsatRasterLayer;
  entry1.ref = QStringLiteral( "landsat@0" );

  QVector<QgsRasterCalculatorEntry> entries;
  entries << entry1;

  QgsCoordinateReferenceSystem crs;
  crs.createFromId( 32633, QgsCoordinateReferenceSystem::EpsgCrsId );
  QgsRectangle extent( 783235, 3348110, 783350, 3347960 );

  QTemporaryFile tmpFile;
  tmpFile.open(); // fileName is not available until open
  QString tmpName = tmpFile.fileName();
  tmpFile.close();

  QgsRasterCalculator rc( QStringLiteral( "\"landsat@0\"" ),
                          tmpName,
                          QStringLiteral( "GTiff" ),
                          extent, crs, 2, 3, entries,
                          QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation() ), 6 );
  QCOMPARE( rc.lastError(), QStringLiteral( "Band number 0 is not valid for entry landsat@0" ) );

  entry1.bandNumber = 10; // bad band
  entries.clear();
  entries << entry1;
  rc = QgsRasterCalculator( QStringLiteral( "\"landsat@0\"" ),
                            tmpName,
                            QStringLiteral( "GTiff" ),
                            extent, crs, 2, 3, entries,
                            QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation() ), 6 );
  QCOMPARE( rc.lastError(), QStringLiteral( "Band number 10 is not valid for entry landsat@0" ) );


  // no raster
  entry1.raster = nullptr;
  entry1.bandNumber = 1;
  entries.clear();
  entries << entry1;
  rc = QgsRasterCalculator( QStringLiteral( "\"landsat@0\"" ),
                            tmpName,
                            QStringLiteral( "GTiff" ),
                            extent, crs, 2, 3, entries,
                            QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation() ), 2 );
  QCOMPARE( rc.lastError(), QStringLiteral( "No raster layer for entry landsat@0" ) );

  // bad driver
  entry1.raster = mpLandsatRasterLayer;
  entry1.bandNumber = 1;
  entries.clear();
  entries << entry1;
  rc = QgsRasterCalculator( QStringLiteral( "\"landsat@0\"" ),
                            tmpName,
                            QStringLiteral( "xxxxx" ),
                            extent, crs, 2, 3, entries,
                            QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation() ), 1 );
  QCOMPARE( rc.lastError(), QStringLiteral( "Could not obtain driver for xxxxx" ) );

  // bad filename
  rc = QgsRasterCalculator( QStringLiteral( "\"landsat@0\"" ),
                            QStringLiteral( "/goodluckwritinghere/blah/blah.tif" ),
                            QStringLiteral( "GTiff" ),
                            extent, crs, 2, 3, entries,
                            QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation() ), 1 );
  QCOMPARE( rc.lastError(), QStringLiteral( "Could not create output /goodluckwritinghere/blah/blah.tif" ) );

  // canceled
  QgsFeedback feedback;
  feedback.cancel();
  rc = QgsRasterCalculator( QStringLiteral( "\"landsat@0\"" ),
                            tmpName,
                            QStringLiteral( "GTiff" ),
                            extent, crs, 2, 3, entries,
                            QgsProject::instance()->transformContext() );
  QCOMPARE( static_cast< int >( rc.processCalculation( &feedback ) ), 3 );
  QVERIFY( rc.lastError().isEmpty() );
}

void TestQgsRasterCalculator::toString()
{
  auto _test = [ ]( QString exp, bool cStyle ) -> QString
  {
    QString error;
    std::unique_ptr< QgsRasterCalcNode > calcNode( QgsRasterCalcNode::parseRasterCalcString( exp, error ) );
    if ( ! error.isEmpty() )
      return error;
    return calcNode->toString( cStyle );
  };
  QCOMPARE( _test( QStringLiteral( "\"raster@1\"  + 2" ), false ), QString( "( \"raster@1\" + 2 )" ) );
  QCOMPARE( _test( QStringLiteral( "\"raster@1\"  +  2" ), true ), QString( "( \"raster@1\" + (float) ( 2 ) )" ) );
  QCOMPARE( _test( QStringLiteral( "\"raster@1\" ^ 3  +  2" ), false ), QString( "( \"raster@1\"^3 + 2 )" ) );
  QCOMPARE( _test( QStringLiteral( "\"raster@1\" ^ 3  +  2" ), true ), QString( "( pow( \"raster@1\", (float) ( 3 ) ) + (float) ( 2 ) )" ) );
  QCOMPARE( _test( QStringLiteral( "atan(\"raster@1\") * cos( 3  +  2 )" ), false ), QString( "atan( \"raster@1\" ) * cos( ( 3 + 2 ) )" ) );
  QCOMPARE( _test( QStringLiteral( "atan(\"raster@1\") * cos( 3  +  2 )" ), true ), QString( "atan( \"raster@1\" ) * cos( ( (float) ( 3 ) + (float) ( 2 ) ) )" ) );
  QCOMPARE( _test( QStringLiteral( "0.5 * ( 1.4 * (\"raster@1\" + 2) )" ), false ), QString( "0.5 * 1.4 * ( \"raster@1\" + 2 )" ) );
  QCOMPARE( _test( QStringLiteral( "0.5 * ( 1.4 * (\"raster@1\" + 2) )" ), true ), QString( "(float) ( 0.5 ) * (float) ( 1.4 ) * ( \"raster@1\" + (float) ( 2 ) )" ) );
}

void TestQgsRasterCalculator::calcFormulasWithReprojectedLayers()
{
  QgsRasterCalculatorEntry entry1;
  entry1.bandNumber = 1;
  entry1.raster = mpLandsatRasterLayer;
  entry1.ref = QStringLiteral( "landsat@1" );

  QgsRasterCalculatorEntry entry2;
  entry2.bandNumber = 2;
  entry2.raster = mpLandsatRasterLayer4326;
  entry2.ref = QStringLiteral( "landsat_4326@2" );

  QVector<QgsRasterCalculatorEntry> entries;
  entries << entry1 << entry2;

  QgsCoordinateReferenceSystem crs;
  crs.createFromId( 32633, QgsCoordinateReferenceSystem::EpsgCrsId );
  QgsRectangle extent( 783235, 3348110, 783350, 3347960 );


  auto _chk = [ = ]( const QString & formula, const std::vector<float> &values, bool useOpenCL )
  {

#ifdef HAVE_OPENCL
    if ( ! QgsOpenClUtils::available() )
      return ;
    QgsOpenClUtils::setEnabled( useOpenCL );
#else
    Q_UNUSED( useOpenCL )
#endif

    QTemporaryFile tmpFile;
    tmpFile.open(); // fileName is not available until open
    QString tmpName = tmpFile.fileName();
    tmpFile.close();
    QgsRasterCalculator rc( formula,
                            tmpName,
                            QStringLiteral( "GTiff" ),
                            extent, crs, 2, 3, entries,
                            QgsProject::instance()->transformContext() );
    QCOMPARE( static_cast< int >( rc.processCalculation() ), 0 );
    //open output file and check results
    QgsRasterLayer *result = new QgsRasterLayer( tmpName, QStringLiteral( "result" ) );
    QCOMPARE( result->width(), 2 );
    QCOMPARE( result->height(), 3 );
    QgsRasterBlock *block = result->dataProvider()->block( 1, extent, 2, 3 );
    qDebug() << block->value( 0, 0 ) << block->value( 0, 1 ) <<  block->value( 1, 0 ) <<  block->value( 1, 1 ) <<  block->value( 2, 0 ) <<  block->value( 2, 1 );
    const float epsilon { 0.0001f };
    QVERIFY( std::abs( ( static_cast<float>( block->value( 0, 0 ) ) - values[0] ) / values[0] ) < epsilon );
    QVERIFY( std::abs( ( static_cast<float>( block->value( 0, 1 ) ) - values[1] ) / values[1] ) < epsilon );
    QVERIFY( std::abs( ( static_cast<float>( block->value( 1, 0 ) ) - values[2] ) / values[2] ) < epsilon );
    QVERIFY( std::abs( ( static_cast<float>( block->value( 1, 1 ) ) - values[3] ) / values[3] ) < epsilon );
    QVERIFY( std::abs( ( static_cast<float>( block->value( 2, 0 ) ) - values[4] ) / values[4] ) < epsilon );
    QVERIFY( std::abs( ( static_cast<float>( block->value( 2, 1 ) ) - values[5] ) / values[5] ) < epsilon );
    delete result;
    delete block;
  };

  _chk( QStringLiteral( "\"landsat@1\" + \"landsat_4326@2\"" ), {264.0, 263.0, 264.0, 264.0, 266.0, 261.0}, false );
  _chk( QStringLiteral( "\"landsat@1\" + \"landsat_4326@2\"" ), {264.0, 263.0, 264.0, 264.0, 266.0, 261.0}, true );
  _chk( QStringLiteral( "\"landsat@1\"^2 + 3 + \"landsat_4326@2\"" ), {15767, 15766, 15519, 15767, 15769, 15516}, false );
  _chk( QStringLiteral( "\"landsat@1\"^2 + 3 + \"landsat_4326@2\"" ), {15767, 15766, 15519, 15767, 15769, 15516}, true );
  _chk( QStringLiteral( "0.5*((2*\"landsat@1\"+1)-sqrt((2*\"landsat@1\"+1)^2-8*(\"landsat@1\"-\"landsat_4326@2\")))" ), {-0.111504f, -0.103543f, -0.128448f, -0.111504f, -0.127425f, -0.104374f}, false );
  _chk( QStringLiteral( "0.5*((2*\"landsat@1\"+1)-sqrt((2*\"landsat@1\"+1)^2-8*(\"landsat@1\"-\"landsat_4326@2\")))" ), {-0.111504f, -0.103543f, -0.128448f, -0.111504f, -0.127425f, -0.104374f}, true );

}


QGSTEST_MAIN( TestQgsRasterCalculator )
#include "testqgsrastercalculator.moc"
