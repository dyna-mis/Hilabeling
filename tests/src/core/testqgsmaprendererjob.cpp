/***************************************************************************
     testqgsvectorfilewriter.cpp
     --------------------------------------
    Date                 : Sun Sep 16 12:22:54 AKDT 2007
    Copyright            : (C) 2007 by Tim Sutton
    Email                : tim @ linfiniti.com
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "qgstest.h"
#include <QObject>
#include <QString>
#include <QStringList>
#include <QPainter>
#include <QTime>
#include <QApplication>
#include <QDesktopServices>

//qgis includes...
#include <qgsvectorlayer.h> //defines QgsFieldMap
#include <qgsvectorfilewriter.h> //logic for writing shpfiles
#include <qgsfeature.h> //we will need to pass a bunch of these for each rec
#include <qgsgeometry.h> //each feature needs a geometry
#include <qgspoint.h> //we will use point geometry
#include <qgscoordinatereferencesystem.h> //needed for creating a srs
#include <qgsapplication.h> //search path for srs.db
#include <qgsfield.h>
#include <qgis.h> //defines GEOWkt
#include "qgsmaprenderersequentialjob.h"
#include <qgsmaplayer.h>
#include <qgsreadwritecontext.h>
#include <qgsvectorlayer.h>
#include <qgsapplication.h>
#include <qgsproviderregistry.h>
#include <qgsproject.h>

//qgs unit test utility class
#include "qgsrenderchecker.h"

/**
 * \ingroup UnitTests
 * This is a unit test for the QgsMapRendererJob class.
 * It will do some performance testing too
 *
 */
class TestQgsMapRendererJob : public QObject
{
    Q_OBJECT

  public:
    TestQgsMapRendererJob() = default;

    ~TestQgsMapRendererJob() override
    {
      delete mMapSettings;
    }

  private slots:
    void initTestCase();// will be called before the first testfunction is executed.
    void cleanupTestCase();// will be called after the last testfunction was executed.
    void init() {} // will be called before each testfunction is executed.
    void cleanup() {} // will be called after every testfunction.

    //! This method tests render performance
    void performanceTest();

    /**
     * This unit test checks if rendering of adjacent tiles (e.g. to render images for tile caches)
     * does not result in border effects
     */
    void testFourAdjacentTiles_data();
    void testFourAdjacentTiles();

  private:
    QString mEncoding;
    QgsVectorFileWriter::WriterError mError =  QgsVectorFileWriter::NoError ;
    QgsCoordinateReferenceSystem mCRS;
    QgsFields mFields;
    QgsMapSettings *mMapSettings = nullptr;
    QgsMapLayer *mpPolysLayer = nullptr;
    QString mReport;
};


void TestQgsMapRendererJob::initTestCase()
{
  //
  // Runs once before any tests are run
  //
  QgsApplication::init();
  QgsApplication::initQgis();
  QgsApplication::showSettings();

  mMapSettings = new QgsMapSettings();

  //create some objects that will be used in all tests...
  mEncoding = QStringLiteral( "UTF-8" );
  QgsField myField1( QStringLiteral( "Value" ), QVariant::Int, QStringLiteral( "int" ), 10, 0, QStringLiteral( "Value on lon" ) );
  mFields.append( myField1 );
  mCRS = QgsCoordinateReferenceSystem( GEOWKT );
  //
  // Create the test dataset if it doesn't exist
  //
  QString myDataDir( TEST_DATA_DIR ); //defined in CmakeLists.txt
  QString myTestDataDir = myDataDir + '/';
  QString myTmpDir = QDir::tempPath() + '/';
  QString myFileName = myTmpDir +  "maprender_testdata.gpkg";
  //copy over the default qml for our generated layer
  QString myQmlFileName = myTestDataDir +  "maprender_testdata.qml";
  QFile::remove( myTmpDir + "maprender_testdata.qml" );
  QVERIFY( QFile::copy( myQmlFileName, myTmpDir + "maprender_testdata.qml" ) );
  qDebug( "Checking test dataset exists...\n%s", myFileName.toLocal8Bit().constData() );
  if ( !QFile::exists( myFileName ) )
  {
    qDebug( "Creating test dataset: " );

    QgsVectorFileWriter myWriter( myFileName,
                                  mEncoding,
                                  mFields,
                                  QgsWkbTypes::Polygon,
                                  mCRS );
    double myInterval = 0.5;
    for ( double i = -180.0; i <= 180.0; i += myInterval )
    {
      for ( double j = -90.0; j <= 90.0; j += myInterval )
      {
        //
        // Create a polygon feature
        //
        QgsPolylineXY myPolyline;
        QgsPointXY myPoint1 = QgsPointXY( i, j );
        QgsPointXY myPoint2 = QgsPointXY( i + myInterval, j );
        QgsPointXY myPoint3 = QgsPointXY( i + myInterval, j + myInterval );
        QgsPointXY myPoint4 = QgsPointXY( i, j + myInterval );
        myPolyline << myPoint1 << myPoint2 << myPoint3 << myPoint4 << myPoint1;
        QgsPolygonXY myPolygon;
        myPolygon << myPolyline;
        //polygon: first item of the list is outer ring,
        // inner rings (if any) start from second item
        //
        // NOTE: don't delete this pointer again -
        // ownership is passed to the feature which will
        // delete it in its dtor!
        QgsGeometry mypPolygonGeometry = QgsGeometry::fromPolygonXY( myPolygon );
        QgsFeature myFeature;
        myFeature.setGeometry( mypPolygonGeometry );
        myFeature.initAttributes( 1 );
        myFeature.setAttribute( 0, i );
        //
        // Write the feature to the filewriter
        // and check for errors
        //
        QVERIFY( myWriter.addFeature( myFeature ) );
        mError = myWriter.hasError();
        if ( mError == QgsVectorFileWriter::ErrDriverNotFound )
        {
          std::cout << "Driver not found error" << std::endl;
        }
        else if ( mError == QgsVectorFileWriter::ErrCreateDataSource )
        {
          std::cout << "Create data source error" << std::endl;
        }
        else if ( mError == QgsVectorFileWriter::ErrCreateLayer )
        {
          std::cout << "Create layer error" << std::endl;
        }
        QVERIFY( mError == QgsVectorFileWriter::NoError );
      }
    }
  } //file exists
  //
  //create a poly layer that will be used in all tests...
  //
  QFileInfo myPolyFileInfo( myFileName );
  mpPolysLayer = new QgsVectorLayer( myPolyFileInfo.filePath(),
                                     myPolyFileInfo.completeBaseName(), QStringLiteral( "ogr" ) );
  QVERIFY( mpPolysLayer->isValid() );
  // Register the layer with the registry
  QgsProject::instance()->addMapLayers( QList<QgsMapLayer *>() << mpPolysLayer );
  // add the test layer to the maprender
  mMapSettings->setLayers( QList<QgsMapLayer *>() << mpPolysLayer );
  mReport += QLatin1String( "<h1>Map Render Tests</h1>\n" );
}

void TestQgsMapRendererJob::cleanupTestCase()
{
  QgsApplication::exitQgis();

  QString myReportFile = QDir::tempPath() + "/qgistest.html";
  QFile myFile( myReportFile );
  if ( myFile.open( QIODevice::WriteOnly | QIODevice::Append ) )
  {
    QTextStream myQTextStream( &myFile );
    myQTextStream << mReport;
    myFile.close();
    //QDesktopServices::openUrl( "file:///" + myReportFile );
  }
}

void TestQgsMapRendererJob::performanceTest()
{
  mMapSettings->setExtent( mpPolysLayer->extent() );
  QgsRenderChecker myChecker;
  myChecker.setControlName( QStringLiteral( "expected_maprender" ) );
  mMapSettings->setFlag( QgsMapSettings::Antialiasing );
  myChecker.setMapSettings( *mMapSettings );
  myChecker.setColorTolerance( 5 );
  bool myResultFlag = myChecker.runTest( QStringLiteral( "maprender" ) );
  mReport += myChecker.report();
  QVERIFY( myResultFlag );
}

void TestQgsMapRendererJob::testFourAdjacentTiles_data()
{
  QTest::addColumn<QStringList>( "bboxList" );
  QTest::addColumn<QString>( "controlName" );
  QTest::addColumn<QString>( "shapeFile" );
  QTest::addColumn<QString>( "qmlFile" );

  QString shapeFile = TEST_DATA_DIR + QStringLiteral( "/france_parts.shp" );
  QString qmlFile = TEST_DATA_DIR + QStringLiteral( "/adjacent_tiles/line_pattern_30_degree.qml" );
  QString controlName = QStringLiteral( "expected_adjacent_line_fill" );

  QStringList bboxList1;
  bboxList1 << QStringLiteral( "-1.5,48,-0.5,49" );
  bboxList1 << QStringLiteral( "-0.5,48,0.5,49" );
  bboxList1 << QStringLiteral( "-1.5,47,-0.5,48" );
  bboxList1 << QStringLiteral( "-0.5,47,0.5,48" );

  QTest::newRow( "adjacent_line_fill" ) << bboxList1 << controlName << shapeFile << qmlFile;

  qmlFile = TEST_DATA_DIR + QStringLiteral( "/adjacent_tiles/point_pattern_simple_marker.qml" );
  controlName = QStringLiteral( "expected_adjacent_marker_fill" );

  QTest::newRow( "adjacent_marker_fill" ) << bboxList1 << controlName << shapeFile << qmlFile;

  shapeFile = TEST_DATA_DIR + QStringLiteral( "/lines.shp" );
  qmlFile = TEST_DATA_DIR + QStringLiteral( "/adjacent_tiles/simple_line_dashed.qml" );
  controlName = QStringLiteral( "expected_adjacent_dashed_line" );

  QStringList bboxList2;
  bboxList2 << QStringLiteral( "-105,35,-95,45" );
  bboxList2 << QStringLiteral( "-95,35,-85,45" );
  bboxList2 << QStringLiteral( "-105,25,-95,35" );
  bboxList2 << QStringLiteral( "-95,25,-85,35" );

  QTest::newRow( "adjacent_dashed_line" ) << bboxList2 << controlName << shapeFile << qmlFile;
}

void TestQgsMapRendererJob::testFourAdjacentTiles()
{
  QFETCH( QStringList, bboxList );
  QFETCH( QString, controlName );
  QFETCH( QString, shapeFile );
  QFETCH( QString, qmlFile );

  QVERIFY( bboxList.size() == 4 );

  //create maplayer, set QML and add to maplayer registry
  QgsVectorLayer *vectorLayer = new QgsVectorLayer( shapeFile, QStringLiteral( "testshape" ), QStringLiteral( "ogr" ) );

  //todo: read QML
  QFile symbologyFile( qmlFile );
  if ( !symbologyFile.open( QIODevice::ReadOnly ) )
  {
    QFAIL( "Open symbology file failed" );
  }

  QDomDocument qmlDoc;
  if ( !qmlDoc.setContent( &symbologyFile ) )
  {
    QFAIL( "QML file not valid" );
  }

  QString errorMsg;
  QgsReadWriteContext context = QgsReadWriteContext();
  if ( !vectorLayer->readSymbology( qmlDoc.documentElement(), errorMsg, context ) )
  {
    QFAIL( errorMsg.toLocal8Bit().data() );
  }

  QgsProject::instance()->addMapLayers( QList<QgsMapLayer *>() << vectorLayer );

  QImage globalImage( 512, 512, QImage::Format_ARGB32_Premultiplied );
  globalImage.fill( Qt::white );
  QPainter globalPainter( &globalImage );

  for ( int i = 0; i < 4; ++i )
  {
    QgsMapSettings mapSettings;

    //extent
    QStringList rectCoords = bboxList.at( i ).split( ',' );
    if ( rectCoords.size() != 4 )
    {
      QFAIL( "bbox string invalid" );
    }
    QgsRectangle rect( rectCoords[0].toDouble(), rectCoords[1].toDouble(), rectCoords[2].toDouble(), rectCoords[3].toDouble() );
    mapSettings.setExtent( rect );
    mapSettings.setOutputSize( QSize( 256, 256 ) );
    mapSettings.setLayers( QList<QgsMapLayer *>() << vectorLayer );
    mapSettings.setFlags( QgsMapSettings::RenderMapTile );
    mapSettings.setOutputDpi( 96 );

    QgsMapRendererSequentialJob renderJob( mapSettings );
    renderJob.start();
    renderJob.waitForFinished();
    QImage img = renderJob.renderedImage();
    int globalImageX = ( i % 2 ) * 256;
    int globalImageY = ( i < 2 ) ? 0 : 256;
    globalPainter.drawImage( globalImageX, globalImageY, img );
  }

  QgsProject::instance()->removeMapLayers( QStringList() << vectorLayer->id() );

  QString renderedImagePath = QDir::tempPath() + "/" + QTest::currentDataTag() + QStringLiteral( ".png" );
  globalImage.save( renderedImagePath );

  QgsRenderChecker checker;

  checker.setControlPathPrefix( QStringLiteral( "adjacent_tiles" ) );
  checker.setControlName( controlName );
  bool result = checker.compareImages( QTest::currentDataTag(), 100, renderedImagePath );
  mReport += checker.report();
  QVERIFY( result );
}


QGSTEST_MAIN( TestQgsMapRendererJob )
#include "testqgsmaprendererjob.moc"


