/***************************************************************************
     testqgssvgcache.cpp
     --------------------
    Date                 : October 2017
    Copyright            : (C) 2017 by Nyall Dawson
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
#include <QObject>
#include <QString>
#include <QStringList>
#include <QApplication>
#include <QFileInfo>
#include <QDir>
#include <QDesktopServices>
#include <QPicture>
#include <QPainter>
#include <QtConcurrent>
#include <QElapsedTimer>
#include "qgssvgcache.h"
#include "qgsmultirenderchecker.h"
#include "qgsapplication.h"

/**
 * \ingroup UnitTests
 * This is a unit test for QgsSvgCache.
 */
class TestQgsSvgCache : public QObject
{
    Q_OBJECT

  private:

    QString mReport;

    bool imageCheck( const QString &testName, QImage &image, int mismatchCount );

  private slots:
    void initTestCase();// will be called before the first testfunction is executed.
    void cleanupTestCase();// will be called after the last testfunction was executed.
    void init() {} // will be called before each testfunction is executed.
    void cleanup() {} // will be called after every testfunction.
    void fillCache();
    void threadSafePicture();
    void threadSafeImage();
    void changeImage(); //check that cache is updated if svg source file changes
    void base64();

};


void TestQgsSvgCache::initTestCase()
{
  QgsApplication::init();
  QgsApplication::initQgis();
  mReport += "<h1>QgsSvgCache Tests</h1>\n";
}

void TestQgsSvgCache::cleanupTestCase()
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

void TestQgsSvgCache::fillCache()
{
  QgsSvgCache cache;
  // flood cache to fill it
  QString svgPath = TEST_DATA_DIR + QStringLiteral( "/sample_svg.svg" );
  bool fitInCache = false;

  // we loop forever, continually increasing the size of the requested
  // svg render. The continually changing image size should quickly fill
  // the svg cache size, forcing use of non-cached images.
  // We break after hitting a certain threshold of non-cached images,
  // (after testing that the result is non-null, i.e. rendered on demand,
  // not from cache).
  int uncached = 0;
  for ( double size = 1000; uncached < 10; size += 100 )
  {
    QImage image = cache.svgAsImage( svgPath, size, QColor( 255, 0, 0 ), QColor( 0, 255, 0 ), 1, 1, fitInCache );
    QVERIFY( !image.isNull() );
    if ( !fitInCache )
      uncached++;
  }
}

struct RenderPictureWrapper
{
  QgsSvgCache &cache;
  QString svgPath;
  double size = 100;
  explicit RenderPictureWrapper( QgsSvgCache &cache, const QString &svgPath )
    : cache( cache )
    , svgPath( svgPath )
  {}
  void operator()( int )
  {
    QPicture pic = cache.svgAsPicture( svgPath, size, QColor( 255, 0, 0 ), QColor( 0, 255, 0 ), 1, 1, true );
    QSize imageSize = pic.boundingRect().size();
    QImage image( imageSize, QImage::Format_ARGB32_Premultiplied );
    image.fill( 0 ); // transparent background
    QPainter p( &image );
    p.drawPicture( 0, 0, pic );
  }
};

void TestQgsSvgCache::threadSafePicture()
{
  // QPicture playback is NOT thread safe with implicitly shared copies - this
  // unit test checks that concurrent drawing of svg as QPicture from QgsSvgCache
  // returns a detached copy which is safe to use across threads

  // refs:
  // https://issues.qgis.org/issues/17077
  // https://issues.qgis.org/issues/17089

  QgsSvgCache cache;
  QString svgPath = TEST_DATA_DIR + QStringLiteral( "/sample_svg.svg" );

  // smash picture rendering over multiple threads
  QVector< int > list;
  list.resize( 100 );
  QtConcurrent::blockingMap( list, RenderPictureWrapper( cache, svgPath ) );
}


struct RenderImageWrapper
{
  QgsSvgCache &cache;
  QString svgPath;
  double size = 100;
  explicit RenderImageWrapper( QgsSvgCache &cache, const QString &svgPath )
    : cache( cache )
    , svgPath( svgPath )
  {}
  void operator()( int )
  {
    bool fitsInCache = false;
    QImage cachedImage = cache.svgAsImage( svgPath, size, QColor( 255, 0, 0 ), QColor( 0, 255, 0 ), 1, 1, fitsInCache );
    QImage image( cachedImage.size(), QImage::Format_ARGB32_Premultiplied );
    image.fill( 0 ); // transparent background
    QPainter p( &image );
    p.drawImage( 0, 0, cachedImage );
  }
};

void TestQgsSvgCache::threadSafeImage()
{
  // This unit test checks that concurrent rendering of svg as QImage from QgsSvgCache
  // works without issues across threads

  QgsSvgCache cache;
  QString svgPath = TEST_DATA_DIR + QStringLiteral( "/sample_svg.svg" );

  // smash image rendering over multiple threads
  QVector< int > list;
  list.resize( 100 );
  QtConcurrent::blockingMap( list, RenderImageWrapper( cache, svgPath ) );
}

void TestQgsSvgCache::changeImage()
{
  bool inCache;
  QgsSvgCache cache;
  // no minimum time between checks
  cache.mFileModifiedCheckTimeout = 0;

  //copy an image to the temp folder
  QString tempImagePath = QDir::tempPath() + "/svg_cache.svg";

  QString originalImage = TEST_DATA_DIR + QStringLiteral( "/test_symbol_svg.svg" );
  if ( QFileInfo::exists( tempImagePath ) )
    QFile::remove( tempImagePath );
  QFile::copy( originalImage, tempImagePath );

  //render it through the cache
  QImage img = cache.svgAsImage( tempImagePath, 200, QColor( 0, 0, 0 ), QColor( 0, 0, 0 ), 1.0,
                                 1.0, inCache );
  QVERIFY( imageCheck( "svgcache_changed_before", img, 30 ) );

  // wait a second so that modified time is different
  QElapsedTimer t;
  t.start();
  while ( !t.hasExpired( 1000 ) )
  {}

  //replace the image in the temp folder
  QString newImage = TEST_DATA_DIR + QStringLiteral( "/test_symbol_svg2.svg" );
  QFile::remove( tempImagePath );
  QFile::copy( newImage, tempImagePath );

  //re-render it
  img = cache.svgAsImage( tempImagePath, 200, QColor( 0, 0, 0 ), QColor( 0, 0, 0 ), 1.0,
                          1.0, inCache );
  QVERIFY( imageCheck( "svgcache_changed_after", img, 30 ) );

  // repeat, with minimum time between checks
  QgsSvgCache cache2;
  QFile::remove( tempImagePath );
  QFile::copy( originalImage, tempImagePath );
  img = cache2.svgAsImage( tempImagePath, 200, QColor( 0, 0, 0 ), QColor( 0, 0, 0 ), 1.0,
                           1.0, inCache );
  QVERIFY( imageCheck( "svgcache_changed_before", img, 30 ) );

  // wait a second so that modified time is different
  t.restart();
  while ( !t.hasExpired( 1000 ) )
  {}

  //replace the image in the temp folder
  QFile::remove( tempImagePath );
  QFile::copy( newImage, tempImagePath );

  //re-render it - not enough time has elapsed between checks, so file modification time will NOT be rechecked and
  // existing cached image should be used
  img = cache2.svgAsImage( tempImagePath, 200, QColor( 0, 0, 0 ), QColor( 0, 0, 0 ), 1.0,
                           1.0, inCache );
  QVERIFY( imageCheck( "svgcache_changed_before", img, 30 ) );
}

void TestQgsSvgCache::base64()
{
  // test rendering svgs encoded in base 64
  QgsSvgCache cache;
  bool inCache = false;

  // invalid base64 strings
  QImage img = cache.svgAsImage( QStringLiteral( "base64:" ), 200, QColor( 0, 0, 0 ), QColor( 0, 0, 0 ), 1.0,
                                 1.0, inCache );
  QVERIFY( imageCheck( QStringLiteral( "null_image" ), img, 0 ) );

  img = cache.svgAsImage( QStringLiteral( "base64:zzzzzzzzzzzzzzzzzzzz" ), 200, QColor( 0, 0, 0 ), QColor( 0, 0, 0 ), 1.0,
                          1.0, inCache );
  QVERIFY( imageCheck( QStringLiteral( "null_image" ), img, 0 ) );

  //valid base 64
  img = cache.svgAsImage( QStringLiteral( "base64:PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gR2VuZXJhdG9yOiBBZG9iZSBJbGx1c3RyYXRvciAxNi4wLjAsIFNWRyBFeHBvcnQgUGx1Zy1JbiAuIFNWRyBWZXJzaW9uOiA2LjAwIEJ1aWxkIDApICAtLT4KCjxzdmcKICAgeG1sbnM6ZGM9Imh0dHA6Ly9wdXJsLm9yZy9kYy9lbGVtZW50cy8xLjEvIgogICB4bWxuczpjYz0iaHR0cDovL2NyZWF0aXZlY29tbW9ucy5vcmcvbnMjIgogICB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiCiAgIHhtbG5zOnN2Zz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciCiAgIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIKICAgeG1sbnM6c29kaXBvZGk9Imh0dHA6Ly9zb2RpcG9kaS5zb3VyY2Vmb3JnZS5uZXQvRFREL3NvZGlwb2RpLTAuZHRkIgogICB4bWxuczppbmtzY2FwZT0iaHR0cDovL3d3dy5pbmtzY2FwZS5vcmcvbmFtZXNwYWNlcy9pbmtzY2FwZSIKICAgdmVyc2lvbj0iMS4xIgogICBpZD0ic3ZnMiIKICAgaW5rc2NhcGU6b3V0cHV0X2V4dGVuc2lvbj0ib3JnLmlua3NjYXBlLm91dHB1dC5zdmcuaW5rc2NhcGUiCiAgIHNvZGlwb2RpOmRvY25hbWU9InNob3BwaW5nX2RpeS5zdmciCiAgIGlua3NjYXBlOnZlcnNpb249IjAuOTEgcjEzNzI1IgogICBzb2RpcG9kaTp2ZXJzaW9uPSIwLjMyIgogICB4PSIwcHgiCiAgIHk9IjBweCIKICAgd2lkdGg9IjMwNi4zMzQ3NSIKICAgaGVpZ2h0PSI0ODQuNzk5OTkiCiAgIHZpZXdCb3g9IjAgMCAzMDYuMzM0NzUgNDg0Ljc5OTk5IgogICBlbmFibGUtYmFja2dyb3VuZD0ibmV3IDAgMCA1ODAgNTgwIgogICB4bWw6c3BhY2U9InByZXNlcnZlIj48bWV0YWRhdGEKICAgICBpZD0ibWV0YWRhdGEyNiI+PHJkZjpSREY+PGNjOldvcmsKICAgICAgICAgcmRmOmFib3V0PSIiPjxkYzpmb3JtYXQ+aW1hZ2Uvc3ZnK3htbDwvZGM6Zm9ybWF0PjxkYzp0eXBlCiAgICAgICAgICAgcmRmOnJlc291cmNlPSJodHRwOi8vcHVybC5vcmcvZGMvZGNtaXR5cGUvU3RpbGxJbWFnZSIgLz48ZGM6dGl0bGU+PC9kYzp0aXRsZT48L2NjOldvcms+PC9yZGY6UkRGPjwvbWV0YWRhdGE+PHNvZGlwb2RpOm5hbWVkdmlldwogICAgIHNob3dncmlkPSJmYWxzZSIKICAgICBpbmtzY2FwZTpjeT0iLTExNS4wNDM3MiIKICAgICBpbmtzY2FwZTpjeD0iMTMwLjIyMTYxIgogICAgIGlua3NjYXBlOnpvb209IjAuNDYwODM4NTYiCiAgICAgcGFnZWNvbG9yPSIjZmZmZmZmIgogICAgIGJvcmRlcmNvbG9yPSIjNjY2NjY2IgogICAgIGd1aWRldG9sZXJhbmNlPSIxMC4wIgogICAgIG9iamVjdHRvbGVyYW5jZT0iMTAuMCIKICAgICBncmlkdG9sZXJhbmNlPSIxMC4wIgogICAgIGJvcmRlcm9wYWNpdHk9IjEuMCIKICAgICBpZD0iYmFzZSIKICAgICBpbmtzY2FwZTpjdXJyZW50LWxheWVyPSJzdmcyIgogICAgIGlua3NjYXBlOndpbmRvdy15PSIzNCIKICAgICBpbmtzY2FwZTp3aW5kb3cteD0iNzUiCiAgICAgaW5rc2NhcGU6cGFnZW9wYWNpdHk9IjAuMCIKICAgICBpbmtzY2FwZTpwYWdlc2hhZG93PSIyIgogICAgIGlua3NjYXBlOndpbmRvdy13aWR0aD0iMTAxNCIKICAgICBpbmtzY2FwZTp3aW5kb3ctaGVpZ2h0PSI3MTEiCiAgICAgZml0LW1hcmdpbi10b3A9IjAiCiAgICAgZml0LW1hcmdpbi1sZWZ0PSIwIgogICAgIGZpdC1tYXJnaW4tcmlnaHQ9IjAiCiAgICAgZml0LW1hcmdpbi1ib3R0b209IjAiCiAgICAgaW5rc2NhcGU6d2luZG93LW1heGltaXplZD0iMCIgLz48ZGVmcwogICAgIGlkPSJkZWZzNCIgLz48ZwogICAgIGlkPSJsYXllcjMiCiAgICAgdHJhbnNmb3JtPSJtYXRyaXgoNDguMTQ5NjksMCwwLDQ4LjE0OTY5LC02NzMuMTA0NTMsLTgwLjkwNTc1MikiCiAgICAgaW5rc2NhcGU6bGFiZWw9IkxheW91dCIKICAgICBkaXNwbGF5PSJub25lIgogICAgIHN0eWxlPSJkaXNwbGF5Om5vbmUiPjxyZWN0CiAgICAgICBpZD0icmVjdDQxMzQiCiAgICAgICB4PSIxIgogICAgICAgeT0iMSIKICAgICAgIGRpc3BsYXk9ImlubGluZSIKICAgICAgIHdpZHRoPSIxMCIKICAgICAgIGhlaWdodD0iMTAiCiAgICAgICBzdHlsZT0iZGlzcGxheTppbmxpbmU7ZmlsbDpub25lO3N0cm9rZTojNzU3NTc1O3N0cm9rZS13aWR0aDowLjEiIC8+PHJlY3QKICAgICAgIGlkPSJyZWN0NDEzNiIKICAgICAgIHg9IjIiCiAgICAgICB5PSIyIgogICAgICAgZGlzcGxheT0iaW5saW5lIgogICAgICAgd2lkdGg9IjgiCiAgICAgICBoZWlnaHQ9IjgiCiAgICAgICBzdHlsZT0iZGlzcGxheTppbmxpbmU7ZmlsbDpub25lO3N0cm9rZTojNzU3NTc1O3N0cm9rZS13aWR0aDowLjEiIC8+PC9nPjxwYXRoCiAgICAgZmlsbD0icGFyYW0oZmlsbCkiCiAgICAgc3Ryb2tlPSJwYXJhbShvdXRsaW5lKSIKICAgICBzdHJva2Utd2lkdGg9InBhcmFtKG91dGxpbmUtd2lkdGgpIgogICAgIGQ9Ik0gMzA2LjI5NTc0LDE5LjU1NCBDIDMwNi4yMTk3NCw4LjU4OSAyOTYuNjg3NzQsMCAyODQuNTk1NzQsMCBsIC0wLjE3OSwwLjAwMSBjIC0xMC4xOTEsMCAtMTguNzksNi40ODggLTIwLjg5MSwxNS41ODYgLTE0LjM2LC0wLjkxNSAtMjEuMTM2LC0zLjYwMiAtMjguMjk4LC02LjQ0MSAtOC44MzQsLTMuNTAzIC0xNy45NjksLTcuMTI1IC00MS40MjEsLTguMjgyIGwgLTc3Ljc5LC0wLjcyIEMgNjcuNzkyNzQyLDAuNjgyIDI4LjgwODc0MiwyOC42MjQgMC4xNDU3NDIwMSw4My4xOTMgYyAtMC4yNjksMC41MTIgLTAuMTU4LDEuMTQgMC4yNywxLjUyNyAwLjQyOCwwLjM4OSAxLjA2Mjk5OTk5LDAuNDM4IDEuNTQ1OTk5OTksMC4xMjMgbCAyLjU2MywtMS42NzggYyAwLjE4OCwtMC4xMjMgMC4zNCwtMC4yOTMgMC40NCwtMC40OTQgOS42MzgsLTE5LjMxNiA0NS43ODUsLTM2LjI2IDc3LjM1NSwtMzYuMjYgMjAuNzk0OTk4LDAgMzUuMjgzOTk4LDcuNDg4IDQwLjgxNDk5OCwyMS4wODggMS41MjYsNy44ODUgOS43OCwxNS4zNDUgMTguNzUzLDE3LjExMSBsIDAuMjMyLDEzNi4xNTkgYyAtOC4wNDcsMC44NTQgLTE2LjI2NSw5LjQ4NCAtMTYuNDI5LDE3LjY1NiBsIC00LjI2NCwyMjQuMjg2IGMgLTAuMTI5LDYuODQ0IDEuNjEyLDEyLjE2OSA1LjE3NSwxNS44MyA1Ljk3NSw2LjEzOCAxNS41NjgsNi4yMTQgMjEuMyw2LjI1OSAwLjAwNCwwIDAuMDA4LDAgMC4wMTIsMCBsIDQxLjA2MSwtMC4wNjIgYyA2LjU1NCwtMC4wMyAxNS43OTQsLTIuNzggMjEuNjYsLTguODU1IDQuMTAxLC00LjI0NiA2LjA4NCwtOS41MjcgNS44OTYsLTE1LjY5MyBsIC02LjE2NCwtMjIxLjgxNSBjIC0wLjIxNiwtOC4wMzUgLTguNDU2LC0xNi42MTYgLTE2LjQ1NCwtMTcuNTA3IGwgLTAuMTc5LC0xMzYuNjA1IGMgMTIuMDU3LC0wLjkzNCAyMC4zNDgsLTYuMDU2IDI2LjAyNywtMTYuMDUzIDQuMjM0LC03LjQ1IDEwLjk5OSwtMTAuNzczIDIxLjkzNSwtMTAuNzczIDQuNzE3LDAgOS43NTgsMC41OTggMTQuNjMzLDEuMTc2IDIuMTE4LDAuMjUxIDQuMjk3LDAuNTEgNi40MjIsMC43MTUgbCAtMC4wMzcsNS42NTggYyAwLjA3NSwxMC44MTMgOS44NjIsMTkuNjA5IDIxLjg1NCwxOS42MDkgMTIuMDE4LC0wLjAxOSAyMS43ODIsLTguODM1IDIxLjc2NywtMTkuNjUyIGwgLTAuMDM5LC00NS4zODkgeiIKICAgICBpZD0icGF0aDI0IgogICAgIGlua3NjYXBlOmNvbm5lY3Rvci1jdXJ2YXR1cmU9IjAiIC8+PC9zdmc+" ),
                          200, QColor( 0, 0, 0 ), QColor( 0, 0, 0 ), 1.0,
                          1.0, inCache );
  QVERIFY( imageCheck( QStringLiteral( "svgcache_base64" ), img, 30 ) );

}

bool TestQgsSvgCache::imageCheck( const QString &testName, QImage &image, int mismatchCount )
{
  //draw background
  QImage imageWithBackground( image.width(), image.height(), QImage::Format_RGB32 );
  QgsRenderChecker::drawBackground( &imageWithBackground );
  QPainter painter( &imageWithBackground );
  painter.drawImage( 0, 0, image );
  painter.end();

  mReport += "<h2>" + testName + "</h2>\n";
  QString tempDir = QDir::tempPath() + '/';
  QString fileName = tempDir + testName + ".png";
  imageWithBackground.save( fileName, "PNG" );
  QgsRenderChecker checker;
  checker.setControlName( "expected_" + testName );
  checker.setRenderedImage( fileName );
  checker.setColorTolerance( 2 );
  bool resultFlag = checker.compareImages( testName, mismatchCount );
  mReport += checker.report();
  return resultFlag;
}
QGSTEST_MAIN( TestQgsSvgCache )
#include "testqgssvgcache.moc"
