/***************************************************************************
               qgsgpsinformationwidget.cpp  -  description
                             -------------------
    begin                : Sat Jan 01 2010
    copyright            : (C) 2010 by Tim Sutton and Marco Hugentobler
    email                : tim@linfiniti.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "qgsgpsinformationwidget.h"

#include "info.h"

#include "qgisapp.h"
#include "qgsapplication.h"
#include "qgscoordinatetransform.h"
#include "qgsfeatureaction.h"
#include "qgsgeometry.h"
#include "qgsgpsconnectionregistry.h"
#include "qgsgpsdetector.h"
#include "qgslayertreeview.h"
#include "qgslogger.h"
#include "qgsmaptooladdfeature.h"
#include "qgsnmeaconnection.h"
#include "qgspointxy.h"
#include "qgsproject.h"
#include "qgsrubberband.h"
#include "qgsvectordataprovider.h"
#include "qgsvectorlayer.h"
#include "qgswkbptr.h"
#include "qgssettings.h"
#include "qgsstatusbar.h"
#include "gmath.h"
#include "qgsmapcanvas.h"

// QWT Charting widget

#include <qwt_global.h>
#include <qwt_legend.h>
#include <qwt_plot.h>
#include <qwt_plot_grid.h>

#ifdef WITH_QWTPOLAR
// QWT Polar plot add on
#include <qwt_symbol.h>
#include <qwt_polar_grid.h>
#include <qwt_polar_curve.h>
#include <qwt_scale_engine.h>
#endif

#include <QMessageBox>
#include <QFileInfo>
#include <QColorDialog>
#include <QFileDialog>
#include <QPixmap>
#include <QPen>


const int MAXACQUISITIONINTERVAL = 3000; // max gps information acquisition suspension interval (in seconds)
const int MAXDISTANCETHRESHOLD = 200; // max gps distance threshold (in meters)


QgsGpsInformationWidget::QgsGpsInformationWidget( QgsMapCanvas *thepCanvas, QWidget *parent, Qt::WindowFlags f )
  : QWidget( parent, f )
  , mpCanvas( thepCanvas )
{
  setupUi( this );
  connect( mConnectButton, &QPushButton::toggled, this, &QgsGpsInformationWidget::mConnectButton_toggled );
  connect( mBtnTrackColor, &QgsColorButton::colorChanged, this, &QgsGpsInformationWidget::trackColorChanged );
  connect( mSpinTrackWidth, static_cast < void ( QSpinBox::* )( int ) > ( &QSpinBox::valueChanged ), this, &QgsGpsInformationWidget::mSpinTrackWidth_valueChanged );
  connect( mBtnPosition, &QToolButton::clicked, this, &QgsGpsInformationWidget::mBtnPosition_clicked );
  connect( mBtnSignal, &QToolButton::clicked, this, &QgsGpsInformationWidget::mBtnSignal_clicked );
  connect( mBtnSatellites, &QToolButton::clicked, this, &QgsGpsInformationWidget::mBtnSatellites_clicked );
  connect( mBtnOptions, &QToolButton::clicked, this, &QgsGpsInformationWidget::mBtnOptions_clicked );
  connect( mBtnDebug, &QToolButton::clicked, this, &QgsGpsInformationWidget::mBtnDebug_clicked );
  connect( mBtnRefreshDevices, &QToolButton::clicked, this, &QgsGpsInformationWidget::mBtnRefreshDevices_clicked );
  connect( mBtnAddVertex, &QPushButton::clicked, this, &QgsGpsInformationWidget::mBtnAddVertex_clicked );
  connect( mBtnCloseFeature, &QPushButton::clicked, this, &QgsGpsInformationWidget::mBtnCloseFeature_clicked );
  connect( mBtnResetFeature, &QToolButton::clicked, this, &QgsGpsInformationWidget::mBtnResetFeature_clicked );
  connect( mBtnLogFile, &QPushButton::clicked, this, &QgsGpsInformationWidget::mBtnLogFile_clicked );

  mpLastLayer = nullptr;

  mLastGpsPosition = QgsPointXY( 0.0, 0.0 );
  mLastNmeaPosition.lat = nmea_degree2radian( 0.0 );
  mLastNmeaPosition.lon = nmea_degree2radian( 0.0 );

  mpMapMarker = nullptr;
  mpRubberBand = nullptr;
  populateDevices();
  QWidget *mpHistogramWidget = mStackedWidget->widget( 1 );
#ifndef WITH_QWTPOLAR
  mBtnSatellites->setVisible( false );
#endif
  //
  // Set up the graph for signal strength
  //
  mpPlot = new QwtPlot( mpHistogramWidget );
  mpPlot->setAutoReplot( false );   // plot on demand
  //mpPlot->setTitle(QObject::tr("Signal Status"));
  //mpPlot->insertLegend(new QwtLegend(), QwtPlot::BottomLegend);
  // Set axis titles
  //mpPlot->setAxisTitle(QwtPlot::xBottom, QObject::tr("Satellite"));
  //mpPlot->setAxisTitle(QwtPlot::yLeft, QObject::tr("Value"));
  mpPlot->setAxisScale( QwtPlot::xBottom, 0, 20 );
  mpPlot->setAxisScale( QwtPlot::yLeft, 0, 100 );  // max is 50dB SNR, I believe - SLM
  // add a grid
  //QwtPlotGrid * mypGrid = new QwtPlotGrid();
  //mypGrid->attach( mpPlot );
  //display satellites first
  mpCurve = new QwtPlotCurve();
  mpCurve->setRenderHint( QwtPlotItem::RenderAntialiased );
  mpCurve->setPen( QPen( Qt::blue ) );
  mpCurve->setBrush( QBrush( Qt::blue ) );
  mpPlot->enableAxis( QwtPlot::yLeft, false );
  mpPlot->enableAxis( QwtPlot::xBottom, false );
  mpCurve->attach( mpPlot );
  //ensure all children get removed
  mpPlot->setAutoDelete( true );
  QVBoxLayout *mpHistogramLayout = new QVBoxLayout( mpHistogramWidget );
  mpHistogramLayout->setContentsMargins( 0, 0, 0, 0 );
  mpHistogramLayout->addWidget( mpPlot );
  mpHistogramWidget->setLayout( mpHistogramLayout );

  //
  // Set up the polar graph for satellite pos
  //
#ifdef WITH_QWTPOLAR
  QWidget *mpPolarWidget = mStackedWidget->widget( 2 );
  mpSatellitesWidget = new QwtPolarPlot( /*QwtText( tr( "Satellite View" ), QwtText::PlainText ),*/ mpPolarWidget );  // possible title for graph removed for now as it is too large in small windows
  mpSatellitesWidget->setAutoReplot( false );   // plot on demand (after all data has been handled)
  mpSatellitesWidget->setPlotBackground( Qt::white );
  // scales
  mpSatellitesWidget->setScale( QwtPolar::ScaleAzimuth,
                                360, //min - reverse the min/max values to get compass orientation - increasing clockwise
                                0, //max
                                90 //interval - just show cardinal and intermediate (NE, N, NW, etc.) compass points (in degrees)
                              );
  mpSatellitesWidget->setAzimuthOrigin( M_PI_2 );    // to get compass orientation - need to rotate 90 deg. ccw; this is in Radians (not indicated in QwtPolarPlot docs)

//  mpSatellitesWidget->setScaleMaxMinor( QwtPolar::ScaleRadius, 2 );  // seems unnecessary
  mpSatellitesWidget->setScale( QwtPolar::ScaleRadius,
                                90, //min - reverse the min/max to get 0 at edge, 90 at center
                                0, //max
                                45 //interval
                              );

  // grids, axes
  mpSatellitesGrid = new QwtPolarGrid();
  mpSatellitesGrid->setGridAttribute( QwtPolarGrid::AutoScaling, false );   // This fixes the issue of autoscaling on the Radius grid. It is ON by default AND is separate from the scaleData.doAutoScale in QwtPolarPlot::setScale(), etc. THIS IS VERY TRICKY!
  mpSatellitesGrid->setPen( QPen( Qt::black ) );
  QPen minorPen( Qt::gray );  // moved outside of for loop; NOTE setting the minor pen isn't necessary if the minor grids aren't shown
  for ( int scaleId = 0; scaleId < QwtPolar::ScaleCount; scaleId++ )
  {
    //mpSatellitesGrid->showGrid( scaleId );
    //mpSatellitesGrid->showMinorGrid(scaleId);
    mpSatellitesGrid->setMinorGridPen( scaleId, minorPen );
  }
//  mpSatellitesGrid->setAxisPen( QwtPolar::AxisAzimuth, QPen( Qt::black ) );

  mpSatellitesGrid->showAxis( QwtPolar::AxisAzimuth, true );
  mpSatellitesGrid->showAxis( QwtPolar::AxisLeft, false ); //alt axis
  mpSatellitesGrid->showAxis( QwtPolar::AxisRight, false );//alt axis
  mpSatellitesGrid->showAxis( QwtPolar::AxisTop, false );//alt axis
  mpSatellitesGrid->showAxis( QwtPolar::AxisBottom, false );//alt axis
  mpSatellitesGrid->showGrid( QwtPolar::ScaleAzimuth, false ); // hide the grid; just show ticks at edge
  mpSatellitesGrid->showGrid( QwtPolar::ScaleRadius, true );
//  mpSatellitesGrid->showMinorGrid( QwtPolar::ScaleAzimuth, true );
  mpSatellitesGrid->showMinorGrid( QwtPolar::ScaleRadius, true );   // for 22.5, 67.5 degree circles
  mpSatellitesGrid->attach( mpSatellitesWidget );

  //QwtLegend *legend = new QwtLegend;
  //mpSatellitesWidget->insertLegend(legend, QwtPolarPlot::BottomLegend);
  QVBoxLayout *mpPolarLayout = new QVBoxLayout( mpPolarWidget );
  mpPolarLayout->setContentsMargins( 0, 0, 0, 0 );
  mpPolarLayout->addWidget( mpSatellitesWidget );
  mpPolarWidget->setLayout( mpPolarLayout );

  // replot on command
  mpSatellitesWidget->replot();
#endif
  mpPlot->replot();

  mBtnTrackColor->setAllowOpacity( true );
  mBtnTrackColor->setColorDialogTitle( tr( "Track Color" ) );
  // Restore state
  QgsSettings mySettings;
  mGroupShowMarker->setChecked( mySettings.value( QStringLiteral( "gps/showMarker" ), "true" ).toBool() );
  mSliderMarkerSize->setValue( mySettings.value( QStringLiteral( "gps/markerSize" ), "12" ).toInt() );
  mSpinTrackWidth->setValue( mySettings.value( QStringLiteral( "gps/trackWidth" ), "2" ).toInt() );
  mBtnTrackColor->setColor( mySettings.value( QStringLiteral( "gps/trackColor" ), QColor( Qt::red ) ).value<QColor>() );
  QString myPortMode = mySettings.value( QStringLiteral( "gps/portMode" ), "scanPorts" ).toString();

  mSpinMapExtentMultiplier->setValue( mySettings.value( QStringLiteral( "gps/mapExtentMultiplier" ), "50" ).toInt() );
  mDateTimeFormat = mySettings.value( QStringLiteral( "gps/dateTimeFormat" ), "" ).toString(); // zero-length string signifies default format

  mGpsdHost->setText( mySettings.value( QStringLiteral( "gps/gpsdHost" ), "localhost" ).toString() );
  mGpsdPort->setText( mySettings.value( QStringLiteral( "gps/gpsdPort" ), 2947 ).toString() );
  mGpsdDevice->setText( mySettings.value( QStringLiteral( "gps/gpsdDevice" ) ).toString() );

  //port mode
  if ( myPortMode == QLatin1String( "scanPorts" ) )
  {
    mRadAutodetect->setChecked( true );
  }
  else if ( myPortMode == QLatin1String( "internalGPS" ) )
  {
    mRadInternal->setChecked( true );
  }
  else if ( myPortMode == QLatin1String( "explicitPort" ) )
  {
    mRadUserPath->setChecked( true );
  }
  else if ( myPortMode == QLatin1String( "gpsd" ) )
  {
    mRadGpsd->setChecked( true );
  }
  //hide the internal port method if build is without QtLocation
#ifndef HAVE_QT_MOBILITY_LOCATION
  if ( mRadInternal->isChecked() )
  {
    mRadAutodetect->setChecked( true );
  }
  mRadInternal->hide();
#endif

  //auto digitizing behavior
  mCbxAutoAddVertices->setChecked( mySettings.value( QStringLiteral( "gps/autoAddVertices" ), "false" ).toBool() );

  mCbxAutoCommit->setChecked( mySettings.value( QStringLiteral( "gps/autoCommit" ), "false" ).toBool() );

  //pan mode
  QString myPanMode = mySettings.value( QStringLiteral( "gps/panMode" ), "recenterWhenNeeded" ).toString();
  if ( myPanMode == QLatin1String( "none" ) )
  {
    radNeverRecenter->setChecked( true );
  }
  else if ( myPanMode == QLatin1String( "recenterAlways" ) )
  {
    radRecenterMap->setChecked( true );
  }
  else
  {
    radRecenterWhenNeeded->setChecked( true );
  }

  mWgs84CRS = QgsCoordinateReferenceSystem::fromOgcWmsCrs( QStringLiteral( "EPSG:4326" ) );

  mBtnDebug->setVisible( mySettings.value( QStringLiteral( "gps/showDebug" ), "false" ).toBool() );  // use a registry setting to control - power users/devs could set it

  // status = unknown
  setStatusIndicator( NoData );

  //SLM - added functionality
  mLogFile = nullptr;

  connect( QgisApp::instance(), &QgisApp::activeLayerChanged,
           this, &QgsGpsInformationWidget::updateCloseFeatureButton );

  mStackedWidget->setCurrentIndex( 3 ); // force to Options
  mBtnPosition->setFocus( Qt::TabFocusReason );

  mAcquisitionIntValidator = new QIntValidator( 0, MAXACQUISITIONINTERVAL, this );
  mDistanceThresholdValidator = new QIntValidator( 0, MAXDISTANCETHRESHOLD, this );
  mAcquisitionTimer = std::unique_ptr<QTimer>( new QTimer( this ) );
  mAcquisitionTimer->setSingleShot( true );
  mCboAcquisitionInterval->addItem( QStringLiteral( "0" ), 0 );
  mCboAcquisitionInterval->addItem( QStringLiteral( "2" ), 2 );
  mCboAcquisitionInterval->addItem( QStringLiteral( "5" ), 5 );
  mCboAcquisitionInterval->addItem( QStringLiteral( "10" ), 10 );
  mCboAcquisitionInterval->addItem( QStringLiteral( "15" ), 15 );
  mCboAcquisitionInterval->addItem( QStringLiteral( "30" ), 30 );
  mCboAcquisitionInterval->addItem( QStringLiteral( "60" ), 60 );
  mCboDistanceThreshold->addItem( QStringLiteral( "0" ), 0 );
  mCboDistanceThreshold->addItem( QStringLiteral( "3" ), 3 );
  mCboDistanceThreshold->addItem( QStringLiteral( "5" ), 5 );
  mCboDistanceThreshold->addItem( QStringLiteral( "10" ), 10 );
  mCboDistanceThreshold->addItem( QStringLiteral( "15" ), 15 );

  mCboAcquisitionInterval->setValidator( mAcquisitionIntValidator );
  mCboDistanceThreshold->setValidator( mDistanceThresholdValidator );
  mCboAcquisitionInterval->setCurrentText( mySettings.value( QStringLiteral( "gps/acquisitionInterval" ), 0 ).toString() );
  mCboDistanceThreshold->setCurrentText( mySettings.value( QStringLiteral( "gps/distanceThreshold" ), 0 ).toString() );

  connect( mAcquisitionTimer.get(), &QTimer::timeout,
           this, &QgsGpsInformationWidget::switchAcquisition );
  connect( mCboAcquisitionInterval, qgis::overload< const QString & >::of( &QComboBox::currentTextChanged ),
           this, &QgsGpsInformationWidget::cboAcquisitionIntervalEdited );
  connect( mCboDistanceThreshold, qgis::overload< const QString & >::of( &QComboBox::currentTextChanged ),
           this, &QgsGpsInformationWidget::cboDistanceThresholdEdited );
}

QgsGpsInformationWidget::~QgsGpsInformationWidget()
{
  if ( mNmea )
  {
    disconnectGps();
  }

  delete mpMapMarker;
  delete mpRubberBand;

#ifdef WITH_QWTPOLAR
  delete mpSatellitesGrid;
#endif

  QgsSettings mySettings;
  mySettings.setValue( QStringLiteral( "gps/lastPort" ), mCboDevices->currentData().toString() );
  mySettings.setValue( QStringLiteral( "gps/trackWidth" ), mSpinTrackWidth->value() );
  mySettings.setValue( QStringLiteral( "gps/trackColor" ), mBtnTrackColor->color() );
  mySettings.setValue( QStringLiteral( "gps/markerSize" ), mSliderMarkerSize->value() );
  mySettings.setValue( QStringLiteral( "gps/showMarker" ), mGroupShowMarker->isChecked() );
  mySettings.setValue( QStringLiteral( "gps/autoAddVertices" ), mCbxAutoAddVertices->isChecked() );
  mySettings.setValue( QStringLiteral( "gps/autoCommit" ), mCbxAutoCommit->isChecked() );
  mySettings.setValue( QStringLiteral( "gps/acquisitionInterval" ), mCboAcquisitionInterval->currentText() );
  mySettings.setValue( QStringLiteral( "gps/distanceThreshold" ), mCboDistanceThreshold->currentText() );

  mySettings.setValue( QStringLiteral( "gps/mapExtentMultiplier" ), mSpinMapExtentMultiplier->value() );

  // scan, explicit port or gpsd
  if ( mRadAutodetect->isChecked() )
  {
    mySettings.setValue( QStringLiteral( "gps/portMode" ), "scanPorts" );
  }
  else if ( mRadInternal->isChecked() )
  {
    mySettings.setValue( QStringLiteral( "gps/portMode" ), "internalGPS" );
  }
  else if ( mRadUserPath->isChecked() )
  {
    mySettings.setValue( QStringLiteral( "gps/portMode" ), "explicitPort" );
  }
  else
  {
    mySettings.setValue( QStringLiteral( "gps/portMode" ), "gpsd" );
  }

  mySettings.setValue( QStringLiteral( "gps/gpsdHost" ), mGpsdHost->text() );
  mySettings.setValue( QStringLiteral( "gps/gpsdPort" ), mGpsdPort->text().toInt() );
  mySettings.setValue( QStringLiteral( "gps/gpsdDevice" ), mGpsdDevice->text() );

  // pan mode
  if ( radRecenterMap->isChecked() )
  {
    mySettings.setValue( QStringLiteral( "gps/panMode" ), "recenterAlways" );
  }
  else if ( radRecenterWhenNeeded->isChecked() )
  {
    mySettings.setValue( QStringLiteral( "gps/panMode" ), "recenterWhenNeeded" );
  }
  else
  {
    mySettings.setValue( QStringLiteral( "gps/panMode" ), "none" );
  }

}

void QgsGpsInformationWidget::mSpinTrackWidth_valueChanged( int value )
{
  if ( mpRubberBand )
  {
    mpRubberBand->setWidth( value );
    mpRubberBand->update();
  }
}

void QgsGpsInformationWidget::trackColorChanged( const QColor &color )
{
  if ( color.isValid() )  // check that a color was picked
  {
    if ( mpRubberBand )
    {
      mpRubberBand->setColor( color );
      mpRubberBand->update();
    }
  }
}

void QgsGpsInformationWidget::mBtnPosition_clicked()
{
  mStackedWidget->setCurrentIndex( 0 );
  if ( mNmea )
    displayGPSInformation( mNmea->currentGPSInformation() );
}

void QgsGpsInformationWidget::mBtnSignal_clicked()
{
  mStackedWidget->setCurrentIndex( 1 );
  if ( mNmea )
    displayGPSInformation( mNmea->currentGPSInformation() );
}

void QgsGpsInformationWidget::mBtnSatellites_clicked()
{
  mStackedWidget->setCurrentIndex( 2 );
  if ( mNmea )
    displayGPSInformation( mNmea->currentGPSInformation() );
}

void QgsGpsInformationWidget::mBtnOptions_clicked()
{
  mStackedWidget->setCurrentIndex( 3 );
}

void QgsGpsInformationWidget::mBtnDebug_clicked()
{
  mStackedWidget->setCurrentIndex( 4 );
}

void QgsGpsInformationWidget::mConnectButton_toggled( bool flag )
{
  if ( flag )
  {
    connectGps();
  }
  else
  {
    disconnectGps();
  }
}

void QgsGpsInformationWidget::connectGps()
{
  // clear position page fields to give better indication that something happened (or didn't happen)
  mTxtLatitude->clear();
  mTxtLongitude->clear();
  mTxtAltitude->clear();
  mTxtDateTime->clear();
  mTxtSpeed->clear();
  mTxtDirection->clear();
  mTxtHdop->clear();
  mTxtVdop->clear();
  mTxtPdop->clear();
  mTxtFixMode->clear();
  mTxtFixType->clear();
  mTxtQuality->clear();
  mTxtSatellitesUsed->clear();
  mTxtStatus->clear();

  mLastGpsPosition = QgsPointXY( 0.0, 0.0 );

  QString port;

  if ( mRadUserPath->isChecked() )
  {
    port = mCboDevices->currentData().toString();

    if ( port.isEmpty() )
    {
      QMessageBox::information( this, tr( "/gps" ), tr( "No path to the GPS port "
                                "is specified. Please enter a path then try again." ) );
      //toggle the button back off
      mConnectButton->setChecked( false );
      return;
    }
  }
  else if ( mRadGpsd->isChecked() )
  {
    port = QStringLiteral( "%1:%2:%3" ).arg( mGpsdHost->text(), mGpsdPort->text(), mGpsdDevice->text() );
  }
  else if ( mRadInternal->isChecked() )
  {
    port = QStringLiteral( "internalGPS" );
  }

  mGPSPlainTextEdit->appendPlainText( tr( "Connecting???" ) );
  showStatusBarMessage( tr( "Connecting to GPS device %1???" ).arg( port ) );

  QgsGpsDetector *detector = new QgsGpsDetector( port );
  connect( detector, static_cast < void ( QgsGpsDetector::* )( QgsGpsConnection * ) > ( &QgsGpsDetector::detected ), this, &QgsGpsInformationWidget::connected );
  connect( detector, &QgsGpsDetector::detectionFailed, this, &QgsGpsInformationWidget::timedout );
  detector->advance();   // start the detection process
}

void QgsGpsInformationWidget::timedout()
{
  mConnectButton->setChecked( false );
  mNmea = nullptr;
  mGPSPlainTextEdit->appendPlainText( tr( "Timed out!" ) );
  showStatusBarMessage( tr( "Failed to connect to GPS device." ) );
}

void QgsGpsInformationWidget::connected( QgsGpsConnection *conn )
{
  mNmea = conn;
  connect( mNmea, &QgsGpsConnection::stateChanged,
           this, &QgsGpsInformationWidget::displayGPSInformation );
  mGPSPlainTextEdit->appendPlainText( tr( "Connected!" ) );
  mConnectButton->setText( tr( "Dis&connect" ) );
  //insert connection into registry such that it can also be used by other dialogs or plugins
  QgsApplication::gpsConnectionRegistry()->registerConnection( mNmea );
  showStatusBarMessage( tr( "Connected to GPS device." ) );

  if ( mLogFileGroupBox->isChecked() && ! mTxtLogFile->text().isEmpty() )
  {
    if ( ! mLogFile )
    {
      mLogFile = new QFile( mTxtLogFile->text() );
    }

    if ( mLogFile->open( QIODevice::Append ) )  // open in binary and explicitly output CR + LF per NMEA
    {
      mLogFileTextStream.setDevice( mLogFile );

      // crude way to separate chunks - use when manually editing file - NMEA parsers should discard
      mLogFileTextStream << "====" << "\r\n";

      connect( mNmea, &QgsGpsConnection::nmeaSentenceReceived, this, &QgsGpsInformationWidget::logNmeaSentence ); // added to handle raw data
    }
    else  // error opening file
    {
      delete mLogFile;
      mLogFile = nullptr;

      // need to indicate why - this just reports that an error occurred
      showStatusBarMessage( tr( "Error opening log file." ) );
    }
  }
}

void QgsGpsInformationWidget::disconnectGps()
{
  if ( mLogFile && mLogFile->isOpen() )
  {
    disconnect( mNmea, &QgsGpsConnection::nmeaSentenceReceived, this, &QgsGpsInformationWidget::logNmeaSentence );
    mLogFile->close();
    delete mLogFile;
    mLogFile = nullptr;
  }

  QgsApplication::gpsConnectionRegistry()->unregisterConnection( mNmea );
  delete mNmea;
  mNmea = nullptr;
  if ( mpMapMarker )  // marker should not be shown on GPS disconnected - not current position
  {
    delete mpMapMarker;
    mpMapMarker = nullptr;
  }
  mGPSPlainTextEdit->appendPlainText( tr( "Disconnected???" ) );
  mConnectButton->setChecked( false );
  mConnectButton->setText( tr( "&Connect" ) );
  showStatusBarMessage( tr( "Disconnected from GPS device." ) );

  setStatusIndicator( NoData );
}

void QgsGpsInformationWidget::displayGPSInformation( const QgsGpsInformation &info )
{
  QVector<QPointF> data;

  // set validity flag and status from GPS data
  // based on GGA, GSA and RMC sentences - the logic does not require all
  bool validFlag = false; // true if GPS indicates position fix
  FixStatus fixStatus = NoData;

  // no fix if any of the three report bad; default values are invalid values and won't be changed if the corresponding NMEA msg is not received
  if ( info.status == 'V' || info.fixType == NMEA_FIX_BAD || info.quality == 0 ) // some sources say that 'V' indicates position fix, but is below acceptable quality
  {
    fixStatus = NoFix;
  }
  else if ( info.fixType == NMEA_FIX_2D ) // 2D indication (from GGA)
  {
    fixStatus = Fix2D;
    validFlag = true;
  }
  else if ( info.status == 'A' || info.fixType == NMEA_FIX_3D || info.quality > 0 ) // good
  {
    fixStatus = Fix3D;
    validFlag = true;
  }
  else  // unknown status (not likely)
  {

  }

  // set visual status indicator -- do only on change of state
  if ( fixStatus != mLastFixStatus )
  {
    setStatusIndicator( fixStatus );
  }

  if ( mStackedWidget->currentIndex() == 1 && info.satInfoComplete ) //signal
  {
    mpPlot->setAxisScale( QwtPlot::xBottom, 0, info.satellitesInView.size() );
  } //signal
#ifdef WITH_QWTPOLAR
  if ( mStackedWidget->currentIndex() == 2 && info.satInfoComplete ) //satellites
  {
    while ( !mMarkerList.isEmpty() )
    {
      delete mMarkerList.takeFirst();
    }
  } //satellites
#endif
  if ( mStackedWidget->currentIndex() == 4 ) //debug
  {
    mGPSPlainTextEdit->clear();
  } //debug

  for ( int i = 0; i < info.satellitesInView.size(); ++i ) //satellite processing loop
  {
    QgsSatelliteInfo currentInfo = info.satellitesInView.at( i );

    if ( mStackedWidget->currentIndex() == 1 && info.satInfoComplete ) //signal
    {
      data << QPointF( i, 0 );
      data << QPointF( i, currentInfo.signal );
      data << QPointF( i + 1, currentInfo.signal );
      data << QPointF( i + 1, 0 );
    } //signal

    if ( mStackedWidget->currentIndex() == 2 && info.satInfoComplete ) //satellites
    {
      QColor bg( Qt::white ); // moved several items outside of the following if block to minimize loop time
      bg.setAlpha( 200 );
      QColor myColor;

      // Add a marker to the polar plot
      if ( currentInfo.id > 0 )       // don't show satellite if id=0 (no satellite indication)
      {
#ifdef WITH_QWTPOLAR
        QwtPolarMarker *mypMarker = new QwtPolarMarker();
#if (QWT_POLAR_VERSION<0x010000)
        mypMarker->setPosition( QwtPolarPoint( currentInfo.azimuth, currentInfo.elevation ) );
#else
        mypMarker->setPosition( QwtPointPolar( currentInfo.azimuth, currentInfo.elevation ) );
#endif
#endif
        if ( currentInfo.signal < 30 ) //weak signal
        {
          myColor = Qt::red;
        }
        else
        {
          myColor = Qt::black; //strong signal
        }
#ifdef WITH_QWTPOLAR
        QBrush symbolBrush( Qt::black );
        QSize markerSize( 9, 9 );
        QBrush textBgBrush( bg );
#if (QWT_POLAR_VERSION<0x010000)
        mypMarker->setSymbol( QwtSymbol( QwtSymbol::Ellipse,
                                         symbolBrush, QPen( myColor ), markerSize ) );
#else
        mypMarker->setSymbol( new QwtSymbol( QwtSymbol::Ellipse,
                                             symbolBrush, QPen( myColor ), markerSize ) );
#endif

        mypMarker->setLabelAlignment( Qt::AlignHCenter | Qt::AlignTop );
        QwtText text( QString::number( currentInfo.id ) );
        text.setColor( myColor );
        text.setBackgroundBrush( textBgBrush );
        mypMarker->setLabel( text );
        mypMarker->attach( mpSatellitesWidget );
        mMarkerList << mypMarker;
#endif
      } // currentInfo.id > 0
    } //satellites
  } //satellite processing loop

  if ( mStackedWidget->currentIndex() == 1 && info.satInfoComplete ) //signal
  {
    mpCurve->setSamples( data );
    mpPlot->replot();
  } //signal
#ifdef WITH_QWTPOLAR
  if ( mStackedWidget->currentIndex() == 2 && info.satInfoComplete ) //satellites
  {
    mpSatellitesWidget->replot();
  } //satellites
#endif
  if ( validFlag )
  {
    validFlag = info.longitude >= -180.0 && info.longitude <= 180.0 && info.latitude >= -90.0 && info.latitude <= 90.0;
  }

  QgsPointXY myNewCenter;
  nmeaPOS newNmeaPosition;
  if ( validFlag )
  {
    myNewCenter = QgsPointXY( info.longitude, info.latitude );
    newNmeaPosition.lat = nmea_degree2radian( info.latitude );
    newNmeaPosition.lon = nmea_degree2radian( info.longitude );
  }
  else
  {
    myNewCenter = mLastGpsPosition;
    newNmeaPosition = mLastNmeaPosition;
  }
  if ( !mAcquisitionEnabled || ( nmea_distance( &newNmeaPosition, &mLastNmeaPosition ) < mDistanceThreshold ) )
  {
    // do not update position if update is disabled by timer or distance is under threshold
    myNewCenter = mLastGpsPosition;

  }
  if ( validFlag && mAcquisitionEnabled )
  {
    // position updated by valid data, reset timer
    switchAcquisition();
  }
  if ( mStackedWidget->currentIndex() == 0 ) //position
  {
    mTxtLatitude->setText( QString::number( info.latitude, 'f', 8 ) );
    mTxtLongitude->setText( QString::number( info.longitude, 'f', 8 ) );
    mTxtAltitude->setText( tr( "%1 m" ).arg( info.elevation, 0, 'f', 1 ) ); // don't know of any GPS receivers that output better than 0.1 m precision
    if ( mDateTimeFormat.isEmpty() )
    {
      mTxtDateTime->setText( info.utcDateTime.toString( Qt::TextDate ) );  // default format
    }
    else
    {
      mTxtDateTime->setText( info.utcDateTime.toString( mDateTimeFormat ) );  //user specified format string for testing the millisecond part of time
    }
    if ( std::isfinite( info.speed ) )
    {
      mTxtSpeed->setEnabled( true );
      mTxtSpeed->setText( tr( "%1 km/h" ).arg( info.speed, 0, 'f', 1 ) );
    }
    else
    {
      mTxtSpeed->setEnabled( false );
      mTxtSpeed->setText( tr( "Not available" ) );
    }
    if ( std::isfinite( info.direction ) )
    {
      mTxtDirection->setEnabled( true );
      mTxtDirection->setText( QString::number( info.direction, 'f', 1 ) + QStringLiteral( "??" ) );
    }
    else
    {
      mTxtDirection->setEnabled( false );
      mTxtDirection->setText( tr( "Not available" ) );
    }
    mTxtHdop->setText( QString::number( info.hdop, 'f', 1 ) );
    mTxtVdop->setText( QString::number( info.vdop, 'f', 1 ) );
    mTxtPdop->setText( QString::number( info.pdop, 'f', 1 ) );
    if ( std::isfinite( info.hacc ) )
    {
      mTxtHacc->setEnabled( true );
      mTxtHacc->setText( QString::number( info.hacc, 'f', 1 ) + "m" );
    }
    else
    {
      mTxtHacc->setEnabled( false );
      mTxtHacc->setText( tr( "Not available" ) );
    }
    if ( std::isfinite( info.vacc ) )
    {
      mTxtVacc->setEnabled( true );
      mTxtVacc->setText( QString::number( info.vacc, 'f', 1 ) + "m" );
    }
    else
    {
      mTxtVacc->setEnabled( false );
      mTxtVacc->setText( tr( "Not available" ) );
    }
    mTxtFixMode->setText( info.fixMode == 'A' ? tr( "Automatic" ) : info.fixMode == 'M' ? tr( "Manual" ) : QString() ); // A=automatic 2d/3d, M=manual; allowing for anything else
    mTxtFixType->setText( info.fixType == 3 ? tr( "3D" ) : info.fixType == 2 ? tr( "2D" ) : info.fixType == 1 ? tr( "No fix" ) : QString::number( info.fixType ) ); // 1=no fix, 2=2D, 3=3D; allowing for anything else
    mTxtQuality->setText( info.quality == 2 ? tr( "Differential" ) : info.quality == 1 ? tr( "Non-differential" ) : info.quality == 0 ? tr( "No position" ) : info.quality > 2 ? QString::number( info.quality ) : QString() ); // allowing for anything else
    mTxtSatellitesUsed->setText( QString::number( info.satellitesUsed ) );
    mTxtStatus->setText( info.status == 'A' ? tr( "Valid" ) : info.status == 'V' ? tr( "Invalid" ) : QString() );
  } //position

  // Avoid refreshing / panning if we haven't moved
  if ( mLastGpsPosition != myNewCenter )
  {
    mLastGpsPosition = myNewCenter;
    mLastNmeaPosition = newNmeaPosition;
    // Pan based on user specified behavior
    if ( radRecenterMap->isChecked() || radRecenterWhenNeeded->isChecked() )
    {
      QgsCoordinateReferenceSystem mypSRS = mpCanvas->mapSettings().destinationCrs();
      QgsCoordinateTransform myTransform( mWgs84CRS, mypSRS, QgsProject::instance() ); // use existing WGS84 CRS

      QgsPointXY myPoint = myTransform.transform( myNewCenter );
      //keep the extent the same just center the map canvas in the display so our feature is in the middle
      QgsRectangle myRect( myPoint, myPoint );  // empty rect can be used to set new extent that is centered on the point used to construct the rect

      // testing if position is outside some proportion of the map extent
      // this is a user setting - useful range: 5% to 100% (0.05 to 1.0)
      QgsRectangle myExtentLimit( mpCanvas->extent() );
      myExtentLimit.scale( mSpinMapExtentMultiplier->value() * 0.01 );

      // only change the extents if the point is beyond the current extents to minimize repaints
      if ( radRecenterMap->isChecked() ||
           ( radRecenterWhenNeeded->isChecked() && !myExtentLimit.contains( myPoint ) ) )
      {
        mpCanvas->setExtent( myRect );
        mpCanvas->refresh();
      }
    } //otherwise never recenter automatically

    if ( mCbxAutoAddVertices->isChecked() )
    {
      addVertex();
    }
  } // mLastGpsPosition != myNewCenter

  // new marker position after recentering
  if ( mGroupShowMarker->isChecked() ) // show marker
  {
    if ( validFlag ) // update cursor position if valid position
    {
      // initially, cursor isn't drawn until first valid fix; remains visible until GPS disconnect
      if ( ! mpMapMarker )
      {
        mpMapMarker = new QgsGpsMarker( mpCanvas );
      }
      mpMapMarker->setSize( mSliderMarkerSize->value() );
      mpMapMarker->setCenter( myNewCenter );
    }
  }
  else
  {
    if ( mpMapMarker )
    {
      delete mpMapMarker;
      mpMapMarker = nullptr;
    }
  } // show marker
}

void QgsGpsInformationWidget::mBtnAddVertex_clicked()
{
  addVertex();
}

void QgsGpsInformationWidget::addVertex()
{
  QgsDebugMsg( QStringLiteral( "Adding Vertex" ) );

  if ( !mpRubberBand )
  {
    createRubberBand();
  }

  // we store the capture list in wgs84 and then transform to layer crs when
  // calling close feature
  mCaptureList.push_back( mLastGpsPosition );

  // we store the rubber band points in map canvas CRS so transform to map crs
  // potential problem with transform errors and wrong coordinates if map CRS is changed after points are stored - SLM
  // should catch map CRS change and transform the points
  QgsPointXY myPoint;
  if ( mpCanvas )
  {
    QgsCoordinateTransform t( mWgs84CRS, mpCanvas->mapSettings().destinationCrs(), QgsProject::instance() );
    myPoint = t.transform( mLastGpsPosition );
  }
  else
  {
    myPoint = mLastGpsPosition;
  }

  mpRubberBand->addPoint( myPoint );
}

void QgsGpsInformationWidget::mBtnResetFeature_clicked()
{
  mNmea->disconnect( this, SLOT( displayGPSInformation( const QgsGpsInformation & ) ) );
  createRubberBand(); //deletes existing rubberband
  mCaptureList.clear();
  connectGpsSlot();
}

void QgsGpsInformationWidget::mBtnCloseFeature_clicked()
{
  QgsVectorLayer *vlayer = qobject_cast<QgsVectorLayer *>( mpCanvas->currentLayer() );
  QgsWkbTypes::Type layerWKBType = vlayer->wkbType();

  // -------------- preconditions ------------------------
  // most of these preconditions are already handled due to the button being enabled/disabled based on layer geom type and editing capabilities, but not on valid GPS data

  //lines: bail out if there are not at least two vertices
  if ( layerWKBType == QgsWkbTypes::LineString  && mCaptureList.size() < 2 )
  {
    QMessageBox::information( nullptr, tr( "Add Feature" ),
                              tr( "Cannot close a line feature until it has at least two vertices." ) );
    return;
  }

  //polygons: bail out if there are not at least three vertices
  if ( layerWKBType == QgsWkbTypes::Polygon && mCaptureList.size() < 3 )
  {
    QMessageBox::information( nullptr, tr( "Add Feature" ),
                              tr( "Cannot close a polygon feature until it has at least three vertices." ) );
    return;
  }
  // -------------- end of preconditions ------------------------

  //
  // POINT CAPTURING
  //
  if ( layerWKBType == QgsWkbTypes::Point )
  {
    QgsFeature *f = new QgsFeature( 0 );

    QgsCoordinateTransform t( mWgs84CRS, vlayer->crs(), QgsProject::instance() );
    QgsPointXY myPoint = t.transform( mLastGpsPosition );
    double x = myPoint.x();
    double y = myPoint.y();

    int size = 1 + sizeof( int ) + 2 * sizeof( double );
    unsigned char *buf = new unsigned char[size];

    QgsWkbPtr wkbPtr( buf, size );
    wkbPtr << ( char ) QgsApplication::endian() << QgsWkbTypes::Point << x << y;

    QgsGeometry g;
    g.fromWkb( buf, size );
    f->setGeometry( g );

    QgsFeatureAction action( tr( "Feature added" ), *f, vlayer, QString(), -1, this );
    if ( action.addFeature() )
    {
      if ( mCbxAutoCommit->isChecked() )
      {
        // should canvas->isDrawing() be checked?
        if ( !vlayer->commitChanges() ) //assumed to be vector layer and is editable and is in editing mode (preconditions have been tested)
        {
          QMessageBox::information( this,
                                    tr( "Save Layer Edits" ),
                                    tr( "Could not commit changes to layer %1\n\nErrors: %2\n" )
                                    .arg( vlayer->name(),
                                          vlayer->commitErrors().join( QStringLiteral( "\n  " ) ) ) );
        }

        vlayer->startEditing();
      }
    }

    delete f;
  } // layerWKBType == QgsWkbTypes::Point
  else // Line or poly
  {
    mNmea->disconnect( this, SLOT( displayGPSInformation( const QgsGpsInformation & ) ) );

    //create QgsFeature with wkb representation
    QgsFeature *f = new QgsFeature( 0 );

    if ( layerWKBType == QgsWkbTypes::LineString )
    {
      int size = 1 + 2 * sizeof( int ) + 2 * mCaptureList.size() * sizeof( double );
      unsigned char *buf = new unsigned char[size];

      QgsWkbPtr wkbPtr( buf, size );
      wkbPtr << ( char ) QgsApplication::endian() << QgsWkbTypes::LineString << mCaptureList.size();

      for ( QList<QgsPointXY>::const_iterator it = mCaptureList.constBegin(); it != mCaptureList.constEnd(); ++it )
      {
        QgsPointXY savePoint = *it;
        // transform the gps point into the layer crs
        QgsCoordinateTransform t( mWgs84CRS, vlayer->crs(), QgsProject::instance() );
        QgsPointXY myPoint = t.transform( savePoint );

        wkbPtr << myPoint.x() << myPoint.y();
      }

      QgsGeometry g;
      g.fromWkb( buf, size );
      f->setGeometry( g );
    }
    else if ( layerWKBType == QgsWkbTypes::Polygon )
    {
      int size = 1 + 3 * sizeof( int ) + 2 * ( mCaptureList.size() + 1 ) * sizeof( double );
      unsigned char *buf = new unsigned char[size];

      QgsWkbPtr wkbPtr( buf, size );
      wkbPtr << ( char ) QgsApplication::endian() << QgsWkbTypes::Polygon << 1 << mCaptureList.size() + 1;

      QList<QgsPointXY>::iterator it;
      for ( it = mCaptureList.begin(); it != mCaptureList.end(); ++it )
      {
        QgsPointXY savePoint = *it;
        // transform the gps point into the layer crs
        QgsCoordinateTransform t( mWgs84CRS, vlayer->crs(), QgsProject::instance() );
        QgsPointXY myPoint = t.transform( savePoint );
        wkbPtr << myPoint.x() << myPoint.y();
      }
      // close the polygon
      it = mCaptureList.begin();
      QgsPointXY savePoint = *it;

      wkbPtr << savePoint.x() << savePoint.y();

      QgsGeometry g;
      g.fromWkb( buf, size );
      f->setGeometry( g );

      QgsGeometry featGeom = f->geometry();
      int avoidIntersectionsReturn = featGeom.avoidIntersections( QgsProject::instance()->avoidIntersectionsLayers() );
      f->setGeometry( featGeom );
      if ( avoidIntersectionsReturn == 1 )
      {
        //not a polygon type. Impossible to get there
      }
      else if ( avoidIntersectionsReturn == 2 )
      {
        //bail out...
        QMessageBox::critical( nullptr, tr( "Add Feature" ), tr( "The feature could not be added because removing the polygon intersections would change the geometry type." ) );
        delete f;
        connectGpsSlot();
        return;
      }
      else if ( avoidIntersectionsReturn == 3 )
      {
        QMessageBox::critical( nullptr, tr( "Add Feature" ), tr( "An error was reported during intersection removal." ) );
        delete f;
        connectGpsSlot();
        return;
      }
    }
    // Should never get here, as preconditions should have removed any that aren't handled
    else // layerWKBType == QgsWkbTypes::Polygon  -  unknown type
    {
      QMessageBox::critical( nullptr, tr( "Add Feature" ), tr( "Cannot add feature. "
                             "Unknown WKB type. Choose a different layer and try again." ) );
      connectGpsSlot();
      delete f;
      return; //unknown wkbtype
    } // layerWKBType == QgsWkbTypes::Polygon

    QgsFeatureAction action( tr( "Feature added" ), *f, vlayer, QString(), -1, this );
    if ( action.addFeature() )
    {
      if ( mCbxAutoCommit->isChecked() )
      {
        if ( !vlayer->commitChanges() ) //swiped... er... appropriated from QgisApp saveEdits()
        {
          QMessageBox::information( this,
                                    tr( "Save Layer Edits" ),
                                    tr( "Could not commit changes to layer %1\n\nErrors: %2\n" )
                                    .arg( vlayer->name(),
                                          vlayer->commitErrors().join( QStringLiteral( "\n  " ) ) ) );
        }

        vlayer->startEditing();
      }
      delete mpRubberBand;
      mpRubberBand = nullptr;

      // delete the elements of mCaptureList
      mCaptureList.clear();
    } // action.addFeature()

    delete f;
    connectGpsSlot();
  } // layerWKBType == QgsWkbTypes::Point
  mpCanvas->refresh();  // NOTE: canceling feature add refreshes canvas, OK does not; this may change, however, so do it anyway

  // force focus back to GPS window/ Add Feature button for ease of use by keyboard
  activateWindow();
  mBtnCloseFeature->setFocus( Qt::OtherFocusReason );
}

void QgsGpsInformationWidget::connectGpsSlot()
{
  connect( mNmea, &QgsGpsConnection::stateChanged,
           this, &QgsGpsInformationWidget::displayGPSInformation );
}

void QgsGpsInformationWidget::mBtnRefreshDevices_clicked()
{
  populateDevices();
}

/* Copied from gps plugin */
void QgsGpsInformationWidget::populateDevices()
{
  QList< QPair<QString, QString> > ports = QgsGpsDetector::availablePorts();

  mCboDevices->clear();

  // add devices to combobox, but skip gpsd which is first.
  for ( int i = 1; i < ports.size(); i++ )
  {
    mCboDevices->addItem( ports[i].second, ports[i].first );
  }

  // remember the last ports used
  QgsSettings settings;
  QString lastPort = settings.value( QStringLiteral( "gps/lastPort" ), "" ).toString();

  int idx = mCboDevices->findData( lastPort );
  mCboDevices->setCurrentIndex( idx < 0 ? 0 : idx );
}

void QgsGpsInformationWidget::createRubberBand()
{
  delete mpRubberBand;

  mpRubberBand = new QgsRubberBand( mpCanvas, QgsWkbTypes::LineGeometry );
  mpRubberBand->setColor( mBtnTrackColor->color() );
  mpRubberBand->setWidth( mSpinTrackWidth->value() );
  mpRubberBand->show();
}

void QgsGpsInformationWidget::mBtnLogFile_clicked()
{
//=========================
  // This does not allow for an extension other than ".nmea"
  // Retrieve last used log file dir from persistent settings
  QgsSettings settings;
  QString settingPath( QStringLiteral( "/gps/lastLogFileDir" ) );
  QString lastUsedDir = settings.value( settingPath, QDir::homePath() ).toString();
  QString saveFilePath = QFileDialog::getSaveFileName( this, tr( "Save GPS log file As" ), lastUsedDir, tr( "NMEA files" ) + " (*.nmea)" );
  if ( saveFilePath.isNull() ) //canceled
  {
    return;
  }
  QFileInfo myFI( saveFilePath );
  QString myPath = myFI.path();
  settings.setValue( settingPath, myPath );

  // make sure the .nmea extension is included in the path name. if not, add it...
  if ( "nmea" != myFI.suffix() )
  {
    saveFilePath = myFI.filePath() + ".nmea";
  }
  mTxtLogFile->setText( saveFilePath );
  mTxtLogFile->setToolTip( saveFilePath );
}

void QgsGpsInformationWidget::logNmeaSentence( const QString &nmeaString )
{
  if ( mLogFileGroupBox->isChecked() && mLogFile && mLogFile->isOpen() )
  {
    mLogFileTextStream << nmeaString << "\r\n"; // specifically output CR + LF (NMEA requirement)
  }
}

void QgsGpsInformationWidget::updateCloseFeatureButton( QgsMapLayer *lyr )
{
  QgsVectorLayer *vlayer = qobject_cast<QgsVectorLayer *>( lyr );

  if ( !( vlayer && vlayer->isValid() ) )
    return;

  // Add feature button tracks edit state of layer
  if ( vlayer != mpLastLayer )
  {
    if ( mpLastLayer )  // disconnect previous layer
    {
      disconnect( mpLastLayer, &QgsVectorLayer::editingStarted,
                  this, &QgsGpsInformationWidget::layerEditStateChanged );
      disconnect( mpLastLayer, &QgsVectorLayer::editingStopped,
                  this, &QgsGpsInformationWidget::layerEditStateChanged );
    }
    if ( vlayer ) // connect new layer
    {
      connect( vlayer, &QgsVectorLayer::editingStarted,
               this, &QgsGpsInformationWidget::layerEditStateChanged );
      connect( vlayer, &QgsVectorLayer::editingStopped,
               this, &QgsGpsInformationWidget::layerEditStateChanged );
    }
    mpLastLayer = vlayer;
  }

  QString buttonLabel = tr( "&Add feature" );
  if ( vlayer ) // must be vector layer
  {
    QgsVectorDataProvider *provider = vlayer->dataProvider();
    QgsWkbTypes::Type layerWKBType = vlayer->wkbType();

    QgsWkbTypes::Type flatType = QgsWkbTypes::flatType( layerWKBType );

    bool enable =
      ( provider->capabilities() & QgsVectorDataProvider::AddFeatures ) &&  // layer can add features
      vlayer->isEditable() && // layer is editing
      ( // layer has geometry type that can be handled
        flatType == QgsWkbTypes::Point ||
        flatType == QgsWkbTypes::LineString ||
        flatType == QgsWkbTypes::Polygon
        // add more types here as they are handled
      )
      ;

    if ( flatType == QgsWkbTypes::Point )
      buttonLabel = tr( "&Add Point" );
    else if ( flatType == QgsWkbTypes::LineString )
      buttonLabel = tr( "&Add Line" );
    else if ( flatType == QgsWkbTypes::Polygon )
      buttonLabel = tr( "&Add Polygon" );
    // TODO: Add multi types

    mBtnCloseFeature->setEnabled( enable );
  }
  else
  {
    mBtnCloseFeature->setEnabled( false );
  }
  mBtnCloseFeature->setText( buttonLabel );
}

void QgsGpsInformationWidget::layerEditStateChanged()
{
  updateCloseFeatureButton( mpLastLayer );
}

void QgsGpsInformationWidget::setStatusIndicator( const FixStatus statusValue )
{
  mLastFixStatus = statusValue;
  // the pixmap will be expanded to the size of the label
  QPixmap status( 4, 4 );
  switch ( statusValue )
  {
    case NoFix:
      status.fill( Qt::red );
      break;
    case Fix2D:
      status.fill( Qt::yellow );
      break;
    case Fix3D:
      status.fill( Qt::green );
      break;
    case NoData:
    default: // anything else - shouldn't happen
      status.fill( Qt::darkGray );
  }
  mLblStatusIndicator->setPixmap( status );
}

void QgsGpsInformationWidget::showStatusBarMessage( const QString &msg )
{
  QgisApp::instance()->statusBarIface()->showMessage( msg );
}
void QgsGpsInformationWidget::setAcquisitionInterval( uint interval )
{
  mAcquisitionInterval = interval * 1000;
  if ( mAcquisitionTimer->isActive() )
    mAcquisitionTimer->stop();
  mAcquisitionEnabled = true;
  switchAcquisition();

}
void QgsGpsInformationWidget::setDistanceThreshold( uint distance )
{
  mDistanceThreshold = distance;
}

void QgsGpsInformationWidget::cboAcquisitionIntervalEdited()
{
  setAcquisitionInterval( mCboAcquisitionInterval->currentText().toUInt() );
}

void QgsGpsInformationWidget::cboDistanceThresholdEdited()
{
  setDistanceThreshold( mCboDistanceThreshold->currentText().toUInt() );
}

void QgsGpsInformationWidget::switchAcquisition()
{
  if ( mAcquisitionInterval > 0 )
  {
    if ( mAcquisitionEnabled )
      mAcquisitionTimer->start( mAcquisitionInterval );
    else
      //wait only acquisitionInterval/10 for new valid data
      mAcquisitionTimer->start( mAcquisitionInterval / 10 );
    // anyway switch to enabled / disabled acquisition
    mAcquisitionEnabled = !mAcquisitionEnabled;
  }
}
