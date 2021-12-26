/***************************************************************************
    qgsowssourceselect.cpp  -  selector for WMS,WFS,WCS
                             -------------------
    begin                : 3 April 2005
    copyright            :
    original             : (C) 2005 by Brendan Morley email  : morb at ozemail dot com dot au
    wms search           : (C) 2009 Mathias Walker <mwa at sourcepole.ch>, Sourcepole AG
    wms-c support        : (C) 2010 Juergen E. Fischer < jef at norbit dot de >, norBIT GmbH

    generalized          : (C) 2012 Radim Blazek, based on qgswmssourceselect.cpp

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgis.h" // GEO_EPSG_CRS_ID
#include "qgscoordinatereferencesystem.h"
#include "qgsdatasourceuri.h"
#include "qgsprojectionselectiondialog.h"
#include "qgslogger.h"
#include "qgsmanageconnectionsdialog.h"
#include "qgsmessageviewer.h"
#include "qgsnewhttpconnection.h"
#include "qgstreewidgetitem.h"
#include "qgsproject.h"
#include "qgsproviderregistry.h"
#include "qgsowsconnection.h"
#include "qgsdataprovider.h"
#include "qgsowssourceselect.h"
#include "qgsnetworkaccessmanager.h"
#include "qgsapplication.h"
#include "qgssettings.h"

#include <QButtonGroup>
#include <QFileDialog>
#include <QRadioButton>
#include <QDomDocument>
#include <QHeaderView>
#include <QImageReader>
#include <QInputDialog>
#include <QMap>
#include <QMessageBox>
#include <QPicture>
#include <QUrl>
#include <QValidator>
#include <QNetworkRequest>
#include <QNetworkReply>

QgsOWSSourceSelect::QgsOWSSourceSelect( const QString &service, QWidget *parent, Qt::WindowFlags fl, QgsProviderRegistry::WidgetMode theWidgetMode )
  : QgsAbstractDataSourceWidget( parent, fl, theWidgetMode )
  , mService( service )

{
  setupUi( this );
  connect( mNewButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mNewButton_clicked );
  connect( mEditButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mEditButton_clicked );
  connect( mDeleteButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mDeleteButton_clicked );
  connect( mSaveButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mSaveButton_clicked );
  connect( mLoadButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mLoadButton_clicked );
  connect( mConnectButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mConnectButton_clicked );
  connect( mChangeCRSButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mChangeCRSButton_clicked );
  connect( mLayersTreeWidget, &QTreeWidget::itemSelectionChanged, this, &QgsOWSSourceSelect::mLayersTreeWidget_itemSelectionChanged );
  connect( mConnectionsComboBox, static_cast<void ( QComboBox::* )( int )>( &QComboBox::activated ), this, &QgsOWSSourceSelect::mConnectionsComboBox_activated );
  connect( mAddDefaultButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mAddDefaultButton_clicked );
  connect( mSearchButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mSearchButton_clicked );
  connect( mSearchTableWidget, &QTableWidget::itemSelectionChanged, this, &QgsOWSSourceSelect::mSearchTableWidget_itemSelectionChanged );
  connect( mTilesetsTableWidget, &QTableWidget::itemClicked, this, &QgsOWSSourceSelect::mTilesetsTableWidget_itemClicked );
  connect( mLayerUpButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mLayerUpButton_clicked );
  connect( mLayerDownButton, &QPushButton::clicked, this, &QgsOWSSourceSelect::mLayerDownButton_clicked );
  setupButtons( buttonBox );


  setWindowTitle( tr( "Add Layer(s) from a %1 Server" ).arg( service ) );

  clearCrs();

  mTileWidthLineEdit->setValidator( new QIntValidator( 0, 9999, this ) );
  mTileHeightLineEdit->setValidator( new QIntValidator( 0, 9999, this ) );
  mFeatureCountLineEdit->setValidator( new QIntValidator( 0, 9999, this ) );

  mCacheComboBox->addItem( tr( "Always cache" ), QNetworkRequest::AlwaysCache );
  mCacheComboBox->addItem( tr( "Prefer cache" ), QNetworkRequest::PreferCache );
  mCacheComboBox->addItem( tr( "Prefer network" ), QNetworkRequest::PreferNetwork );
  mCacheComboBox->addItem( tr( "Always network" ), QNetworkRequest::AlwaysNetwork );

  // 'Prefer network' is the default noted in the combobox's tool tip
  mCacheComboBox->setCurrentIndex( mCacheComboBox->findData( QNetworkRequest::PreferNetwork ) );

  if ( widgetMode() != QgsProviderRegistry::WidgetMode::Manager )
  {
    //set the current project CRS if available
    QgsCoordinateReferenceSystem currentRefSys = QgsProject::instance()->crs();
    //convert CRS id to epsg
    if ( currentRefSys.isValid() )
    {
      mSelectedCRS = currentRefSys.authid();
    }
  }
  else
  {
    mTabWidget->removeTab( mTabWidget->indexOf( mLayerOrderTab ) );
    mTabWidget->removeTab( mTabWidget->indexOf( mTilesetsTab ) );
    mTimeWidget->hide();
    mFormatWidget->hide();
    mCRSWidget->hide();
    mCacheWidget->hide();
  }

  // set up the WMS connections we already know about
  populateConnectionList();

  QgsSettings settings;
  QgsDebugMsg( QStringLiteral( "restoring geometry" ) );
  restoreGeometry( settings.value( QStringLiteral( "Windows/WMSSourceSelect/geometry" ) ).toByteArray() );
}

QgsOWSSourceSelect::~QgsOWSSourceSelect()
{
  QgsSettings settings;
  QgsDebugMsg( QStringLiteral( "saving geometry" ) );
  settings.setValue( QStringLiteral( "Windows/WMSSourceSelect/geometry" ), saveGeometry() );
}

void QgsOWSSourceSelect::refresh()
{
  populateConnectionList();
}

void QgsOWSSourceSelect::clearFormats()
{
  mFormatComboBox->clear();
  mFormatComboBox->setEnabled( false );
}

void QgsOWSSourceSelect::populateFormats()
{

  // A server may offer more similar formats, which are mapped
  // to the same GDAL format, e.g. GeoTIFF and TIFF
  // -> recreate always buttons for all available formats, enable supported

  clearFormats();

  if ( mProviderFormats.isEmpty() )
  {
    mProviderFormats = providerFormats();
    for ( int i = 0; i < mProviderFormats.size(); i++ )
    {
      // GDAL mime types may be image/tiff, image/png, ...
      mMimeLabelMap.insert( mProviderFormats[i].format, mProviderFormats[i].label );
    }
  }

  // selectedLayersFormats may come in various forms:
  // image/tiff, GTiff, GeoTIFF, TIFF, geotiff_int16, geotiff_rgb,
  // PNG, GTOPO30, ARCGRID, IMAGEMOSAIC
  // and even any string defined in server configuration, for example the
  // value used in UMN Mapserver for OUTPUTFORMAT->NAME is used in
  // WCS 1.0.0 SupportedFormats/Format

  // TODO: It is impossible to cover all possible formats coming from server
  //       -> enabled all formats, GDAL may be able to open them

  QMap<QString, QString> formatsMap;
  formatsMap.insert( QStringLiteral( "geotiff" ), QStringLiteral( "tiff" ) );
  formatsMap.insert( QStringLiteral( "gtiff" ), QStringLiteral( "tiff" ) );
  formatsMap.insert( QStringLiteral( "tiff" ), QStringLiteral( "tiff" ) );
  formatsMap.insert( QStringLiteral( "tif" ), QStringLiteral( "tiff" ) );
  formatsMap.insert( QStringLiteral( "gif" ), QStringLiteral( "gif" ) );
  formatsMap.insert( QStringLiteral( "jpeg" ), QStringLiteral( "jpeg" ) );
  formatsMap.insert( QStringLiteral( "jpg" ), QStringLiteral( "jpeg" ) );
  formatsMap.insert( QStringLiteral( "png" ), QStringLiteral( "png" ) );

  int preferred = -1;
  QStringList layersFormats = selectedLayersFormats();
  for ( int i = 0; i < layersFormats.size(); i++ )
  {
    QString format = layersFormats.value( i );
    QgsDebugMsg( "server format = " + format );
    QString simpleFormat = format.toLower().remove( QStringLiteral( "image/" ) ).remove( QRegExp( "_.*" ) );
    QgsDebugMsg( "server simpleFormat = " + simpleFormat );
    QString mimeFormat = "image/" + formatsMap.value( simpleFormat );
    QgsDebugMsg( "server mimeFormat = " + mimeFormat );

    QString label = format;

    if ( mMimeLabelMap.contains( mimeFormat ) )
    {
      if ( format != mMimeLabelMap.value( mimeFormat ) )
      {
        // Append name of GDAL driver
        label += " / " + mMimeLabelMap.value( mimeFormat );
      }

      if ( simpleFormat.contains( QLatin1String( "tif" ) ) ) // prefer *tif*
      {
        if ( preferred < 0 || simpleFormat.startsWith( 'g' ) ) // prefer geotiff
        {
          preferred = i;
        }
      }
    }
    else
    {
      // We cannot always say that the format is not supported by GDAL because
      // server can use strange names, but format itself is supported
      QgsDebugMsg( QStringLiteral( "format %1 unknown" ).arg( format ) );
    }

    mFormatComboBox->insertItem( i, label );
  }
  // Set preferred
  // TODO: all enabled for now, see above
  preferred = preferred >= 0 ? preferred : 0;
  mFormatComboBox->setCurrentIndex( preferred );

  mFormatComboBox->setEnabled( true );
}

void QgsOWSSourceSelect::populateTimes()
{
  mTimeComboBox->clear();
  mTimeComboBox->insertItems( 0, selectedLayersTimes() );
  mTimeComboBox->setEnabled( !selectedLayersTimes().isEmpty() );
}

void QgsOWSSourceSelect::clearTimes()
{
  mTimeComboBox->clear();
  mTimeComboBox->setEnabled( false );
}

void QgsOWSSourceSelect::populateConnectionList()
{
  mConnectionsComboBox->clear();
  mConnectionsComboBox->addItems( QgsOwsConnection::connectionList( mService ) );

  setConnectionListPosition();
}

QgsNewHttpConnection::ConnectionType connectionTypeFromServiceString( const QString &string )
{
  if ( string.compare( QLatin1String( "wms" ), Qt::CaseInsensitive ) == 0 )
    return QgsNewHttpConnection::ConnectionWms;
  else if ( string.compare( QLatin1String( "wfs" ), Qt::CaseInsensitive ) == 0 )
    return QgsNewHttpConnection::ConnectionWfs;
  else if ( string.compare( QLatin1String( "wcs" ), Qt::CaseInsensitive ) == 0 )
    return QgsNewHttpConnection::ConnectionWcs;
  else
    return QgsNewHttpConnection::ConnectionWms;
}

void QgsOWSSourceSelect::mNewButton_clicked()
{
  QgsNewHttpConnection::ConnectionType type = connectionTypeFromServiceString( mService );
  QgsNewHttpConnection *nc = new QgsNewHttpConnection( this, type, "/qgis/connections-" + mService.toLower() + '/' );

  if ( nc->exec() )
  {
    populateConnectionList();
    emit connectionsChanged();
  }

  delete nc;
}

void QgsOWSSourceSelect::mEditButton_clicked()
{
  QgsNewHttpConnection::ConnectionType type = connectionTypeFromServiceString( mService );
  QgsNewHttpConnection *nc = new QgsNewHttpConnection( this, type, "/qgis/connections-" + mService.toLower() + '/', mConnectionsComboBox->currentText() );

  if ( nc->exec() )
  {
    populateConnectionList();
    emit connectionsChanged();
  }

  delete nc;
}

void QgsOWSSourceSelect::mDeleteButton_clicked()
{
  QString msg = tr( "Are you sure you want to remove the %1 connection and all associated settings?" )
                .arg( mConnectionsComboBox->currentText() );
  QMessageBox::StandardButton result = QMessageBox::question( this, tr( "Delete Connection" ), msg, QMessageBox::Yes | QMessageBox::No );
  if ( result == QMessageBox::Yes )
  {
    QgsOwsConnection::deleteConnection( mService, mConnectionsComboBox->currentText() );
    mConnectionsComboBox->removeItem( mConnectionsComboBox->currentIndex() );  // populateConnectionList();
    setConnectionListPosition();
    emit connectionsChanged();
  }
}

void QgsOWSSourceSelect::mSaveButton_clicked()
{
  QgsManageConnectionsDialog dlg( this, QgsManageConnectionsDialog::Export, QgsManageConnectionsDialog::WCS );
  dlg.exec();
}

void QgsOWSSourceSelect::mLoadButton_clicked()
{
  QString fileName = QFileDialog::getOpenFileName( this, tr( "Load Connections" ), QDir::homePath(),
                     tr( "XML files (*.xml *.XML)" ) );
  if ( fileName.isEmpty() )
  {
    return;
  }

  QgsManageConnectionsDialog dlg( this, QgsManageConnectionsDialog::Import, QgsManageConnectionsDialog::WCS, fileName );
  dlg.exec();
  populateConnectionList();
  emit connectionsChanged();
}

QgsTreeWidgetItem *QgsOWSSourceSelect::createItem(
  int id,
  const QStringList &names,
  QMap<int, QgsTreeWidgetItem *> &items,
  int &layerAndStyleCount,
  const QMap<int, int> &layerParents,
  const QMap<int, QStringList> &layerParentNames )
{
  QgsDebugMsg( QStringLiteral( "id = %1 layerAndStyleCount = %2 names = %3 " ).arg( id ).arg( layerAndStyleCount ).arg( names.join( "," ) ) );
  if ( items.contains( id ) )
    return items[id];


  QgsTreeWidgetItem *item = nullptr;
  if ( layerParents.contains( id ) )
  {
    // it has parent -> create first its parent
    int parent = layerParents[ id ];
    item = new QgsTreeWidgetItem( createItem( parent, layerParentNames[ parent ], items, layerAndStyleCount, layerParents, layerParentNames ) );
  }
  else
    item = new QgsTreeWidgetItem( mLayersTreeWidget );

  item->setText( 0, QString::number( ++layerAndStyleCount ) );
  item->setText( 1, names[0].simplified() );
  item->setText( 2, names[1].simplified() );
  item->setText( 3, names[2].simplified() );
  item->setToolTip( 3, "<font color=black>" + names[2].simplified()  + "</font>" );

  items[ id ] = item;

  return item;
}

void QgsOWSSourceSelect::populateLayerList()
{
}

void QgsOWSSourceSelect::mConnectButton_clicked()
{

  mLayersTreeWidget->clear();
  clearFormats();
  clearTimes();
  clearCrs();

  mConnName = mConnectionsComboBox->currentText();

  QgsOwsConnection connection( mService, mConnectionsComboBox->currentText() );
  mUri = connection.uri();

  QApplication::setOverrideCursor( Qt::WaitCursor );

  QgsDebugMsg( QStringLiteral( "call populateLayerList" ) );
  populateLayerList();

  QApplication::restoreOverrideCursor();
}

void QgsOWSSourceSelect::enableLayersForCrs( QTreeWidgetItem * )
{
}

void QgsOWSSourceSelect::mChangeCRSButton_clicked()
{
  QStringList layers;
  const auto constSelectedItems = mLayersTreeWidget->selectedItems();
  for ( QTreeWidgetItem *item : constSelectedItems )
  {
    QString layer = item->data( 0, Qt::UserRole + 0 ).toString();
    if ( !layer.isEmpty() )
      layers << layer;
  }

  QgsProjectionSelectionDialog *mySelector = new QgsProjectionSelectionDialog( this );
  mySelector->setMessage( QString() );
  mySelector->setOgcWmsCrsFilter( mSelectedLayersCRSs );

  QgsCoordinateReferenceSystem defaultCRS = QgsProject::instance()->crs();
  if ( defaultCRS.isValid() )
  {
    mySelector->setCrs( defaultCRS );
  }

  if ( !mySelector->exec() )
    return;

  mSelectedCRS = mySelector->crs().authid();
  delete mySelector;

  mSelectedCRSLabel->setText( descriptionForAuthId( mSelectedCRS ) );

  for ( int i = 0; i < mLayersTreeWidget->topLevelItemCount(); i++ )
  {
    enableLayersForCrs( mLayersTreeWidget->topLevelItem( i ) );
  }

  updateButtons();
}

void QgsOWSSourceSelect::mLayersTreeWidget_itemSelectionChanged()
{
}

void QgsOWSSourceSelect::populateCrs()
{
  clearCrs();
  mSelectedLayersCRSs = selectedLayersCrses().toSet();
  mCRSLabel->setText( tr( "Coordinate Reference System (%n available)", "crs count", mSelectedLayersCRSs.count() ) + ':' );

  mChangeCRSButton->setDisabled( mSelectedLayersCRSs.isEmpty() );

  if ( !mSelectedLayersCRSs.isEmpty() )
  {
    // check whether current CRS is supported
    // if not, use one of the available CRS
    QString defaultCRS;
    QSet<QString>::const_iterator it = mSelectedLayersCRSs.constBegin();
    for ( ; it != mSelectedLayersCRSs.constEnd(); ++it )
    {
      if ( it->compare( mSelectedCRS, Qt::CaseInsensitive ) == 0 )
        break;

      // save first CRS in case the current CRS is not available
      if ( it == mSelectedLayersCRSs.constBegin() )
        defaultCRS = *it;

      // prefer value of DEFAULT_GEO_EPSG_CRS_ID if available
      if ( *it == GEO_EPSG_CRS_AUTHID )
        defaultCRS = *it;
    }

    if ( it == mSelectedLayersCRSs.constEnd() )
    {
      // not found
      mSelectedCRS = defaultCRS;
    }
    mSelectedCRSLabel->setText( descriptionForAuthId( mSelectedCRS ) );
    mChangeCRSButton->setEnabled( true );
  }
  QgsDebugMsg( "mSelectedCRS = " + mSelectedCRS );
}

void QgsOWSSourceSelect::clearCrs()
{
  mCRSLabel->setText( tr( "Coordinate Reference System" ) + ':' );
  mSelectedCRS.clear();
  mSelectedCRSLabel->clear();
  mChangeCRSButton->setEnabled( false );
}

void QgsOWSSourceSelect::mTilesetsTableWidget_itemClicked( QTableWidgetItem *item )
{
  Q_UNUSED( item );

  QTableWidgetItem *rowItem = mTilesetsTableWidget->item( mTilesetsTableWidget->currentRow(), 0 );
  bool wasSelected = mCurrentTileset == rowItem;

  mTilesetsTableWidget->blockSignals( true );
  mTilesetsTableWidget->clearSelection();
  if ( !wasSelected )
  {
    QgsDebugMsg( QStringLiteral( "selecting current row %1" ).arg( mTilesetsTableWidget->currentRow() ) );
    mTilesetsTableWidget->selectRow( mTilesetsTableWidget->currentRow() );
    mCurrentTileset = rowItem;
  }
  else
  {
    mCurrentTileset = nullptr;
  }
  mTilesetsTableWidget->blockSignals( false );

  updateButtons();
}



QString QgsOWSSourceSelect::connName()
{
  return mConnName;
}

QString QgsOWSSourceSelect::connectionInfo()
{
  return mConnectionInfo;
}

QString QgsOWSSourceSelect::selectedFormat()
{
  return selectedLayersFormats().value( mFormatComboBox->currentIndex() );
}

QNetworkRequest::CacheLoadControl QgsOWSSourceSelect::selectedCacheLoadControl()
{
  int cache = mCacheComboBox->currentData().toInt();
  return static_cast<QNetworkRequest::CacheLoadControl>( cache );
}

QString QgsOWSSourceSelect::selectedCrs()
{
  return mSelectedCRS;
}

QString QgsOWSSourceSelect::selectedTime()
{
  return mTimeComboBox->currentText();
}

void QgsOWSSourceSelect::setConnectionListPosition()
{
  QString toSelect = QgsOwsConnection::selectedConnection( mService );

  mConnectionsComboBox->setCurrentIndex( mConnectionsComboBox->findText( toSelect ) );

  if ( mConnectionsComboBox->currentIndex() < 0 )
  {
    if ( toSelect.isNull() )
      mConnectionsComboBox->setCurrentIndex( 0 );
    else
      mConnectionsComboBox->setCurrentIndex( mConnectionsComboBox->count() - 1 );
  }

  if ( mConnectionsComboBox->count() == 0 )
  {
    // No connections - disable various buttons
    mConnectButton->setEnabled( false );
    mEditButton->setEnabled( false );
    mDeleteButton->setEnabled( false );
    mSaveButton->setEnabled( false );
  }
  else
  {
    // Connections - enable various buttons
    mConnectButton->setEnabled( true );
    mEditButton->setEnabled( true );
    mDeleteButton->setEnabled( true );
    mSaveButton->setEnabled( true );
  }

  QgsOwsConnection::setSelectedConnection( mService, mConnectionsComboBox->currentText() );
}

void QgsOWSSourceSelect::showStatusMessage( QString const &message )
{
  mStatusLabel->setText( message );

  // update the display of this widget
  update();
}


void QgsOWSSourceSelect::showError( QString const &title, QString const &format, QString const &error )
{
  QgsMessageViewer *mv = new QgsMessageViewer( this );
  mv->setWindowTitle( title );

  if ( format == QLatin1String( "text/html" ) )
  {
    mv->setMessageAsHtml( error );
  }
  else
  {
    mv->setMessageAsPlainText( tr( "Could not understand the response:\n%1" ).arg( error ) );
  }
  mv->showMessage( true ); // Is deleted when closed
}

void QgsOWSSourceSelect::mConnectionsComboBox_activated( int )
{
  // Remember which server was selected.
  QgsOwsConnection::setSelectedConnection( mService, mConnectionsComboBox->currentText() );
}

void QgsOWSSourceSelect::mAddDefaultButton_clicked()
{
  addDefaultServers();
}

QString QgsOWSSourceSelect::descriptionForAuthId( const QString &authId )
{
  if ( mCrsNames.contains( authId ) )
    return mCrsNames[ authId ];

  QgsCoordinateReferenceSystem qgisSrs = QgsCoordinateReferenceSystem::fromOgcWmsCrs( authId );
  mCrsNames.insert( authId, qgisSrs.description() );
  return qgisSrs.description();
}

void QgsOWSSourceSelect::addDefaultServers()
{
  QMap<QString, QString> exampleServers;
  exampleServers[QStringLiteral( "DM Solutions GMap" )] = QStringLiteral( "http://www2.dmsolutions.ca/cgi-bin/mswms_gmap" );
  exampleServers[QStringLiteral( "Lizardtech server" )] = QStringLiteral( "http://wms.lizardtech.com/lizardtech/iserv/ows" );
  // Nice to have the qgis users map, but I'm not sure of the URL at the moment.
  //  exampleServers["Qgis users map"] = "http://qgis.org/wms.cgi";

  QgsSettings settings;
  settings.beginGroup( "/qgis/connections-" + mService.toLower() );
  QMap<QString, QString>::const_iterator i = exampleServers.constBegin();
  for ( ; i != exampleServers.constEnd(); ++i )
  {
    // Only do a server if it's name doesn't already exist.
    QStringList keys = settings.childGroups();
    if ( !keys.contains( i.key() ) )
    {
      QString path = i.key();
      settings.setValue( path + "/url", i.value() );
    }
  }
  settings.endGroup();
  populateConnectionList();

  QMessageBox::information( this, tr( "Add WMS Servers" ), "<p>" + tr( "Several WMS servers have "
                            "been added to the server list. Note that if "
                            "you access the Internet via a web proxy, you will "
                            "need to set the proxy settings in the QGIS options dialog." ) + "</p>" );
}

void QgsOWSSourceSelect::addWmsListRow( const QDomElement &item, int row )
{
  QDomElement title = item.firstChildElement( QStringLiteral( "title" ) );
  addWmsListItem( title, row, 0 );
  QDomElement description = item.firstChildElement( QStringLiteral( "description" ) );
  addWmsListItem( description, row, 1 );
  QDomElement link = item.firstChildElement( QStringLiteral( "link" ) );
  addWmsListItem( link, row, 2 );
}

void QgsOWSSourceSelect::addWmsListItem( const QDomElement &el, int row, int column )
{
  if ( !el.isNull() )
  {
    QTableWidgetItem *tableItem = new QTableWidgetItem( el.text() );
    // TODO: add linebreaks to long tooltips?
    tableItem->setToolTip( el.text() );
    mSearchTableWidget->setItem( row, column, tableItem );
  }
}

void QgsOWSSourceSelect::mSearchButton_clicked()
{
  // clear results
  mSearchTableWidget->clearContents();
  mSearchTableWidget->setRowCount( 0 );

  // disable Add WMS button
  mSearchAddButton->setEnabled( false );

  QApplication::setOverrideCursor( Qt::WaitCursor );

  QgsSettings settings;
  // geopole.org (geopole.ch) 25.4.2012 : 503 Service Unavailable, archive: Recently added 20 Jul 2011
  QString mySearchUrl = settings.value( QStringLiteral( "qgis/WMSSearchUrl" ), "http://geopole.org/wms/search?search=%1&type=rss" ).toString();
  QUrl url( mySearchUrl.arg( mSearchTermLineEdit->text() ) );
  QgsDebugMsg( url.toString() );

  QNetworkReply *r = QgsNetworkAccessManager::instance()->get( QNetworkRequest( url ) );
  connect( r, &QNetworkReply::finished, this, &QgsOWSSourceSelect::searchFinished );
}

void QgsOWSSourceSelect::searchFinished()
{
  QApplication::restoreOverrideCursor();

  QNetworkReply *r = qobject_cast<QNetworkReply *>( sender() );
  if ( !r )
    return;

  if ( r->error() == QNetworkReply::NoError )
  {
    // parse results
    QDomDocument doc( QStringLiteral( "RSS" ) );
    QByteArray res = r->readAll();
    QString error;
    int line, column;
    if ( doc.setContent( res, &error, &line, &column ) )
    {
      QDomNodeList list = doc.elementsByTagName( QStringLiteral( "item" ) );
      mSearchTableWidget->setRowCount( list.size() );
      for ( int i = 0; i < list.size(); i++ )
      {
        if ( list.item( i ).isElement() )
        {
          QDomElement item = list.item( i ).toElement();
          addWmsListRow( item, i );
        }
      }

      mSearchTableWidget->resizeColumnsToContents();
    }
    else
    {
      QgsDebugMsg( QStringLiteral( "setContent failed" ) );
      showStatusMessage( tr( "parse error at row %1, column %2: %3" ).arg( line ).arg( column ).arg( error ) );
    }
  }
  else
  {
    showStatusMessage( tr( "network error: %1" ).arg( r->error() ) );
  }

  r->deleteLater();
}

void QgsOWSSourceSelect::mSearchTableWidget_itemSelectionChanged()
{
  mSearchAddButton->setEnabled( mSearchTableWidget->currentRow() != -1 );
}

void QgsOWSSourceSelect::mLayerUpButton_clicked()
{
  QList<QTreeWidgetItem *> selectionList = mLayerOrderTreeWidget->selectedItems();
  if ( selectionList.empty() )
  {
    return;
  }
  int selectedIndex = mLayerOrderTreeWidget->indexOfTopLevelItem( selectionList[0] );
  if ( selectedIndex < 1 )
  {
    return; //item not existing or already on top
  }

  QTreeWidgetItem *selectedItem = mLayerOrderTreeWidget->takeTopLevelItem( selectedIndex );
  mLayerOrderTreeWidget->insertTopLevelItem( selectedIndex - 1, selectedItem );
  mLayerOrderTreeWidget->clearSelection();
  selectedItem->setSelected( true );
}

void QgsOWSSourceSelect::mLayerDownButton_clicked()
{
  QList<QTreeWidgetItem *> selectionList = mLayerOrderTreeWidget->selectedItems();
  if ( selectionList.empty() )
  {
    return;
  }
  int selectedIndex = mLayerOrderTreeWidget->indexOfTopLevelItem( selectionList[0] );
  if ( selectedIndex < 0 || selectedIndex > mLayerOrderTreeWidget->topLevelItemCount() - 2 )
  {
    return; //item not existing or already at bottom
  }

  QTreeWidgetItem *selectedItem = mLayerOrderTreeWidget->takeTopLevelItem( selectedIndex );
  mLayerOrderTreeWidget->insertTopLevelItem( selectedIndex + 1, selectedItem );
  mLayerOrderTreeWidget->clearSelection();
  selectedItem->setSelected( true );
}

QList<QgsOWSSourceSelect::SupportedFormat> QgsOWSSourceSelect::providerFormats()
{
  return QList<SupportedFormat>();
}

QStringList QgsOWSSourceSelect::selectedLayersFormats()
{
  return QStringList();
}

QStringList QgsOWSSourceSelect::selectedLayersCrses()
{
  return QStringList();
}

QStringList QgsOWSSourceSelect::selectedLayersTimes()
{
  return QStringList();
}

void QgsOWSSourceSelect::updateButtons()
{
}
