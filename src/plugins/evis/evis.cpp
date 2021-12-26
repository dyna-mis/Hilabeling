/*
** File: evis.cpp
** Author: Peter J. Ersts ( ersts at amnh.org )
** Creation Date: 2007-03-06
**
** Copyright ( c ) 2007, American Museum of Natural History. All rights reserved.
**
** This library/program is free software; you can redistribute it
** and/or modify it under the terms of the GNU Library General Public
** License as published by the Free Software Foundation; either
** version 2 of the License, or ( at your option ) any later version.
**
** This library/program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
** Library General Public License for more details.
**
** This work was made possible through a grant by the John D. and
** Catherine T. MacArthur Foundation. Additionally, this program was prepared by
** the American Museum of Natural History under award No. NA05SEC46391002
** from the National Oceanic and Atmospheric Administration, U.S. Department
** of Commerce.  The statements, findings, conclusions, and recommendations
** are those of the author( s ) and do not necessarily reflect the views of the
** National Oceanic and Atmospheric Administration or the Department of Commerce.
**
**/

//This file was created using the plugin generator distributed with QGIS evis.h
//is based on and a modification of the default plugin.h file which carried the
//following header
/***************************************************************************
  evis.cpp
  An event visualization plugin for QGIS

  -------------------
         begin                : [PluginDate]
         copyright            : [( C ) Your Name and Date]
         email                : [Your Email]

 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   ( at your option ) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "evis.h"

//
// QGIS Specific includes
//
#include "qgsapplication.h"
#include "qgsrasterlayer.h"
#include "qgisinterface.h"
#include "qgsmaplayer.h"
#include "qgsguiutils.h"

//the gui subclass
#include "evisdatabaseconnectiongui.h"
#include "evisgenericeventbrowsergui.h"
#include "eviseventidtool.h"

//
// Qt Related Includes
//
#include <QMessageBox>
#include <QToolBar>
#include <QMenuBar>
#include <QMenu>
#include <QLineEdit>
#include <QAction>
#include <QApplication>
#include <QCursor>

static const QString sName = QObject::tr( "eVis" );
static const QString sDescription = QObject::tr( "An event visualization tool - view images associated with vector features" );
static const QString sCategory = QObject::tr( "Database" );
static const QString sPluginVersion = QObject::tr( "Version 1.1.0" );
static const QgisPlugin::PluginType sPluginType = QgisPlugin::UI;
static const QString sIcon = QStringLiteral( ":/evis/eVisEventBrowser.png" );



eVis::eVis( QgisInterface *qgisInterface )
  : QgisPlugin( sName, sDescription, sCategory, sPluginVersion, sPluginType )
  , mQGisIface( qgisInterface )
{
  mIdTool = nullptr;
}

void eVis::initGui()
{
  delete mDatabaseConnectionActionPointer;
  delete mEventIdToolActionPointer;
  delete mEventBrowserActionPointer;

  // Create the action for tool
  mDatabaseConnectionActionPointer = new QAction( QIcon( ":/evis/eVisDatabaseConnection.png" ), tr( "eVis Database Connection" ), this );
  mDatabaseConnectionActionPointer->setObjectName( QStringLiteral( "mDatabaseConnectionActionPointer" ) );
  mEventIdToolActionPointer = new QAction( QIcon( ":/evis/eVisEventIdTool.png" ), tr( "eVis Event Id Tool" ), this );
  mEventIdToolActionPointer->setObjectName( QStringLiteral( "mEventIdToolActionPointer" ) );
  mEventBrowserActionPointer = new QAction( QIcon( ":/evis/eVisEventBrowser.png" ), tr( "eVis Event Browser" ), this );
  mEventBrowserActionPointer->setObjectName( QStringLiteral( "mEventBrowserActionPointer" ) );

  // Set the what's this text
  mDatabaseConnectionActionPointer->setWhatsThis( tr( "Create layer from a database query" ) );
  mEventIdToolActionPointer->setWhatsThis( tr( "Open an Event Browser and display the selected feature" ) );
  mEventBrowserActionPointer->setWhatsThis( tr( "Open an Event Browser to explore the current layer's features" ) );

  // Connect the action to the runmQGisIface->mapCanvas()
  connect( mDatabaseConnectionActionPointer, &QAction::triggered, this, &eVis::launchDatabaseConnection );
  connect( mEventIdToolActionPointer, &QAction::triggered, this, &eVis::launchEventIdTool );
  connect( mEventBrowserActionPointer, &QAction::triggered, this, &eVis::launchEventBrowser );


  // Add the icon to the toolbar
  mQGisIface->addDatabaseToolBarIcon( mDatabaseConnectionActionPointer );
  mQGisIface->addDatabaseToolBarIcon( mEventIdToolActionPointer );
  mQGisIface->addDatabaseToolBarIcon( mEventBrowserActionPointer );

  mQGisIface->addPluginToDatabaseMenu( QStringLiteral( "&eVis" ), mDatabaseConnectionActionPointer );
  mQGisIface->addPluginToDatabaseMenu( QStringLiteral( "&eVis" ), mEventIdToolActionPointer );
  mQGisIface->addPluginToDatabaseMenu( QStringLiteral( "&eVis" ), mEventBrowserActionPointer );

  mEventIdToolActionPointer->setCheckable( true );
}

//method defined in interface
void eVis::help()
{
  //implement me!
}

void eVis::launchDatabaseConnection()
{
  eVisDatabaseConnectionGui *myPluginGui = new eVisDatabaseConnectionGui( &mTemporaryFileList, mQGisIface->mainWindow(), QgsGuiUtils::ModalDialogFlags );
  myPluginGui->setAttribute( Qt::WA_DeleteOnClose );

  connect( myPluginGui, &eVisDatabaseConnectionGui::drawVectorLayer, this, &eVis::drawVectorLayer );
  myPluginGui->show();
}

void eVis::launchEventIdTool()
{
  if ( !mIdTool )
  {
    mIdTool = new eVisEventIdTool( mQGisIface->mapCanvas() );
    mIdTool->setAction( mEventIdToolActionPointer );
  }
  else
  {
    mQGisIface->mapCanvas()->setMapTool( mIdTool );
  }
}

void eVis::launchEventBrowser()
{
  eVisGenericEventBrowserGui *myPluginGui = new eVisGenericEventBrowserGui( mQGisIface->mainWindow(), mQGisIface, QgsGuiUtils::ModalDialogFlags );
  myPluginGui->setAttribute( Qt::WA_DeleteOnClose );
}

void eVis::unload()
{
  // remove the GUI
  mQGisIface->removePluginDatabaseMenu( QStringLiteral( "&eVis" ), mDatabaseConnectionActionPointer );
  mQGisIface->removeDatabaseToolBarIcon( mDatabaseConnectionActionPointer );
  delete mDatabaseConnectionActionPointer;

  mQGisIface->removePluginDatabaseMenu( QStringLiteral( "&eVis" ), mEventIdToolActionPointer );
  mQGisIface->removeDatabaseToolBarIcon( mEventIdToolActionPointer );
  delete mEventIdToolActionPointer;

  mQGisIface->removePluginDatabaseMenu( QStringLiteral( "&eVis" ), mEventBrowserActionPointer );
  mQGisIface->removeDatabaseToolBarIcon( mEventBrowserActionPointer );
  delete mEventBrowserActionPointer;

  while ( !mTemporaryFileList.isEmpty() )
  {
    delete ( mTemporaryFileList.takeLast() );
  }

  delete mIdTool;
}

void eVis::drawVectorLayer( const QString &pathNameQString, const QString &baseNameQString, const QString &providerQString )
{
  mQGisIface->addVectorLayer( pathNameQString, baseNameQString, providerQString );
}


//////////////////////////////////////////////////////////////////////////
//
//
//  THE FOLLOWING CODE IS AUTOGENERATED BY THE PLUGIN BUILDER SCRIPT
//    YOU WOULD NORMALLY NOT NEED TO MODIFY THIS, AND YOUR PLUGIN
//      MAY NOT WORK PROPERLY IF YOU MODIFY THIS INCORRECTLY
//
//
//////////////////////////////////////////////////////////////////////////

/**
 * Required extern functions needed  for every plugin
 * These functions can be called prior to creating an instance
 * of the plugin class
 */
// Class factory to return a new instance of the plugin class
QGISEXTERN QgisPlugin *classFactory( QgisInterface *qgisInterfacePointer )
{
  return new eVis( qgisInterfacePointer );
}
// Return the name of the plugin - note that we do not user class members as
// the class may not yet be insantiated when this method is called.
QGISEXTERN QString name()
{
  return sName;
}

// Return the description
QGISEXTERN QString description()
{
  return sDescription;
}

// Return the category
QGISEXTERN QString category()
{
  return sCategory;
}

// Return the type ( either UI or MapLayer plugin )
QGISEXTERN int type()
{
  return sPluginType;
}

// Return the icon
QGISEXTERN QString icon()
{
  return sIcon;
}

// Return the version number for the plugin
QGISEXTERN QString version()
{
  return sPluginVersion;
}

// Delete ourself
QGISEXTERN void unload( QgisPlugin *pluginPointer )
{
  delete pluginPointer;
}
