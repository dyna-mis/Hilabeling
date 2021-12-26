/***************************************************************************
                         qgisappstylesheet.cpp
                         ----------------------
    begin                : Jan 18, 2013
    copyright            : (C) 2013 by Larry Shaffer
    email                : larrys at dakotacarto dot com

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgisappstylesheet.h"
#include "qgsapplication.h"
#include "qgisapp.h"
#include "qgslogger.h"
#include "qgssettings.h"

#include <QFont>
#include <QStyle>

QgisAppStyleSheet::QgisAppStyleSheet( QObject *parent )
  : QObject( parent )
{
  setActiveValues();
}

QMap<QString, QVariant> QgisAppStyleSheet::defaultOptions()
{
  QMap<QString, QVariant> opts;

  // the following default values, before insertion in opts, can be
  // configured using the platforms and window servers defined in the
  // constructor to set reasonable non-Qt defaults for the app stylesheet
  QgsSettings settings;
  // handle move from old QgsSettings group (/) to new (/qgis/stylesheet)
  // NOTE: don't delete old QgsSettings keys, in case user is also running older QGIS
  QVariant oldFontPointSize = settings.value( QStringLiteral( "fontPointSize" ) );
  QVariant oldFontFamily = settings.value( QStringLiteral( "fontFamily" ) );

  settings.beginGroup( QStringLiteral( "qgis/stylesheet" ) );

  int fontSize = mDefaultFont.pointSize();
  if ( mAndroidOS )
  {
    // TODO: find a better default fontsize maybe using DPI detection or so (from Marco Bernasocchi commit)
    fontSize = 8;
  }
  if ( oldFontPointSize.isValid() && !settings.value( QStringLiteral( "fontPointSize" ) ).isValid() )
  {
    fontSize = oldFontPointSize.toInt();
  }
  QgsDebugMsg( QStringLiteral( "fontPointSize: %1" ).arg( fontSize ) );
  opts.insert( QStringLiteral( "fontPointSize" ), settings.value( QStringLiteral( "fontPointSize" ), QVariant( fontSize ) ) );

  QString fontFamily = mDefaultFont.family();
  if ( oldFontFamily.isValid() && !settings.value( QStringLiteral( "fontFamily" ) ).isValid() )
  {
    fontFamily = oldFontFamily.toString();
  }
  fontFamily = settings.value( QStringLiteral( "fontFamily" ), QVariant( fontFamily ) ).toString();
  // make sure family exists on system
  if ( fontFamily != mDefaultFont.family() )
  {
    QFont tempFont( fontFamily );
    if ( tempFont.family() != fontFamily )
    {
      // missing from system, drop back to default
      fontFamily = mDefaultFont.family();
    }
  }
  QgsDebugMsg( QStringLiteral( "fontFamily: %1" ).arg( fontFamily ) );
  opts.insert( QStringLiteral( "fontFamily" ), QVariant( fontFamily ) );

  opts.insert( QStringLiteral( "toolbarSpacing" ), settings.value( QStringLiteral( "toolbarSpacing" ), QString() ) );

  settings.endGroup(); // "qgis/stylesheet"

  opts.insert( QStringLiteral( "iconSize" ), settings.value( QStringLiteral( "/qgis/iconSize" ), QGIS_ICON_SIZE ) );

  return opts;
}

void QgisAppStyleSheet::buildStyleSheet( const QMap<QString, QVariant> &opts )
{
  QgsSettings settings;
  QString ss;

  // QgisApp-wide font
  QString fontSize = opts.value( QStringLiteral( "fontPointSize" ) ).toString();
  QgsDebugMsg( QStringLiteral( "fontPointSize: %1" ).arg( fontSize ) );
  if ( fontSize.isEmpty() ) { return; }

  QString fontFamily = opts.value( QStringLiteral( "fontFamily" ) ).toString();
  QgsDebugMsg( QStringLiteral( "fontFamily: %1" ).arg( fontFamily ) );
  if ( fontFamily.isEmpty() ) { return; }

  const QString defaultSize = QString::number( mDefaultFont.pointSize() );
  const QString defaultFamily = mDefaultFont.family();
  if ( fontSize != defaultSize || fontFamily != defaultFamily )
    ss += QStringLiteral( "* { font: %1pt \"%2\"} " ).arg( fontSize, fontFamily );

  // Fix for macOS Qt 5.9+, where close boxes do not show on document mode tab bar tabs
  // See: https://bugreports.qt.io/browse/QTBUG-61092
  //      https://bugreports.qt.io/browse/QTBUG-61742
  // Setting any stylesheet makes the default close button disappear.
  // Specifically setting a custom close button temporarily works around issue.
  // TODO: Remove when regression is fixed (Qt 5.9.3 or 5.10?); though hard to tell,
  //       since we are overriding the default close button image now.
  if ( mMacStyle )
  {
    ss += QLatin1String( "QTabBar::close-button{ image: url(:/images/themes/default/mIconCloseTab.svg); }" );
    ss += QLatin1String( "QTabBar::close-button:hover{ image: url(:/images/themes/default/mIconCloseTabHover.svg); }" );
  }

  ss += QLatin1String( "QGroupBox{ font-weight: 600; }" );

  QString themeName = settings.value( QStringLiteral( "UI/UITheme" ), "default" ).toString();
  if ( themeName == QStringLiteral( "default" ) || !QgsApplication::uiThemes().contains( themeName ) )
  {
    //sidebar style
    QString style = "QListWidget#mOptionsListWidget {"
                    "    background-color: rgb(69, 69, 69, 0);"
                    "    outline: 0;"
                    "}"
                    "QFrame#mOptionsListFrame {"
                    "    background-color: rgb(69, 69, 69, 220);"
                    "}"
                    "QListWidget#mOptionsListWidget::item {"
                    "    color: white;"
                    "    padding: 3px;"
                    "}"
                    "QListWidget#mOptionsListWidget::item::selected {"
                    "    color: black;"
                    "    background-color:palette(Window);"
                    "    padding-right: 0px;"
                    "}";

    QString toolbarSpacing = opts.value( QStringLiteral( "toolbarSpacing" ), QString() ).toString();
    if ( !toolbarSpacing.isEmpty() )
    {
      bool ok = false;
      int toolbarSpacingInt = toolbarSpacing.toInt( &ok );
      if ( ok )
      {
        style += QStringLiteral( "QToolBar > QToolButton { padding: %1px; } " ).arg( toolbarSpacingInt );
      }
    }
    ss += style;

    // Fix selection color on losing focus (Windows)
    const QPalette palette = qApp->palette();

    ss += QString( "QTableView {"
                   "selection-background-color: %1;"
                   "selection-color: %2;"
                   "}" )
          .arg( palette.highlight().color().name(),
                palette.highlightedText().color().name() );

    ss += QStringLiteral( "QgsPropertyOverrideButton { background: none; border: 1px solid rgba(0, 0, 0, 0%); } QgsPropertyOverrideButton:focus { border: 1px solid palette(highlight); }" );
#ifdef Q_OS_MACX
    ss += QStringLiteral( "QgsPropertyOverrideButton::menu-indicator { width: 5px; }" );
#endif
  }

  QgsDebugMsg( QStringLiteral( "Stylesheet built: %1" ).arg( ss ) );

  emit appStyleSheetChanged( ss );
}

void QgisAppStyleSheet::saveToSettings( const QMap<QString, QVariant> &opts )
{
  QgsSettings settings;
  settings.beginGroup( QStringLiteral( "qgis/stylesheet" ) );

  QMap<QString, QVariant>::const_iterator opt = opts.constBegin();
  while ( opt != opts.constEnd() )
  {
    settings.setValue( QString( opt.key() ), opt.value() );
    ++opt;
  }
  settings.endGroup(); // "qgis/stylesheet"
}

void QgisAppStyleSheet::setActiveValues()
{
  mStyle = qApp->style()->objectName(); // active style name (lowercase)
  QgsDebugMsg( QStringLiteral( "Style name: %1" ).arg( mStyle ) );

  mMacStyle = mStyle.contains( QLatin1String( "macintosh" ) ); // macintosh (aqua)
  mOxyStyle = mStyle.contains( QLatin1String( "oxygen" ) ); // oxygen

  mDefaultFont = qApp->font(); // save before it is changed in any way

  // platforms, specific
#ifdef Q_OS_WIN
  mWinOS = true;
#else
  mWinOS = false;
#endif
#ifdef ANDROID
  mAndroidOS = true;
#else
  mAndroidOS = false;
#endif

}
