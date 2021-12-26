/***************************************************************************
                              qgswmsrendercontext.h
                              ---------------------
  begin                : March 22, 2019
  copyright            : (C) 2019 by Paul Blottiere
  email                : paul.blottiere@oslandia.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QGSWMSRENDERCONTEXT_H
#define QGSWMSRENDERCONTEXT_H

#include "qgswmsparameters.h"
#include "qgsproject.h"
#include "qgsserverinterface.h"

namespace QgsWms
{

  /**
   * \ingroup server
   * \class QgsWmsRenderContext
   * \brief Rendering context for the WMS renderer
   * \since QGIS 3.8
   */
  class QgsWmsRenderContext
  {
    public:
      //! Available rendering options
      enum Flag
      {
        UseScaleDenominator    = 0x01,
        UseOpacity             = 0x02,
        UseFilter              = 0x04,
        UseSelection           = 0x08,
        AddHighlightLayers     = 0x10,
        UpdateExtent           = 0x20,
        SetAccessControl       = 0x40,
        AddQueryLayers         = 0x80,
        UseWfsLayersOnly       = 0x100,
        AddExternalLayers      = 0x200
      };
      Q_DECLARE_FLAGS( Flags, Flag )

      /**
       * Default constructor for QgsWmsRenderContext.
       */
      QgsWmsRenderContext() = default;

      /**
       * Constructor for QgsWmsRenderContext.
       * \param project The project to use for the rendering
       * \param interface The server interface
       */
      QgsWmsRenderContext( const QgsProject *project, QgsServerInterface *interface );

      /**
       * Sets WMS parameters.
       */
      void setParameters( const QgsWmsParameters &parameters );

      /**
       * Returns WMS parameters.
       */
      QgsWmsParameters parameters() const;

      /**
       * Returns settings of the server.
       */
      const QgsServerSettings &settings() const;

      /**
       * Returns the project.
       */
      const QgsProject *project() const;

      /**
       * Sets or unsets a rendering flag according to the \a on value.
       */
      void setFlag( Flag flag, bool on = true );

      /**
       * Returns the status of a rendering flag.
       * \param flag The flag to test
       * \returns true if the rendering option is activated, false otherwise
       */
      bool testFlag( Flag flag ) const;

      /**
       * Returns a list of all layers read from the project.
       */
      QList<QgsMapLayer *> layers() const;

      /**
       * Returns a list of all layers to actually render according to the
       * current configuration.
       */
      QList<QgsMapLayer *> layersToRender() const;

      /**
       * Returns a SLD document for a specific layer. An empty document is
       * returned if not available.
       */
      QDomElement sld( const QgsMapLayer &layer ) const;

      /**
       * Returns a style's name for a specific layer. An empty string is
       * returned if not available.
       */
      QString style( const QgsMapLayer &layer ) const;

      /**
       * Returns the scale denominator to use for rendering according to the
       * current configuration.
       */
      double scaleDenominator() const;

      /**
       * Sets a custom scale denominator. In this case, layers to render are
       * updated according to their scale visibility.
       */
      void setScaleDenominator( double scaleDenominator );

      /**
       * Returns true if the extent has to be updated before the rendering,
       * false otherwise.
       */
      bool updateExtent() const;

      /**
       * Returns WMS parameters for a specific layer. An empty instance is
       * returned if not available.
       */
      QgsWmsParametersLayer parameters( const QgsMapLayer &layer ) const;

      /**
       * Returns the image quality to use for rendering according to the
       * current configuration.
       */
      int imageQuality() const;

      /**
       * Returns the precision to use according to the current configuration.
       */
      int precision() const;

      /**
       * Returns the nickname (short name, id or name) of the layer according
       * to the current configuration.
       */
      QString layerNickname( const QgsMapLayer &layer ) const;

      /**
       * Returns the layer corresponding to the nickname, or a nullptr if not
       * found or if the layer do not need to be rendered.
       */
      QgsMapLayer *layer( const QString &nickname ) const;

      /**
       * Returns true if the layer has to be rendered, false otherwise.
       */
      bool isValidLayer( const QString &nickname ) const;

      /**
       * Returns true if \a name is a group.
       */
      bool isValidGroup( const QString &name ) const;

      /**
       * Returns default dots per mm according to the current configuration.
       */
      qreal dotsPerMm() const;

      /**
       * Returns a list of query layer names where group names are replaced by the names of their layer components.
       * \since QGIS 3.8
       */
      QStringList flattenedQueryLayers() const;

#ifdef HAVE_SERVER_PYTHON_PLUGINS

      /**
       * Returns the access control interface.
       */
      QgsAccessControl *accessControl() const;
#endif

      /**
       * Returns a map having layer group names as keys and a list of layers as values.
       * \since QGIS 3.8
       */
      QMap<QString, QList<QgsMapLayer *> > layerGroups() const;

    private:
      void initNicknameLayers();
      void initRestrictedLayers();
      void initLayerGroupsRecursive( const QgsLayerTreeGroup *group, const QString &groupName );

      void searchLayersToRender();
      void searchLayersToRenderSld();
      void searchLayersToRenderStyle();
      void removeUnwantedLayers();

      void checkLayerReadPermissions();

      bool layerScaleVisibility( const QString &name ) const;

      const QgsProject *mProject = nullptr;
      QgsServerInterface *mInterface = nullptr;
      QgsWmsParameters mParameters;
      Flags mFlags = nullptr;
      double mScaleDenominator = -1.0;

      // nickname of all layers defined within the project
      QMap<QString, QgsMapLayer *> mNicknameLayers;

      // map of layers to use for rendering
      QList<QgsMapLayer *> mLayersToRender;

      // list of layers which are not usable
      QStringList mRestrictedLayers;

      QMap<QString, QList<QgsMapLayer *> > mLayerGroups;

      QMap<QString, QDomElement> mSlds;
      QMap<QString, QString> mStyles;
  };
};

#endif
