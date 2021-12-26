/***************************************************************************
  qgsquickmapsettings.h
  --------------------------------------
  Date                 : 27.12.2014
  Copyright            : (C) 2014 by Matthias Kuhn
  Email                : matthias (at) opengis.ch
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QGSQUICKMAPSETTINGS_H
#define QGSQUICKMAPSETTINGS_H

#include <QObject>

#include "qgscoordinatetransformcontext.h"
#include "qgsmapsettings.h"
#include "qgsmapthemecollection.h"
#include "qgspoint.h"
#include "qgsrectangle.h"

#include "qgis_quick.h"

class QgsProject;

/**
 * \ingroup quick
 * The QgsQuickMapSettings class encapsulates QgsMapSettings class to offer
 * settings of configuration of map rendering via QML properties.
 *
 * On top of QgsMapSettings functionality, when QgsProject is attached,
 * it automatically loads its default settings from the project.
 * QgsProject should be attached before it is read.
 *
 * \note QML Type: MapSettings
 *
 * \sa QgsMapCanvas
 *
 * \since QGIS 3.2
 */
class QUICK_EXPORT QgsQuickMapSettings : public QObject
{
    Q_OBJECT

    /**
     * A project property should be used as a primary source of project all other components
     * in the application. QgsProject should be attached to QgsQuickMapSettings before
     * it is read (QgsProject::read)
     *
     * When project is read, map settings (CRS, extent, ...) are automatically set from its DOM.
     */
    Q_PROPERTY( QgsProject *project READ project WRITE setProject NOTIFY projectChanged )

    /**
     * Geographical coordinates of the rectangle that should be rendered.
     * The actual visible extent used for rendering could be slightly different
     * since the given extent may be expanded in order to fit the aspect ratio
     * of output size. Use QgsQuickMapSettings::visibleExtent to get the resulting extent.
     *
     * Automatically loaded from project on QgsProject::readProject
     */
    Q_PROPERTY( QgsRectangle extent READ extent WRITE setExtent NOTIFY extentChanged )
    //! \copydoc QgsMapSettings::visibleExtent()
    Q_PROPERTY( QgsRectangle visibleExtent READ visibleExtent NOTIFY visibleExtentChanged )
    //! \copydoc QgsMapSettings::mapUnitsPerPixel()
    Q_PROPERTY( double mapUnitsPerPixel READ mapUnitsPerPixel NOTIFY mapUnitsPerPixelChanged )

    /**
     * The rotation of the resulting map image, in degrees clockwise.
     * Map canvas rotation support is not implemented, 0 is used
     */
    Q_PROPERTY( double rotation READ rotation WRITE setRotation NOTIFY rotationChanged )

    /**
     * The size of the resulting map image
     *
     * Automatically loaded from project on QgsProject::readProject
     */
    Q_PROPERTY( QSize outputSize READ outputSize WRITE setOutputSize NOTIFY outputSizeChanged )

    /**
     * Output DPI used for conversion between real world units (e.g. mm) and pixels
     *
     * Automatically loaded from project on QgsProject::readProject
     */
    Q_PROPERTY( double outputDpi READ outputDpi WRITE setOutputDpi NOTIFY outputDpiChanged )

    /**
      * CRS of destination coordinate reference system.
      *
      * Automatically loaded from project on QgsProject::readProject
      */
    Q_PROPERTY( QgsCoordinateReferenceSystem destinationCrs READ destinationCrs WRITE setDestinationCrs NOTIFY destinationCrsChanged )

    /**
     * Set list of layers for map rendering. The layers must be registered in QgsProject.
     * The layers are stored in the reverse order of how they are rendered (layer with index 0 will be on top)
     *
     * \note Any non-spatial layers will be automatically stripped from the list (since they cannot be rendered!).
     *
     * Not loaded automatically from the associated project
     */
    Q_PROPERTY( QList<QgsMapLayer *> layers READ layers WRITE setLayers NOTIFY layersChanged )

  public:
    //! Create new map settings
    QgsQuickMapSettings( QObject *parent = nullptr );
    ~QgsQuickMapSettings() = default;

    //! Clone map settings
    QgsMapSettings mapSettings() const;

    //! \copydoc QgsMapSettings::extent()
    QgsRectangle extent() const;

    //! \copydoc QgsMapSettings::setExtent()
    void setExtent( const QgsRectangle &extent );

    //! \copydoc QgsQuickMapSettings::project
    void setProject( QgsProject *project );

    //! \copydoc QgsQuickMapSettings::project
    QgsProject *project() const;

    //! Move current map extent to have center point defined by \a center
    Q_INVOKABLE void setCenter( const QgsPoint &center );

    //! \copydoc QgsMapSettings::mapUnitsPerPixel()
    double mapUnitsPerPixel() const;

    //! \copydoc QgsMapSettings::visibleExtent()
    QgsRectangle visibleExtent() const;

    //! \copydoc QgsMapSettings::transformContext()
    Q_INVOKABLE QgsCoordinateTransformContext transformContext() const;

    /**
     * Convert a map coordinate to screen pixel coordinates
     *
     * \param point A coordinate in map coordinates
     *
     * \return A coordinate in pixel / screen space
     */
    Q_INVOKABLE QPointF coordinateToScreen( const QgsPoint &point ) const;


    /**
     * Convert a screen coordinate to a map coordinate
     *
     * \param point A coordinate in pixel / screen coordinates
     *
     * \return A coordinate in map coordinates
     */
    Q_INVOKABLE QgsPoint screenToCoordinate( const QPointF &point ) const;

    //! \copydoc QgsMapSettings::setTransformContext()
    void setTransformContext( const QgsCoordinateTransformContext &context );

    //! \copydoc QgsQuickMapSettings::rotation
    double rotation() const;

    //! \copydoc QgsQuickMapSettings::rotation
    void setRotation( double rotation );

    //! \copydoc QgsMapSettings::outputSize()
    QSize outputSize() const;

    //! \copydoc QgsMapSettings::setOutputSize()
    void setOutputSize( const QSize &outputSize );

    //! \copydoc QgsMapSettings::outputDpi()
    double outputDpi() const;

    //! \copydoc QgsMapSettings::setOutputDpi()
    void setOutputDpi( double outputDpi );

    //! \copydoc QgsMapSettings::destinationCrs()
    QgsCoordinateReferenceSystem destinationCrs() const;

    //! \copydoc QgsMapSettings::setDestinationCrs()
    void setDestinationCrs( const QgsCoordinateReferenceSystem &destinationCrs );

    //! \copydoc QgsMapSettings::layers()
    QList<QgsMapLayer *> layers() const;

    //! \copydoc QgsMapSettings::setLayers()
    void setLayers( const QList<QgsMapLayer *> &layers );

  signals:
    //! \copydoc QgsQuickMapSettings::project
    void projectChanged();

    //! \copydoc QgsQuickMapSettings::extent
    void extentChanged();

    //! \copydoc QgsQuickMapSettings::destinationCrs
    void destinationCrsChanged();

    //! \copydoc QgsQuickMapSettings::mapUnitsPerPixel
    void mapUnitsPerPixelChanged();

    //! \copydoc QgsQuickMapSettings::rotation
    void rotationChanged();

    //! \copydoc QgsQuickMapSettings::visibleExtent
    void visibleExtentChanged();

    //! \copydoc QgsQuickMapSettings::outputSize
    void outputSizeChanged();

    //! \copydoc QgsQuickMapSettings::outputDpi
    void outputDpiChanged();

    //! \copydoc QgsQuickMapSettings::layers
    void layersChanged();

  private slots:

    /**
     * Read map canvas settings stored in a QGIS project file
     *
     * \param doc parsed DOM of a QgsProject
     */
    void onReadProject( const QDomDocument &doc );

  private:
    QgsProject *mProject = nullptr;
    QgsMapSettings mMapSettings;

};

#endif // QGSQUICKMAPSETTINGS_H
