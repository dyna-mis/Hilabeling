/***************************************************************************
    qgshighlight.h - widget to highlight geometries
     --------------------------------------
    Date                 : 02-Mar-2011
    Copyright            : (C) 2011 by Juergen E. Fischer, norBIT GmbH
    Email                : jef at norbit dot de
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef QGSHIGHLIGHT_H
#define QGSHIGHLIGHT_H

#include "qgsmapcanvasitem.h"
#include "qgsgeometry.h"
#include "qgssymbol.h"
#include <QBrush>
#include <QColor>
#include <QList>
#include <QPen>
#include <QPainter>
#include <QPainterPath>
#include "qgis_gui.h"

class QgsMapLayer;
class QgsVectorLayer;
class QgsSymbol;

/**
 * \ingroup gui
 * A class for highlight features on the map.
 *
 * The QgsHighlight class provides a transparent overlay canvas item
 * for highlighting features or geometries on a map canvas.
 *
 * \code{.py}
 *   color = QColor(Qt.red)
 *   highlight = QgsHighlight(mapCanvas, feature, layer)
 *   highlight.setColor(color)
 *   color.setAlpha(50)
 *   highlight.setFillColor(color)
 *   highlight.show()
 * \endcode
 */
class GUI_EXPORT QgsHighlight: public QObject, public QgsMapCanvasItem
{

    Q_OBJECT
    Q_PROPERTY( QColor color READ color WRITE setColor )
    Q_PROPERTY( QColor fillColor READ fillColor WRITE setFillColor )
    Q_PROPERTY( int width READ width WRITE setWidth )
    Q_PROPERTY( int buffer READ buffer WRITE setBuffer )

  public:

    /**
     * Constructor for QgsHighlight
     * \param mapCanvas associated map canvas
     * \param geom initial geometry of highlight
     * \param layer associated map layer
     */
    QgsHighlight( QgsMapCanvas *mapCanvas, const QgsGeometry &geom, QgsMapLayer *layer );

    /**
     * Constructor for highlighting TRUE feature shape using feature attributes
     * and renderer.
     * \param mapCanvas map canvas
     * \param feature
     * \param layer vector layer
     */
    QgsHighlight( QgsMapCanvas *mapCanvas, const QgsFeature &feature, QgsVectorLayer *layer );
    ~QgsHighlight() override;

    /**
     * Returns the line/stroke color
     * \since QGIS 3.4
     */
    QColor color( ) const { return  mColor; }

    /**
     * Set line/stroke to color, polygon fill to color with alpha = 63.
     *  This is legacy function, use setFillColor() after setColor() if different fill color is required. */
    void setColor( const QColor &color );

    /**
     * Returns the fill color
     * \since QGIS 3.4
     */
    QColor fillColor( ) const { return mFillColor; }

    /**
     * Fill color for the highlight.
     * Will be used for polygons and points.
     *
     * \since QGIS 2.4
     */
    void setFillColor( const QColor &fillColor );

    /**
     * Returns the stroke width
     * \since QGIS 3.4
     */
    int width( ) const { return mWidth; }

    /**
     * Set stroke width.
     *
     * \note Ignored in feature mode.
     */
    void setWidth( int width );

    /**
     * Returns the buffer
     * \since QGIS 3.4
     */
    double buffer( ) const { return mBuffer; }

    /**
     * Set line / stroke buffer in millimeters.
     *
     * \since QGIS 2.4
     */
    void setBuffer( double buffer ) { mBuffer = buffer; }

    /**
     * Set minimum line / stroke width in millimeters.
     *
     * \since QGIS 2.4
     */
    void setMinWidth( double width ) { mMinWidth = width; }

    /**
     * Returns the layer for which this highlight has been created.
     */
    QgsMapLayer *layer() const { return mLayer; }

    void updatePosition() override;

  protected:
    void paint( QPainter *p ) override;

    //! recalculates needed rectangle
    void updateRect();

  private:
    void init();
    void setSymbol( QgsSymbol *symbol, const QgsRenderContext &context, const QColor &color, const QColor &fillColor );
    double getSymbolWidth( const QgsRenderContext &context, double width, QgsUnitTypes::RenderUnit unit );
    //! Gets renderer for current color mode and colors. The renderer should be freed by caller.
    std::unique_ptr< QgsFeatureRenderer > createRenderer( QgsRenderContext &context, const QColor &color, const QColor &fillColor );
    void paintPoint( QPainter *p, const QgsPointXY &point );
    void paintLine( QPainter *p, QgsPolylineXY line );
    void paintPolygon( QPainter *p, const QgsPolygonXY &polygon );

    int mWidth = 1; // line / stroke width property
    QColor mColor; // line / stroke color property
    QColor mFillColor; // line / stroke fillColor property
    QBrush mBrush;
    QPen mPen;
    QgsGeometry *mGeometry = nullptr;
    QgsMapLayer *mLayer = nullptr;
    QgsFeature mFeature;
    double mBuffer = 0; // line / stroke buffer in pixels
    double mMinWidth = 0; // line / stroke minimum width in pixels
};

#endif
