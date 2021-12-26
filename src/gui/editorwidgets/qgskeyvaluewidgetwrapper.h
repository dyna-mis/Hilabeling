/***************************************************************************
    qgskeyvaluewidgetwrapper.h
     --------------------------------------
    Date                 : 08.2016
    Copyright            : (C) 2016 Patrick Valsecchi
    Email                : patrick.valsecchi@camptocamp.com
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QGSKEYVALUEWIDGETWRAPPER_H
#define QGSKEYVALUEWIDGETWRAPPER_H

#include "qgseditorwidgetwrapper.h"
#include "qgis_gui.h"

SIP_NO_FILE

class QgsKeyValueWidget;

/**
 * \ingroup gui
 * Wraps a key/value widget.
 * \note not available in Python bindings
 * \since QGIS 3.0
 */
class GUI_EXPORT QgsKeyValueWidgetWrapper : public QgsEditorWidgetWrapper
{
    Q_OBJECT
  public:

    /**
     * Constructor for QgsKeyValueWidgetWrapper.
     *
     * The \a layer and \a fieldIdx arguments specify the vector layer field associated with the wrapper.
     *
     * The \a editor argument indicates the editor widget to use with the wrapper. This can be NULLPTR if a
     * new widget should be autogenerated.
     *
     * A \a parent widget for this widget wrapper and the created widget can also be specified.
     */
    explicit QgsKeyValueWidgetWrapper( QgsVectorLayer *layer, int fieldIdx, QWidget *editor = nullptr, QWidget *parent = nullptr );

    // QgsEditorWidgetWrapper interface
  public:
    QVariant value() const override;
    void showIndeterminateState() override;

  protected:
    QWidget *createWidget( QWidget *parent ) override;
    void initWidget( QWidget *editor ) override;
    bool valid() const override;

  public slots:
    void setValue( const QVariant &value ) override;

  private slots:
    void onValueChanged();

  private:
    void updateConstraintWidgetStatus() override;

    QgsKeyValueWidget *mWidget = nullptr;
};

#endif // QGSKEYVALUEWIDGETWRAPPER_H
