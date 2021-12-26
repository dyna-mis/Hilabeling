/***************************************************************************
                            qgslayoutscalebarwidget.h
                            ---------------------------
    begin                : 11 June 2008
    copyright            : (C) 2008 by Marco Hugentobler
    email                : marco dot hugentobler at karto dot baug dot ethz dot ch
 ***************************************************************************/
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QGSLAYOUTSCALEBARWIDGET_H
#define QGSLAYOUTSCALEBARWIDGET_H

#include "ui_qgslayoutscalebarwidgetbase.h"
#include "qgslayoutitemwidget.h"

#include <QButtonGroup>

class QgsLayoutItemScaleBar;

/**
 * \ingroup app
 * A widget to define the properties of a QgsLayoutItemScaleBar.
 */
class QgsLayoutScaleBarWidget: public QgsLayoutItemBaseWidget, private Ui::QgsLayoutScaleBarWidgetBase
{
    Q_OBJECT

  public:
    explicit QgsLayoutScaleBarWidget( QgsLayoutItemScaleBar *scaleBar );

  protected:

    bool setNewItem( QgsLayoutItem *item ) override;

  public slots:

    void mHeightSpinBox_valueChanged( double d );
    void mLineWidthSpinBox_valueChanged( double d );
    void mSegmentSizeSpinBox_valueChanged( double d );
    void mSegmentsLeftSpinBox_valueChanged( int i );
    void mNumberOfSegmentsSpinBox_valueChanged( int i );
    void mUnitLabelLineEdit_textChanged( const QString &text );
    void mMapUnitsPerBarUnitSpinBox_valueChanged( double d );
    void mFillColorButton_colorChanged( const QColor &newColor );
    void mFillColor2Button_colorChanged( const QColor &newColor );
    void mStrokeColorButton_colorChanged( const QColor &newColor );
    void mStyleComboBox_currentIndexChanged( const QString &text );
    void mLabelBarSpaceSpinBox_valueChanged( double d );
    void mBoxSizeSpinBox_valueChanged( double d );
    void mAlignmentComboBox_currentIndexChanged( int index );
    void mUnitsComboBox_currentIndexChanged( int index );
    void mLineJoinStyleCombo_currentIndexChanged( int index );
    void mLineCapStyleCombo_currentIndexChanged( int index );
    void mMinWidthSpinBox_valueChanged( double d );
    void mMaxWidthSpinBox_valueChanged( double d );

  private slots:
    void setGuiElements();
    void segmentSizeRadioChanged( QAbstractButton *radio );
    void mapChanged( QgsLayoutItem *item );
    void textFormatChanged();

  private:
    QPointer< QgsLayoutItemScaleBar > mScalebar;
    QgsLayoutItemPropertiesWidget *mItemPropertiesWidget = nullptr;

    QButtonGroup mSegmentSizeRadioGroup;

    //! Enables/disables the signals of the input gui elements
    void blockMemberSignals( bool enable );

    //! Enables/disables controls based on scale bar style
    void toggleStyleSpecificControls( const QString &style );

    void connectUpdateSignal();
    void disconnectUpdateSignal();
};

#endif //QGSLAYOUTSCALEBARWIDGET_H
