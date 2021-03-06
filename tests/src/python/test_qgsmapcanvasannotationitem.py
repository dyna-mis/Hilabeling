# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsMapCanvasAnnotationItem.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '24/1/2017'
__copyright__ = 'Copyright 2017, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

from qgis.core import (QgsTextAnnotation,
                       QgsCoordinateReferenceSystem,
                       QgsRectangle,
                       QgsPointXY,
                       QgsVectorLayer,
                       QgsFeature,
                       QgsGeometry,
                       QgsFillSymbol)
from qgis.gui import QgsMapCanvas, QgsMapCanvasAnnotationItem

from qgis.PyQt.QtCore import QPointF, QSizeF

from qgis.testing import start_app, unittest
from utilities import unitTestDataPath

start_app()
TEST_DATA_DIR = unitTestDataPath()


class TestQgsMapCanvasAnnotationItem(unittest.TestCase):

    def testPosition(self):
        """ test that map canvas annotation item syncs position correctly """
        a = QgsTextAnnotation()
        a.setFrameSize(QSizeF(300, 200))
        a.setFrameOffsetFromReferencePoint(QPointF(40, 50))
        a.setMapPosition(QgsPointXY(12, 34))
        a.setMapPositionCrs(QgsCoordinateReferenceSystem(4326))

        canvas = QgsMapCanvas()
        canvas.setDestinationCrs(QgsCoordinateReferenceSystem(4326))
        canvas.setFrameStyle(0)
        canvas.resize(600, 400)
        canvas.show()

        canvas.setExtent(QgsRectangle(10, 30, 20, 35))

        i = QgsMapCanvasAnnotationItem(a, canvas)
        self.assertEqual(canvas.width(), 600)
        self.assertEqual(canvas.height(), 400)

        # test that correct initial position is set
        self.assertAlmostEqual(i.pos().x(), 120, 1)
        self.assertAlmostEqual(i.pos().y(), 110, 1)

        # shift annotation map position, check that item is moved
        a.setMapPosition(QgsPointXY(14, 32))
        self.assertAlmostEqual(i.pos().x(), 240, 1)
        self.assertAlmostEqual(i.pos().y(), 230, 1)

        # check relative position
        a.setHasFixedMapPosition(False)
        a.setRelativePosition(QPointF(0.8, 0.4))
        self.assertAlmostEqual(i.pos().x(), 480, 1)
        self.assertAlmostEqual(i.pos().y(), 160, 1)

        # flicking between relative and fixed position
        a.setHasFixedMapPosition(True)
        self.assertAlmostEqual(i.pos().x(), 240, 1)
        self.assertAlmostEqual(i.pos().y(), 230, 1)
        a.setHasFixedMapPosition(False)
        self.assertAlmostEqual(i.pos().x(), 480, 1)
        self.assertAlmostEqual(i.pos().y(), 160, 1)

    def testSize(self):
        """ test that map canvas annotation item size is correct """
        a = QgsTextAnnotation()
        a.setFrameSize(QSizeF(300, 200))
        a.setHasFixedMapPosition(False)
        a.setFillSymbol(QgsFillSymbol.createSimple({'color': 'blue', 'width_border': '0'}))

        canvas = QgsMapCanvas()
        canvas.setDestinationCrs(QgsCoordinateReferenceSystem(4326))
        canvas.setFrameStyle(0)
        canvas.resize(600, 400)
        canvas.show()

        canvas.setExtent(QgsRectangle(10, 30, 20, 35))

        i = QgsMapCanvasAnnotationItem(a, canvas)
        self.assertAlmostEqual(i.boundingRect().width(), 300, 1)
        self.assertAlmostEqual(i.boundingRect().height(), 200, 1)

        a.setHasFixedMapPosition(True)
        a.setFrameOffsetFromReferencePoint(QPointF(0, 0))
        self.assertAlmostEqual(i.boundingRect().width(), 300, -1)
        self.assertAlmostEqual(i.boundingRect().height(), 200, -1)
        a.setFrameOffsetFromReferencePoint(QPointF(10, 20))
        self.assertAlmostEqual(i.boundingRect().width(), 310, -1)
        self.assertAlmostEqual(i.boundingRect().height(), 220, -1)

    def testVisibility(self):
        """ test that map canvas annotation item visibility follows layer"""
        a = QgsTextAnnotation()
        canvas = QgsMapCanvas()
        i = QgsMapCanvasAnnotationItem(a, canvas)

        self.assertTrue(i.isVisible())

        layer = QgsVectorLayer("Point?crs=EPSG:3111&field=fldtxt:string",
                               'test', "memory")
        a.setMapLayer(layer)
        self.assertFalse(i.isVisible())
        canvas.setLayers([layer])
        self.assertTrue(i.isVisible())

    def testSettingFeature(self):
        """ test that feature is set when item moves """
        a = QgsTextAnnotation()
        a.setFrameSize(QSizeF(300, 200))
        a.setFrameOffsetFromReferencePoint(QPointF(40, 50))
        a.setHasFixedMapPosition(True)
        a.setMapPosition(QgsPointXY(12, 34))
        a.setMapPositionCrs(QgsCoordinateReferenceSystem(4326))

        canvas = QgsMapCanvas()
        canvas.setDestinationCrs(QgsCoordinateReferenceSystem(4326))
        canvas.setFrameStyle(0)
        canvas.resize(600, 400)

        canvas.setExtent(QgsRectangle(10, 30, 20, 35))

        i = QgsMapCanvasAnnotationItem(a, canvas)  # NOQA

        layer = QgsVectorLayer("Point?crs=EPSG:4326&field=station:string&field=suburb:string",
                               'test', "memory")
        canvas.setLayers([layer])
        f = QgsFeature(layer.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(14, 31)))
        f.setValid(True)
        f.setAttributes(['hurstbridge', 'somewhere'])
        self.assertTrue(layer.dataProvider().addFeatures([f]))
        a.setMapLayer(layer)
        self.assertFalse(a.associatedFeature().isValid())

        a.setMapPosition(QgsPointXY(14, 31))
        self.assertTrue(a.associatedFeature().isValid())
        self.assertEqual(a.associatedFeature().attributes()[0], 'hurstbridge')

        a.setMapPosition(QgsPointXY(17, 31))
        self.assertFalse(a.associatedFeature().isValid())


if __name__ == '__main__':
    unittest.main()
