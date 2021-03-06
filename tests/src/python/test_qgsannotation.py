# -*- coding: utf-8 -*-
"""QGIS Unit tests for annotations.

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
                       QgsSvgAnnotation,
                       QgsHtmlAnnotation,
                       QgsMapSettings,
                       QgsRenderContext,
                       QgsCoordinateReferenceSystem,
                       QgsRectangle,
                       QgsMultiRenderChecker,
                       QgsRenderChecker,
                       QgsVectorLayer,
                       QgsFeature,
                       QgsMargins,
                       QgsFillSymbol)
from qgis.PyQt.QtCore import (QDir,
                              QPointF,
                              QSizeF)
from qgis.PyQt.QtGui import (QColor,
                             QPainter,
                             QImage,
                             QTextDocument)

from qgis.testing import start_app, unittest
from utilities import unitTestDataPath

start_app()
TEST_DATA_DIR = unitTestDataPath()


class TestQgsAnnotation(unittest.TestCase):

    def setUp(self):
        self.report = "<h1>Python QgsAnnotation Tests</h1>\n"

    def tearDown(self):
        report_file_path = "%s/qgistest.html" % QDir.tempPath()
        with open(report_file_path, 'a') as report_file:
            report_file.write(self.report)

    def testTextAnnotation(self):
        """ test rendering a text annotation"""
        a = QgsTextAnnotation()
        a.fillSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.markerSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.setFrameSize(QSizeF(300, 200))
        a.setFrameOffsetFromReferencePoint(QPointF(40, 50))
        doc = QTextDocument()
        doc.setHtml('<p style="font-family: arial; font-weight: bold; font-size: 40px;">test annotation</p>')
        a.setDocument(doc)
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('text_annotation', 'text_annotation', im))

        # check clone
        clone = a.clone()
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('text_annotation', 'text_annotation', im))

    def testSvgAnnotation(self):
        """ test rendering a svg annotation"""
        a = QgsSvgAnnotation()
        a.fillSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.markerSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.setFrameSize(QSizeF(300, 200))
        a.setFrameOffsetFromReferencePoint(QPointF(40, 50))
        svg = TEST_DATA_DIR + "/sample_svg.svg"
        a.setFilePath(svg)
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('svg_annotation', 'svg_annotation', im))

        # check clone
        clone = a.clone()
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('svg_annotation', 'svg_annotation', im))

    def testHtmlAnnotation(self):
        """ test rendering a html annotation"""
        a = QgsHtmlAnnotation()
        a.fillSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.markerSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.setFrameSize(QSizeF(400, 250))
        a.setFrameOffsetFromReferencePoint(QPointF(70, 90))
        html = TEST_DATA_DIR + "/test_html.html"
        a.setSourceFile(html)
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('html_annotation', 'html_annotation', im))

        # check clone
        clone = a.clone()
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('html_annotation', 'html_annotation', im))

    def testHtmlAnnotationWithFeature(self):
        """ test rendering a html annotation with a feature"""
        layer = QgsVectorLayer("Point?crs=EPSG:3111&field=station:string&field=suburb:string",
                               'test', "memory")

        a = QgsHtmlAnnotation()
        a.fillSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.markerSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.setFrameSize(QSizeF(400, 250))
        a.setFrameOffsetFromReferencePoint(QPointF(70, 90))
        a.setMapLayer(layer)
        html = TEST_DATA_DIR + "/test_html_feature.html"
        a.setSourceFile(html)
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('html_nofeature', 'html_nofeature', im))
        f = QgsFeature(layer.fields())
        f.setValid(True)
        f.setAttributes(['hurstbridge', 'somewhere'])
        a.setAssociatedFeature(f)
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('html_feature', 'html_feature', im))

    def testRelativePosition(self):
        """ test rendering an annotation without map point"""
        a = QgsHtmlAnnotation()
        a.fillSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.setFrameSize(QSizeF(400, 250))
        a.setHasFixedMapPosition(False)
        html = TEST_DATA_DIR + "/test_html.html"
        a.setSourceFile(html)
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('relative_style', 'relative_style', im))

    def testMargins(self):
        """ test rendering an annotation with margins"""
        a = QgsHtmlAnnotation()
        a.fillSymbol().symbolLayer(0).setStrokeColor(QColor(0, 0, 0))
        a.setFrameSize(QSizeF(400, 250))
        a.setHasFixedMapPosition(False)
        a.setContentsMargin(QgsMargins(15, 10, 30, 20))
        html = TEST_DATA_DIR + "/test_html.html"
        a.setSourceFile(html)
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('annotation_margins', 'annotation_margins', im))

    def testFillSymbol(self):
        """ test rendering an annotation with fill symbol"""
        a = QgsTextAnnotation()
        a.setFrameSize(QSizeF(400, 250))
        a.setHasFixedMapPosition(False)
        a.setFillSymbol(QgsFillSymbol.createSimple({'color': 'blue', 'width_border': '5', 'outline_color': 'black'}))
        im = self.renderAnnotation(a, QPointF(20, 30))
        self.assertTrue(self.imageCheck('annotation_fillstyle', 'annotation_fillstyle', im))

    def renderAnnotation(self, annotation, offset):
        image = QImage(600, 400, QImage.Format_RGB32)
        image.fill(QColor(0, 0, 0, 0))
        QgsRenderChecker.drawBackground(image)

        painter = QPainter()
        ms = QgsMapSettings()
        ms.setDestinationCrs(QgsCoordinateReferenceSystem(4326))
        extent = QgsRectangle(0, 5, 40, 30)

        ms.setExtent(extent)
        ms.setOutputSize(image.size())
        context = QgsRenderContext.fromMapSettings(ms)
        context.setPainter(painter)
        context.setScaleFactor(96 / 25.4)  # 96 DPI

        painter.begin(image)
        painter.translate(offset.x(), offset.y())
        annotation.render(context)
        painter.end()
        return image

    def imageCheck(self, name, reference_image, image):
        self.report += "<h2>Render {}</h2>\n".format(name)
        temp_dir = QDir.tempPath() + '/'
        file_name = temp_dir + 'annotation_' + name + ".png"
        image.save(file_name, "PNG")
        checker = QgsMultiRenderChecker()
        checker.setControlPathPrefix("annotations")
        checker.setControlName("expected_" + reference_image)
        checker.setRenderedImage(file_name)
        checker.setColorTolerance(2)
        result = checker.runTest(name, 20)
        self.report += checker.report()
        print((self.report))
        return result


if __name__ == '__main__':
    unittest.main()
