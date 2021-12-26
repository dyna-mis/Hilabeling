# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsFeatureIterator.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Matthias Kuhn'
__date__ = '18/09/2013'
__copyright__ = 'Copyright 2013, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

import os

from qgis.core import (QgsVectorLayer,
                       QgsFeatureRequest,
                       QgsFeature,
                       QgsField,
                       NULL,
                       QgsProject,
                       QgsVectorLayerJoinInfo,
                       QgsGeometry)
from qgis.testing import start_app, unittest
from qgis.PyQt.QtCore import QVariant

from utilities import unitTestDataPath
start_app()
TEST_DATA_DIR = unitTestDataPath()


class TestQgsFeatureIterator(unittest.TestCase):

    def __init__(self, methodName):
        """Run once on class initialization."""
        unittest.TestCase.__init__(self, methodName)

    def test_FilterExpression(self):
        # create point layer
        myShpFile = os.path.join(TEST_DATA_DIR, 'points.shp')
        pointLayer = QgsVectorLayer(myShpFile, 'Points', 'ogr')

        ids = [feat.id() for feat in pointLayer.getFeatures(QgsFeatureRequest().setFilterExpression('Staff > 3'))]
        expectedIds = [1, 5, 6, 7, 8]
        myMessage = '\nExpected: {0} features\nGot: {1} features'.format(repr(expectedIds), repr(ids))
        assert ids == expectedIds, myMessage

        pointLayer.startEditing()
        self.addFeatures(pointLayer)

        ids = [feat.id() for feat in pointLayer.getFeatures(QgsFeatureRequest().setFilterExpression('Staff > 3'))]
        expectedIds = [-2, 1, 5, 6, 7, 8]
        myMessage = '\nExpected: {0} features\nGot: {1} features'.format(repr(expectedIds), repr(ids))
        assert ids == expectedIds, myMessage

        pointLayer.rollBack()

        ids = [feat.id() for feat in pointLayer.getFeatures(QgsFeatureRequest().setFilterExpression('Staff > 3'))]
        expectedIds = [1, 5, 6, 7, 8]
        myMessage = '\nExpected: {0} features\nGot: {1} features'.format(repr(expectedIds), repr(ids))
        assert ids == expectedIds, myMessage

    def test_FilterExpressionWithAccents(self):
        myShpFile = os.path.join(TEST_DATA_DIR, 'france_parts.shp')
        layer = QgsVectorLayer(myShpFile, 'poly', 'ogr')

        layer.setProviderEncoding("ISO-8859-1")
        ids = [feat.id() for feat in layer.getFeatures(QgsFeatureRequest().setFilterExpression("TYPE_1 = 'Région'"))]
        expectedIds = [0, 1, 2, 3]
        myMessage = '\nExpected: {0} features\nGot: {1} features'.format(repr(expectedIds), repr(ids))
        assert ids == expectedIds, myMessage

        layer.setProviderEncoding("UTF-8")
        ids = [feat.id() for feat in layer.getFeatures(QgsFeatureRequest().setFilterExpression("TYPE_1 = 'Région'"))]
        expectedIds = []
        myMessage = '\nExpected: {0} features\nGot: {1} features'.format(repr(expectedIds), repr(ids))
        assert ids == expectedIds, myMessage

    def test_FilterFids(self):
        # create point layer
        myShpFile = os.path.join(TEST_DATA_DIR, 'points.shp')
        pointLayer = QgsVectorLayer(myShpFile, 'Points', 'ogr')

        ids = [feat.id() for feat in pointLayer.getFeatures(QgsFeatureRequest().setFilterFids([7, 8, 12, 30]))]
        expectedIds = [7, 8, 12]
        self.assertEqual(set(ids), set(expectedIds))

        pointLayer.startEditing()
        self.addFeatures(pointLayer)

        ids = [feat.id() for feat in pointLayer.getFeatures(QgsFeatureRequest().setFilterFids([-4, 7, 8, 12, 30]))]
        expectedIds = [-4, 7, 8, 12]
        self.assertEqual(set(ids), set(expectedIds))

        pointLayer.rollBack()

        ids = [feat.id() for feat in pointLayer.getFeatures(QgsFeatureRequest().setFilterFids([-2, 7, 8, 12, 30]))]
        expectedIds = [7, 8, 12]
        self.assertEqual(set(ids), set(expectedIds))

    def addFeatures(self, vl):
        feat = QgsFeature()
        fields = vl.fields()
        feat.setFields(fields)
        feat['Staff'] = 4
        vl.addFeature(feat)

        feat = QgsFeature()
        fields = vl.fields()
        feat.setFields(fields)
        feat['Staff'] = 2
        vl.addFeature(feat)

    def test_ExpressionFieldNested(self):
        myShpFile = os.path.join(TEST_DATA_DIR, 'points.shp')
        layer = QgsVectorLayer(myShpFile, 'Points', 'ogr')
        self.assertTrue(layer.isValid())

        idx = layer.addExpressionField('"Staff"*2', QgsField('exp1', QVariant.LongLong))  # NOQA
        idx = layer.addExpressionField('"exp1"-1', QgsField('exp2', QVariant.LongLong))  # NOQA

        fet = next(layer.getFeatures(QgsFeatureRequest().setSubsetOfAttributes(['exp2'], layer.fields())))
        self.assertEqual(fet['Class'], NULL)
        # nested virtual fields should make all these attributes be fetched
        self.assertEqual(fet['Staff'], 2)
        self.assertEqual(fet['exp2'], 3)
        self.assertEqual(fet['exp1'], 4)

    def test_ExpressionFieldNestedGeometry(self):
        myShpFile = os.path.join(TEST_DATA_DIR, 'points.shp')
        layer = QgsVectorLayer(myShpFile, 'Points', 'ogr')
        self.assertTrue(layer.isValid())

        idx = layer.addExpressionField('$x*2', QgsField('exp1', QVariant.LongLong))  # NOQA
        idx = layer.addExpressionField('"exp1"/1.5', QgsField('exp2', QVariant.LongLong))  # NOQA

        fet = next(layer.getFeatures(QgsFeatureRequest().setFlags(QgsFeatureRequest.NoGeometry).setSubsetOfAttributes(['exp2'], layer.fields())))
        # nested virtual fields should have made geometry be fetched
        self.assertEqual(fet['exp2'], -156)
        self.assertEqual(fet['exp1'], -234)

    def test_ExpressionFieldNestedCircular(self):
        """ test circular virtual field definitions """

        myShpFile = os.path.join(TEST_DATA_DIR, 'points.shp')
        layer = QgsVectorLayer(myShpFile, 'Points', 'ogr')
        self.assertTrue(layer.isValid())

        cnt = layer.fields().count()  # NOQA
        idx = layer.addExpressionField('"exp3"*2', QgsField('exp1', QVariant.LongLong))  # NOQA
        idx = layer.addExpressionField('"exp1"-1', QgsField('exp2', QVariant.LongLong))  # NOQA
        idx = layer.addExpressionField('"exp2"*3', QgsField('exp3', QVariant.LongLong))  # NOQA

        # really just testing that this doesn't hang/crash... there's no good result here!
        fet = next(layer.getFeatures(QgsFeatureRequest().setSubsetOfAttributes(['exp2'], layer.fields())))
        self.assertEqual(fet['Class'], NULL)

    def test_JoinUsingExpression(self):
        """ test joining a layer using a virtual field """
        joinLayer = QgsVectorLayer(
            "Point?field=x:string&field=y:integer&field=z:integer",
            "joinlayer", "memory")
        pr = joinLayer.dataProvider()
        f1 = QgsFeature()
        f1.setAttributes(["foo", 246, 321])
        f2 = QgsFeature()
        f2.setAttributes(["bar", 456, 654])
        self.assertTrue(pr.addFeatures([f1, f2]))

        layer = QgsVectorLayer("Point?field=fldtxt:string&field=fldint:integer",
                               "addfeat", "memory")
        pr = layer.dataProvider()
        f = QgsFeature()
        f.setAttributes(["test", 123])
        self.assertTrue(pr.addFeatures([f]))
        layer.addExpressionField('"fldint"*2', QgsField('exp1', QVariant.LongLong))

        QgsProject.instance().addMapLayers([layer, joinLayer])

        join = QgsVectorLayerJoinInfo()
        join.setTargetFieldName("exp1")
        join.setJoinLayer(joinLayer)
        join.setJoinFieldName("y")
        join.setUsingMemoryCache(True)
        layer.addJoin(join)

        f = QgsFeature()
        fi = layer.getFeatures()
        self.assertTrue(fi.nextFeature(f))
        attrs = f.attributes()
        self.assertEqual(attrs[0], "test")
        self.assertEqual(attrs[1], 123)
        self.assertEqual(attrs[2], "foo")
        self.assertEqual(attrs[3], 321)
        self.assertEqual(attrs[4], 246)
        self.assertFalse(fi.nextFeature(f))

        QgsProject.instance().removeMapLayers([layer.id(), joinLayer.id()])

    def test_JoinUsingExpression2(self):
        """ test joining a layer using a virtual field (the other way!) """
        joinLayer = QgsVectorLayer(
            "Point?field=x:string&field=y:integer&field=z:integer",
            "joinlayer", "memory")
        pr = joinLayer.dataProvider()
        f1 = QgsFeature()
        f1.setAttributes(["foo", 246, 321])
        f2 = QgsFeature()
        f2.setAttributes(["bar", 456, 654])
        self.assertTrue(pr.addFeatures([f1, f2]))
        joinLayer.addExpressionField('"y"/2', QgsField('exp1', QVariant.LongLong))

        layer = QgsVectorLayer("Point?field=fldtxt:string&field=fldint:integer",
                               "addfeat", "memory")
        pr = layer.dataProvider()
        f = QgsFeature()
        f.setAttributes(["test", 123])
        self.assertTrue(pr.addFeatures([f]))

        QgsProject.instance().addMapLayers([layer, joinLayer])

        join = QgsVectorLayerJoinInfo()
        join.setTargetFieldName("fldint")
        join.setJoinLayer(joinLayer)
        join.setJoinFieldName("exp1")
        join.setUsingMemoryCache(True)
        layer.addJoin(join)

        f = QgsFeature()
        fi = layer.getFeatures()
        self.assertTrue(fi.nextFeature(f))
        attrs = f.attributes()
        self.assertEqual(attrs[0], "test")
        self.assertEqual(attrs[1], 123)
        self.assertEqual(attrs[2], "foo")
        self.assertEqual(attrs[3], 246)
        self.assertEqual(attrs[4], 321)
        self.assertFalse(fi.nextFeature(f))

        QgsProject.instance().removeMapLayers([layer.id(), joinLayer.id()])

    def test_JoinUsingFeatureRequestExpression(self):
        """ test requesting features using a filter expression which requires joined columns """
        joinLayer = QgsVectorLayer(
            "Point?field=x:string&field=y:integer&field=z:integer",
            "joinlayer", "memory")
        pr = joinLayer.dataProvider()
        f1 = QgsFeature()
        f1.setAttributes(["foo", 123, 321])
        f2 = QgsFeature()
        f2.setAttributes(["bar", 124, 654])
        self.assertTrue(pr.addFeatures([f1, f2]))

        layer = QgsVectorLayer("Point?field=fldtxt:string&field=fldint:integer",
                               "addfeat", "memory")
        pr = layer.dataProvider()
        f1 = QgsFeature()
        f1.setAttributes(["test", 123])
        f2 = QgsFeature()
        f2.setAttributes(["test", 124])
        self.assertTrue(pr.addFeatures([f1, f2]))

        QgsProject.instance().addMapLayers([layer, joinLayer])

        join = QgsVectorLayerJoinInfo()
        join.setTargetFieldName("fldint")
        join.setJoinLayer(joinLayer)
        join.setJoinFieldName("y")
        join.setUsingMemoryCache(True)
        layer.addJoin(join)

        f = QgsFeature()
        fi = layer.getFeatures(QgsFeatureRequest().setFlags(QgsFeatureRequest.SubsetOfAttributes).setFilterExpression('joinlayer_z=654'))
        self.assertTrue(fi.nextFeature(f))
        self.assertEqual(f['fldint'], 124)
        self.assertEqual(f['joinlayer_z'], 654)

        QgsProject.instance().removeMapLayers([layer.id(), joinLayer.id()])

    def test_invalidGeometryFilter(self):
        layer = QgsVectorLayer(
            "Polygon?field=x:string",
            "joinlayer", "memory")

        # add some features, one has invalid geometry
        pr = layer.dataProvider()
        f1 = QgsFeature(1)
        f1.setAttributes(["a"])
        f1.setGeometry(QgsGeometry.fromWkt('Polygon((0 0, 1 0, 1 1, 0 1, 0 0))')) # valid
        f2 = QgsFeature(2)
        f2.setAttributes(["b"])
        f2.setGeometry(QgsGeometry.fromWkt('Polygon((0 0, 1 0, 0 1, 1 1, 0 0))'))  # invalid
        f3 = QgsFeature(3)
        f3.setAttributes(["c"])
        f3.setGeometry(QgsGeometry.fromWkt('Polygon((0 0, 1 0, 1 1, 0 1, 0 0))'))  # valid
        self.assertTrue(pr.addFeatures([f1, f2, f3]))

        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck))]
        self.assertEqual(res, ['a', 'b', 'c'])
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometrySkipInvalid))]
        self.assertEqual(res, ['a', 'c'])
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometryAbortOnInvalid))]
        self.assertEqual(res, ['a'])

        # with callback
        self.callback_feature_val = None

        def callback(feature):
            self.callback_feature_val = feature['x']

        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(
                   QgsFeatureRequest.GeometryAbortOnInvalid).setInvalidGeometryCallback(callback))]
        self.assertEqual(res, ['a'])
        self.assertEqual(self.callback_feature_val, 'b')
        # clear callback
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(
                   QgsFeatureRequest.GeometryAbortOnInvalid).setInvalidGeometryCallback(None))]
        self.assertEqual(res, ['a'])

        # check with filter fids
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setFilterFid(f2.id()).setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck))]
        self.assertEqual(res, ['b'])
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setFilterFid(f2.id()).setInvalidGeometryCheck(QgsFeatureRequest.GeometrySkipInvalid))]
        self.assertEqual(res, [])
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setFilterFid(f2.id()).setInvalidGeometryCheck(QgsFeatureRequest.GeometryAbortOnInvalid))]
        self.assertEqual(res, [])

        f4 = QgsFeature(4)
        f4.setAttributes(["d"])
        f4.setGeometry(QgsGeometry.fromWkt('Polygon((0 0, 1 0, 0 1, 1 1, 0 0))'))  # invalid

        # check with added features
        layer.startEditing()
        self.assertTrue(layer.addFeatures([f4]))
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck))]
        self.assertEqual(set(res), {'a', 'b', 'c', 'd'})
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometrySkipInvalid))]
        self.assertEqual(set(res), {'a', 'c'})
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometryAbortOnInvalid))]
        self.assertEqual(res, ['a'])

        # check with features with changed geometry
        layer.rollBack()
        layer.startEditing()
        layer.changeGeometry(2, QgsGeometry.fromWkt('Polygon((0 0, 1 0, 1 1, 0 1, 0 0))'))  # valid
        layer.changeGeometry(3, QgsGeometry.fromWkt('Polygon((0 0, 1 0, 0 1, 1 1, 0 0))'))# invalid
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck))]
        self.assertEqual(set(res), {'a', 'b', 'c'})
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometrySkipInvalid))]
        self.assertEqual(set(res), {'a', 'b'})
        res = [f['x'] for f in
               layer.getFeatures(QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometryAbortOnInvalid))]
        self.assertEqual(res, ['a', 'b'])
        layer.rollBack()


if __name__ == '__main__':
    unittest.main()
