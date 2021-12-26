# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsExpressionBuilderWidget

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '30/07/2017'
__copyright__ = 'Copyright 2017, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

from qgis.PyQt.QtCore import Qt
from qgis.testing import start_app, unittest
from qgis.gui import QgsExpressionBuilderWidget
from qgis.core import (QgsExpressionContext,
                       QgsExpressionContextScope,
                       QgsProject,
                       QgsVectorLayer,
                       QgsRelation,
                       QgsFeature,
                       QgsGeometry)
start_app()


def createReferencingLayer():
    layer = QgsVectorLayer("Point?field=fldtxt:string&field=foreignkey:integer",
                           "referencinglayer", "memory")
    pr = layer.dataProvider()
    f1 = QgsFeature()
    f1.setFields(layer.fields())
    f1.setAttributes(["test1", 123])
    f2 = QgsFeature()
    f2.setFields(layer.fields())
    f2.setAttributes(["test2", 123])
    f3 = QgsFeature()
    f3.setFields(layer.fields())
    f3.setAttributes(["foobar'bar", 124])
    assert pr.addFeatures([f1, f2, f3])
    return layer


def createReferencedLayer():
    layer = QgsVectorLayer(
        "Point?field=x:string&field=y:integer&field=z:integer",
        "referencedlayer", "memory")
    pr = layer.dataProvider()
    f1 = QgsFeature()
    f1.setFields(layer.fields())
    f1.setAttributes(["foo", 123, 321])
    f2 = QgsFeature()
    f2.setFields(layer.fields())
    f2.setAttributes(["bar", 456, 654])
    f3 = QgsFeature()
    f3.setFields(layer.fields())
    f3.setAttributes(["foobar'bar", 789, 554])
    assert pr.addFeatures([f1, f2, f3])
    return layer


class TestQgsExpressionBuilderWidget(unittest.TestCase):

    def setUp(self):
        self.referencedLayer = createReferencedLayer()
        self.referencingLayer = createReferencingLayer()
        QgsProject.instance().addMapLayers([self.referencedLayer, self.referencingLayer])

    def testFunctionPresent(self):
        """ check through widget model to ensure it is initially populated with functions """
        w = QgsExpressionBuilderWidget()
        m = w.model()
        # check that some standard expression functions are shown
        items = m.findItems('lower', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)
        items = m.findItems('upper', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)
        items = m.findItems('asdasdasda#$@#$', Qt.MatchRecursive)
        self.assertEqual(len(items), 0)

    def testVariables(self):
        """ check through widget model to ensure it is populated with variables """
        w = QgsExpressionBuilderWidget()
        m = w.model()

        s = QgsExpressionContextScope()
        s.setVariable('my_var1', 'x')
        s.setVariable('my_var2', 'y')
        c = QgsExpressionContext()
        c.appendScope(s)

        # check that variables are added when setting context
        w.setExpressionContext(c)
        items = m.findItems('my_var1', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)
        items = m.findItems('my_var2', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)
        items = m.findItems('not_my_var', Qt.MatchRecursive)
        self.assertEqual(len(items), 0)
        # double check that functions are still only there once
        items = m.findItems('lower', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)
        items = m.findItems('upper', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)

    def testLayers(self):
        """ check that layers are shown in widget model"""
        p = QgsProject.instance()
        layer = QgsVectorLayer("Point", "layer1", "memory")
        layer2 = QgsVectorLayer("Point", "layer2", "memory")
        p.addMapLayers([layer, layer2])

        w = QgsExpressionBuilderWidget()
        m = w.model()

        # check that layers are shown
        items = m.findItems('layer1', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)
        items = m.findItems('layer2', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)

        # change project
        p2 = QgsProject()
        layer3 = QgsVectorLayer("Point", "layer3", "memory")
        p2.addMapLayers([layer3])
        w.setProject(p2)
        m = w.model()
        items = m.findItems('layer1', Qt.MatchRecursive)
        self.assertEqual(len(items), 0)
        items = m.findItems('layer2', Qt.MatchRecursive)
        self.assertEqual(len(items), 0)
        items = m.findItems('layer3', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)

    def testRelations(self):
        """ check that layers are shown in widget model"""
        p = QgsProject.instance()

        # not valid, but doesn't matter for test....
        rel = QgsRelation()
        rel.setId('rel1')
        rel.setName('Relation Number One')
        rel.setReferencingLayer(self.referencingLayer.id())
        rel.setReferencedLayer(self.referencedLayer.id())
        rel.addFieldPair('foreignkey', 'y')

        rel2 = QgsRelation()
        rel2.setId('rel2')
        rel2.setName('Relation Number Two')
        rel2.setReferencingLayer(self.referencingLayer.id())
        rel2.setReferencedLayer(self.referencedLayer.id())
        rel2.addFieldPair('foreignkey', 'y')

        p.relationManager().addRelation(rel)
        p.relationManager().addRelation(rel2)

        w = QgsExpressionBuilderWidget()
        m = w.model()

        # check that relations are shown
        items = m.findItems('Relation Number One', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)
        items = m.findItems('Relation Number Two', Qt.MatchRecursive)
        self.assertEqual(len(items), 1)


if __name__ == '__main__':
    unittest.main()
