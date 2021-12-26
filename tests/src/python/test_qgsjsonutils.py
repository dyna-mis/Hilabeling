# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsJsonUtils.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '3/05/2016'
__copyright__ = 'Copyright 2016, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

from qgis.testing import unittest, start_app
from qgis.core import (QgsJsonUtils,
                       QgsJsonExporter,
                       QgsCoordinateReferenceSystem,
                       QgsProject,
                       QgsFeature,
                       QgsField,
                       QgsFields,
                       QgsWkbTypes,
                       QgsGeometry,
                       QgsPoint,
                       QgsLineString,
                       NULL,
                       QgsVectorLayer,
                       QgsRelation,
                       QgsEditorWidgetSetup
                       )
from qgis.PyQt.QtCore import QVariant, QTextCodec, QLocale

start_app()
codec = QTextCodec.codecForName("System")


class TestQgsJsonUtils(unittest.TestCase):

    def testStringToFeatureList(self):
        """Test converting json string to features"""

        fields = QgsFields()
        fields.append(QgsField("name", QVariant.String))

        # empty string
        features = QgsJsonUtils.stringToFeatureList("", fields, codec)
        self.assertEqual(features, [])

        # bad string
        features = QgsJsonUtils.stringToFeatureList("asdasdas", fields, codec)
        self.assertEqual(features, [])

        # geojson string with 1 feature
        features = QgsJsonUtils.stringToFeatureList('{\n"type": "Feature","geometry": {"type": "Point","coordinates": [125, 10]},"properties": {"name": "Dinagat Islands"}}', fields, codec)
        self.assertEqual(len(features), 1)
        self.assertFalse(features[0].geometry().isNull())
        self.assertEqual(features[0].geometry().wkbType(), QgsWkbTypes.Point)
        point = features[0].geometry().constGet()
        self.assertEqual(point.x(), 125.0)
        self.assertEqual(point.y(), 10.0)
        self.assertEqual(features[0]['name'], "Dinagat Islands")

        # geojson string with 2 features
        features = QgsJsonUtils.stringToFeatureList('{ "type": "FeatureCollection","features":[{\n"type": "Feature","geometry": {"type": "Point","coordinates": [125, 10]},"properties": {"name": "Dinagat Islands"}}, {\n"type": "Feature","geometry": {"type": "Point","coordinates": [110, 20]},"properties": {"name": "Henry Gale Island"}}]}', fields, codec)
        self.assertEqual(len(features), 2)
        self.assertFalse(features[0].geometry().isNull())
        self.assertEqual(features[0].geometry().wkbType(), QgsWkbTypes.Point)
        point = features[0].geometry().constGet()
        self.assertEqual(point.x(), 125.0)
        self.assertEqual(point.y(), 10.0)
        self.assertEqual(features[0]['name'], "Dinagat Islands")
        self.assertFalse(features[1].geometry().isNull())
        self.assertEqual(features[1].geometry().wkbType(), QgsWkbTypes.Point)
        point = features[1].geometry().constGet()
        self.assertEqual(point.x(), 110.0)
        self.assertEqual(point.y(), 20.0)
        self.assertEqual(features[1]['name'], "Henry Gale Island")

    def testStringToFields(self):
        """test retrieving fields from GeoJSON strings"""

        # empty string
        fields = QgsJsonUtils.stringToFields("", codec)
        self.assertEqual(fields.count(), 0)

        # bad string
        fields = QgsJsonUtils.stringToFields("asdasdas", codec)
        self.assertEqual(fields.count(), 0)

        # geojson string
        fields = QgsJsonUtils.stringToFields('{\n"type": "Feature","geometry": {"type": "Point","coordinates": [125, 10]},"properties": {"name": "Dinagat Islands","height":5.5}}', codec)
        self.assertEqual(fields.count(), 2)
        self.assertEqual(fields[0].name(), "name")
        self.assertEqual(fields[0].type(), QVariant.String)
        self.assertEqual(fields[1].name(), "height")
        self.assertEqual(fields[1].type(), QVariant.Double)

    def testEncodeValue(self):
        """ test encoding various values for use in GeoJSON strings """
        self.assertEqual(QgsJsonUtils.encodeValue(NULL), 'null')
        self.assertEqual(QgsJsonUtils.encodeValue(5), '5')
        self.assertEqual(QgsJsonUtils.encodeValue(5.9), '5.9')
        self.assertEqual(QgsJsonUtils.encodeValue(5999999999), '5999999999')
        self.assertEqual(QgsJsonUtils.encodeValue('string'), '"string"')
        self.assertEqual(QgsJsonUtils.encodeValue('str\ning'), '"str\\ning"')
        self.assertEqual(QgsJsonUtils.encodeValue('str\ring'), '"str\\ring"')
        self.assertEqual(QgsJsonUtils.encodeValue('str\bing'), '"str\\bing"')
        self.assertEqual(QgsJsonUtils.encodeValue('str\ting'), '"str\\ting"')
        self.assertEqual(QgsJsonUtils.encodeValue('str\\ing'), '"str\\\\ing"')
        self.assertEqual(QgsJsonUtils.encodeValue('str\\ning'), '"str\\\\ning"')
        self.assertEqual(QgsJsonUtils.encodeValue('str\n\\\\ing'), '"str\\n\\\\\\\\ing"')
        self.assertEqual(QgsJsonUtils.encodeValue('str/ing'), '"str\\/ing"')
        self.assertEqual(QgsJsonUtils.encodeValue([5, 6]), '[5,6]')
        self.assertEqual(QgsJsonUtils.encodeValue(['a', 'b', 'c']), '["a","b","c"]')
        self.assertEqual(QgsJsonUtils.encodeValue(['a', 3, 'c']), '["a",3,"c"]')
        self.assertEqual(QgsJsonUtils.encodeValue(['a', 'c\nd']), '["a","c\\nd"]')
        # handle differences due to Qt5 version, where compact output now lacks \n
        enc_str = QgsJsonUtils.encodeValue({'key': 'value', 'key2': 5})
        self.assertTrue(enc_str == '{"key":"value",\n"key2":5}' or enc_str == '{"key":"value","key2":5}')
        enc_str = QgsJsonUtils.encodeValue({'key': [1, 2, 3], 'key2': {'nested': 'nested\\result'}})
        self.assertTrue(enc_str == '{"key":[1,2,3],\n"key2":{"nested":"nested\\\\result"}}' or enc_str == '{"key":[1,2,3],"key2":{"nested":"nested\\\\result"}}')

    def testExportAttributes(self):
        """ test exporting feature's attributes to JSON object """
        fields = QgsFields()

        # test empty attributes
        feature = QgsFeature(fields, 5)
        expected = "{}"
        self.assertEqual(QgsJsonUtils.exportAttributes(feature), expected)

        # test feature with attributes
        fields.append(QgsField("name", QVariant.String))
        fields.append(QgsField("cost", QVariant.Double))
        fields.append(QgsField("population", QVariant.Int))

        feature = QgsFeature(fields, 5)
        feature.setGeometry(QgsGeometry(QgsPoint(5, 6)))
        feature.setAttributes(['Valsier Peninsula', 6.8, 198])

        expected = """{"name":"Valsier Peninsula",
"cost":6.8,
"population":198}"""
        self.assertEqual(QgsJsonUtils.exportAttributes(feature), expected)

        # test using field formatters
        source = QgsVectorLayer("Point?field=fldtxt:string&field=fldint:integer",
                                "parent", "memory")
        pf1 = QgsFeature()
        pf1.setFields(source.fields())
        pf1.setAttributes(["test1", 1])

        setup = QgsEditorWidgetSetup('ValueMap', {"map": {"one": 1, "two": 2, "three": 3}})
        source.setEditorWidgetSetup(1, setup)

        expected = """{"fldtxt":"test1",
"fldint":"one"}"""
        self.assertEqual(QgsJsonUtils.exportAttributes(pf1, source), expected)

    def testJSONExporter(self):
        """ test converting features to GeoJSON """
        fields = QgsFields()
        fields.append(QgsField("name", QVariant.String))
        fields.append(QgsField("cost", QVariant.Double))
        fields.append(QgsField("population", QVariant.Int))

        feature = QgsFeature(fields, 5)
        feature.setGeometry(QgsGeometry(QgsPoint(5, 6)))
        feature.setAttributes(['Valsier Peninsula', 6.8, 198])

        exporter = QgsJsonExporter()

        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        # test with linestring for bbox inclusion
        l = QgsLineString()
        l.setPoints([QgsPoint(5, 6), QgsPoint(15, 16)])
        feature.setGeometry(QgsGeometry(QgsLineString(l)))

        expected = """{
   "type":"Feature",
   "id":5,
   "bbox":[5, 6, 15, 16],
   "geometry":
   {"type": "LineString", "coordinates": [ [5, 6], [15, 16]]},
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        # test that precision is respected
        feature.setGeometry(QgsGeometry(QgsPoint(5.444444444, 6.333333333)))
        exporter.setPrecision(3)
        self.assertEqual(exporter.precision(), 3)
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5.444, 6.333]},
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)
        feature.setGeometry(QgsGeometry(QgsPoint(5, 6)))
        exporter.setPrecision(17)

        # test that attribute subset is respected
        exporter.setAttributes([0, 2])
        self.assertEqual(exporter.attributes(), [0, 2])
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "name":"Valsier Peninsula",
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        exporter.setAttributes([1])
        self.assertEqual(exporter.attributes(), [1])
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "cost":6.8
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)
        exporter.setAttributes([])

        # text excluding attributes

        exporter.setExcludedAttributes([1])
        self.assertEqual(exporter.excludedAttributes(), [1])
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "name":"Valsier Peninsula",
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        exporter.setExcludedAttributes([1, 2])
        self.assertEqual(exporter.excludedAttributes(), [1, 2])
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "name":"Valsier Peninsula"
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        exporter.setExcludedAttributes([0, 1, 2])
        self.assertEqual(exporter.excludedAttributes(), [0, 1, 2])
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":null
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        # test that excluded attributes take precedence over included

        exporter.setAttributes([1, 2])
        exporter.setExcludedAttributes([0, 1])
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        exporter.setAttributes([])
        exporter.setExcludedAttributes([])

        # test excluding geometry
        exporter.setIncludeGeometry(False)
        self.assertEqual(exporter.includeGeometry(), False)
        feature.setGeometry(QgsGeometry(QgsLineString(l)))

        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":null,
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)
        exporter.setIncludeGeometry(True)

        feature.setGeometry(QgsGeometry(QgsPoint(5, 6)))

        # test excluding attributes
        exporter.setIncludeAttributes(False)
        self.assertEqual(exporter.includeAttributes(), False)
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":null
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

        exporter.setIncludeGeometry(False)
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":null,
   "properties":null
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)
        exporter.setIncludeAttributes(True)

        # test overriding ID
        expected = """{
   "type":"Feature",
   "id":29,
   "geometry":null,
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198
   }
}"""
        self.assertEqual(exporter.exportFeature(feature, id=29), expected)

        # test injecting extra attributes
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":null,
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198,
      "extra":"val1",
      "extra2":2
   }
}"""
        self.assertEqual(exporter.exportFeature(feature, extraProperties={"extra": "val1", "extra2": 2}), expected)

        exporter.setIncludeAttributes(False)
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":null,
   "properties":{
      "extra":"val1",
      "extra2":{"nested_map":5,
"nested_map2":"val"},
      "extra3":[1,2,3]
   }
}"""
        expected2 = """{
   "type":"Feature",
   "id":5,
   "geometry":null,
   "properties":{
      "extra":"val1",
      "extra2":{"nested_map":5,"nested_map2":"val"},
      "extra3":[1,2,3]
   }
}"""
        exp_f = exporter.exportFeature(feature, extraProperties={"extra": "val1", "extra2": {"nested_map": 5, "nested_map2": "val"}, "extra3": [1, 2, 3]})
        self.assertTrue(exp_f == expected or exp_f == expected2)
        exporter.setIncludeGeometry(True)

    def testExportFeatureFieldFormatter(self):
        """ Test exporting a feature with formatting fields """

        # source layer
        source = QgsVectorLayer("Point?field=fldtxt:string&field=fldint:integer",
                                "parent", "memory")
        pr = source.dataProvider()
        pf1 = QgsFeature()
        pf1.setFields(source.fields())
        pf1.setAttributes(["test1", 1])
        pf2 = QgsFeature()
        pf2.setFields(source.fields())
        pf2.setAttributes(["test2", 2])
        assert pr.addFeatures([pf1, pf2])

        setup = QgsEditorWidgetSetup('ValueMap', {"map": {"one": 1, "two": 2, "three": 3}})
        source.setEditorWidgetSetup(1, setup)

        exporter = QgsJsonExporter()
        exporter.setVectorLayer(source)

        expected = """{
   "type":"Feature",
   "id":0,
   "geometry":null,
   "properties":{
      "fldtxt":"test1",
      "fldint":"one"
   }
}"""
        self.assertEqual(exporter.exportFeature(pf1), expected)

    def testExportFeatureCrs(self):
        """ Test CRS transform when exporting features """

        exporter = QgsJsonExporter()
        self.assertFalse(exporter.sourceCrs().isValid())

        #test layer
        layer = QgsVectorLayer("Point?crs=epsg:3111&field=fldtxt:string",
                               "parent", "memory")
        exporter = QgsJsonExporter(layer)
        self.assertTrue(exporter.sourceCrs().isValid())
        self.assertEqual(exporter.sourceCrs().authid(), 'EPSG:3111')

        exporter.setSourceCrs(QgsCoordinateReferenceSystem(3857, QgsCoordinateReferenceSystem.EpsgCrsId))
        self.assertTrue(exporter.sourceCrs().isValid())
        self.assertEqual(exporter.sourceCrs().authid(), 'EPSG:3857')

        # vector layer CRS should override
        exporter.setVectorLayer(layer)
        self.assertEqual(exporter.sourceCrs().authid(), 'EPSG:3111')

        # test that exported feature is reprojected
        feature = QgsFeature(layer.fields(), 5)
        feature.setGeometry(QgsGeometry(QgsPoint(2502577, 2403869)))
        feature.setAttributes(['test point'])

        # low precision, only need rough coordinate to check and don't want to deal with rounding errors
        exporter.setPrecision(1)
        expected = """{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [145, -37.9]},
   "properties":{
      "fldtxt":"test point"
   }
}"""
        self.assertEqual(exporter.exportFeature(feature), expected)

    def testExportFeatureRelations(self):
        """ Test exporting a feature with relations """

        #parent layer
        parent = QgsVectorLayer("Point?field=fldtxt:string&field=fldint:integer&field=foreignkey:integer",
                                "parent", "memory")
        pr = parent.dataProvider()
        pf1 = QgsFeature()
        pf1.setFields(parent.fields())
        pf1.setAttributes(["test1", 67, 123])
        pf2 = QgsFeature()
        pf2.setFields(parent.fields())
        pf2.setAttributes(["test2", 68, 124])
        assert pr.addFeatures([pf1, pf2])

        #child layer
        child = QgsVectorLayer(
            "Point?field=x:string&field=y:integer&field=z:integer",
            "referencedlayer", "memory")
        pr = child.dataProvider()
        f1 = QgsFeature()
        f1.setFields(child.fields())
        f1.setAttributes(["foo", 123, 321])
        f2 = QgsFeature()
        f2.setFields(child.fields())
        f2.setAttributes(["bar", 123, 654])
        f3 = QgsFeature()
        f3.setFields(child.fields())
        f3.setAttributes(["foobar", 124, 554])
        assert pr.addFeatures([f1, f2, f3])

        QgsProject.instance().addMapLayers([child, parent])

        rel = QgsRelation()
        rel.setId('rel1')
        rel.setName('relation one')
        rel.setReferencingLayer(child.id())
        rel.setReferencedLayer(parent.id())
        rel.addFieldPair('y', 'foreignkey')

        QgsProject.instance().relationManager().addRelation(rel)

        exporter = QgsJsonExporter()

        exporter.setVectorLayer(parent)
        self.assertEqual(exporter.vectorLayer(), parent)
        exporter.setIncludeRelated(True)
        self.assertEqual(exporter.includeRelated(), True)

        expected = """{
   "type":"Feature",
   "id":0,
   "geometry":null,
   "properties":{
      "fldtxt":"test1",
      "fldint":67,
      "foreignkey":123,
      "relation one":[{"x":"foo",
"y":123,
"z":321},
{"x":"bar",
"y":123,
"z":654}]
   }
}"""
        self.assertEqual(exporter.exportFeature(pf1), expected)

        expected = """{
   "type":"Feature",
   "id":0,
   "geometry":null,
   "properties":{
      "fldtxt":"test2",
      "fldint":68,
      "foreignkey":124,
      "relation one":[{"x":"foobar",
"y":124,
"z":554}]
   }
}"""
        self.assertEqual(exporter.exportFeature(pf2), expected)

        # with field formatter
        setup = QgsEditorWidgetSetup('ValueMap', {"map": {"apples": 123, "bananas": 124}})
        child.setEditorWidgetSetup(1, setup)
        expected = """{
   "type":"Feature",
   "id":0,
   "geometry":null,
   "properties":{
      "fldtxt":"test1",
      "fldint":67,
      "foreignkey":123,
      "relation one":[{"x":"foo",
"y":"apples",
"z":321},
{"x":"bar",
"y":"apples",
"z":654}]
   }
}"""
        self.assertEqual(exporter.exportFeature(pf1), expected)

        # test excluding related attributes
        exporter.setIncludeRelated(False)
        self.assertEqual(exporter.includeRelated(), False)

        expected = """{
   "type":"Feature",
   "id":0,
   "geometry":null,
   "properties":{
      "fldtxt":"test2",
      "fldint":68,
      "foreignkey":124
   }
}"""
        self.assertEqual(exporter.exportFeature(pf2), expected)

        # test without vector layer set
        exporter.setIncludeRelated(True)
        exporter.setVectorLayer(None)

        expected = """{
   "type":"Feature",
   "id":0,
   "geometry":null,
   "properties":{
      "fldtxt":"test2",
      "fldint":68,
      "foreignkey":124
   }
}"""
        self.assertEqual(exporter.exportFeature(pf2), expected)

    def testExportFeatures(self):
        """ Test exporting feature collections """

        fields = QgsFields()
        fields.append(QgsField("name", QVariant.String))
        fields.append(QgsField("cost", QVariant.Double))
        fields.append(QgsField("population", QVariant.Int))

        feature = QgsFeature(fields, 5)
        feature.setGeometry(QgsGeometry(QgsPoint(5, 6)))
        feature.setAttributes(['Valsier Peninsula', 6.8, 198])

        exporter = QgsJsonExporter()

        # single feature
        expected = """{ "type": "FeatureCollection",
    "features":[
{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198
   }
}
]}"""
        self.assertEqual(exporter.exportFeatures([feature]), expected)

        # multiple features
        feature2 = QgsFeature(fields, 6)
        feature2.setGeometry(QgsGeometry(QgsPoint(7, 8)))
        feature2.setAttributes(['Henry Gale Island', 9.7, 38])

        expected = """{ "type": "FeatureCollection",
    "features":[
{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198
   }
},
{
   "type":"Feature",
   "id":6,
   "geometry":
   {"type": "Point", "coordinates": [7, 8]},
   "properties":{
      "name":"Henry Gale Island",
      "cost":9.7,
      "population":38
   }
}
]}"""
        self.assertEqual(exporter.exportFeatures([feature, feature2]), expected)

    def testExportFeaturesWithLocale_regression20053(self):
        """ Test exporting feature export with range widgets and locale different than C
        Regression: https://issues.qgis.org/issues/20053 - decimal separator in csv files
        """

        source = QgsVectorLayer("Point?field=name:string&field=cost:double&field=population:int&field=date:date",
                                "parent", "memory")
        self.assertTrue(source.isValid())
        fields = source.fields()

        feature = QgsFeature(fields, 5)
        feature.setGeometry(QgsGeometry(QgsPoint(5, 6)))
        feature.setAttributes(['Valsier Peninsula', 6.8, 198000, '2018-09-10'])

        exporter = QgsJsonExporter()

        # single feature
        expected = """{ "type": "FeatureCollection",
    "features":[
{
   "type":"Feature",
   "id":5,
   "geometry":
   {"type": "Point", "coordinates": [5, 6]},
   "properties":{
      "name":"Valsier Peninsula",
      "cost":6.8,
      "population":198000,
      "date":"2018-09-10"
   }
}
]}"""
        self.assertEqual(exporter.exportFeatures([feature]), expected)

        setup = QgsEditorWidgetSetup('Range', {
            'AllowNull': True,
            'Max': 2147483647,
            'Min': -2147483648,
            'Precision': 4,
            'Step': 1,
            'Style': 'SpinBox'
        }
        )
        source.setEditorWidgetSetup(1, setup)
        source.setEditorWidgetSetup(2, setup)

        QLocale.setDefault(QLocale('it'))
        exporter.setVectorLayer(source)
        self.assertEqual(exporter.exportFeatures([feature]), expected)

    def testExportFieldAlias(self):
        """ Test exporting a feature with fields' alias """

        # source layer
        source = QgsVectorLayer("Point?field=fldtxt:string&field=fldint:integer",
                                "parent", "memory")
        pr = source.dataProvider()
        pf1 = QgsFeature()
        pf1.setFields(source.fields())
        pf1.setAttributes(["test1", 1])
        pf2 = QgsFeature()
        pf2.setFields(source.fields())
        pf2.setAttributes(["test2", 2])
        assert pr.addFeatures([pf1, pf2])

        source.setFieldAlias(0, "alias_fldtxt")
        source.setFieldAlias(1, "alias_fldint")

        exporter = QgsJsonExporter()
        exporter.setAttributeDisplayName(True)
        exporter.setVectorLayer(source)

        expected = """{
   "type":"Feature",
   "id":0,
   "geometry":null,
   "properties":{
      "alias_fldtxt":"test1",
      "alias_fldint":1
   }
}"""
        self.assertEqual(exporter.exportFeature(pf1), expected)


if __name__ == "__main__":
    unittest.main()
