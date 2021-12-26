# -*- coding: utf-8 -*-
"""QGIS Unit tests for the Oracle provider.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '2016-07-06'
__copyright__ = 'Copyright 2016, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

import os

from qgis.core import QgsSettings, QgsVectorLayer, QgsFeatureRequest, NULL

from qgis.PyQt.QtCore import QDate, QTime, QDateTime, QVariant
from qgis.PyQt.QtSql import QSqlDatabase, QSqlQuery

from utilities import unitTestDataPath, compareWkt
from qgis.testing import start_app, unittest
from providertestbase import ProviderTestCase

start_app()
TEST_DATA_DIR = unitTestDataPath()


class TestPyQgsOracleProvider(unittest.TestCase, ProviderTestCase):

    @classmethod
    def setUpClass(cls):
        """Run before all tests"""
        cls.dbconn = "host=localhost port=1521 user='QGIS' password='qgis'"
        if 'QGIS_ORACLETEST_DB' in os.environ:
            cls.dbconn = os.environ['QGIS_ORACLETEST_DB']
        # Create test layers
        cls.vl = QgsVectorLayer(
            cls.dbconn + ' sslmode=disable key=\'pk\' srid=4326 type=POINT table="QGIS"."SOME_DATA" (GEOM) sql=', 'test', 'oracle')
        assert(cls.vl.isValid())
        cls.source = cls.vl.dataProvider()
        cls.poly_vl = QgsVectorLayer(
            cls.dbconn + ' sslmode=disable key=\'pk\' srid=4326 type=POLYGON table="QGIS"."SOME_POLY_DATA" (GEOM) sql=', 'test', 'oracle')
        assert(cls.poly_vl.isValid())
        cls.poly_provider = cls.poly_vl.dataProvider()

        cls.conn = QSqlDatabase.addDatabase('QOCISPATIAL', "oracletest")
        cls.conn.setDatabaseName('10.0.0.2/orcl')
        if 'QGIS_ORACLETEST_DBNAME' in os.environ:
            cls.conn.setDatabaseName('QGIS_ORACLETEST_DBNAME')
        cls.conn.setUserName('QGIS')
        cls.conn.setPassword('qgis')
        assert cls.conn.open()

    @classmethod
    def tearDownClass(cls):
        """Run after all tests"""

    def execSQLCommand(self, sql, ignore_errors=False):
        self.assertTrue(self.conn)
        query = QSqlQuery(self.conn)
        self.assertTrue(query.exec_(sql), sql + ': ' + query.lastError().text())
        query.finish()

    def getSource(self):
        # create temporary table for edit tests
        self.execSQLCommand('DROP TABLE "QGIS"."EDIT_DATA"', ignore_errors=True)
        self.execSQLCommand("""CREATE TABLE QGIS.EDIT_DATA ("pk" INTEGER GENERATED by default ON null as IDENTITY(START WITH 1 INCREMENT BY 1) PRIMARY KEY, "cnt" INTEGER, "name" VARCHAR2(100), "name2" VARCHAR2(100), "num_char" VARCHAR2(100), GEOM SDO_GEOMETRY)""")
        self.execSQLCommand("""INSERT INTO QGIS.EDIT_DATA ("pk", "cnt", "name", "name2", "num_char", GEOM)
      SELECT 5, -200, NULL, 'NuLl', '5', SDO_GEOMETRY( 2001,4326,SDO_POINT_TYPE(-71.123, 78.23, NULL), NULL, NULL) from dual
  UNION ALL SELECT 3,  300, 'Pear', 'PEaR', '3', NULL from dual
  UNION ALL SELECT 1,  100, 'Orange', 'oranGe', '1', SDO_GEOMETRY( 2001,4326,SDO_POINT_TYPE(-70.332, 66.33, NULL), NULL, NULL) from dual
  UNION ALL SELECT 2,  200, 'Apple', 'Apple', '2', SDO_GEOMETRY( 2001,4326,SDO_POINT_TYPE(-68.2, 70.8, NULL), NULL, NULL) from dual
  UNION ALL SELECT 4,  400, 'Honey', 'Honey', '4', SDO_GEOMETRY( 2001,4326,SDO_POINT_TYPE(-65.32, 78.3, NULL), NULL, NULL) from dual""")
        vl = QgsVectorLayer(
            self.dbconn + ' sslmode=disable key=\'pk\' srid=4326 type=POINT table="QGIS"."EDIT_DATA" (GEOM) sql=',
            'test', 'oracle')
        return vl

    def getEditableLayer(self):
        return self.getSource()

    def enableCompiler(self):
        QgsSettings().setValue('/qgis/compileExpressions', True)
        return True

    def disableCompiler(self):
        QgsSettings().setValue('/qgis/compileExpressions', False)

    def uncompiledFilters(self):
        filters = set([
            '(name = \'Apple\') is not null',
            '"name" || \' \' || "name" = \'Orange Orange\'',
            '"name" || \' \' || "cnt" = \'Orange 100\'',
            '\'x\' || "name" IS NOT NULL',
            '\'x\' || "name" IS NULL',
            'false and NULL',
            'true and NULL',
            'NULL and false',
            'NULL and true',
            'NULL and NULL',
            'false or NULL',
            'true or NULL',
            'NULL or false',
            'NULL or true',
            'NULL or NULL',
            'not null',
            'radians(cnt) < 2',
            'degrees(pk) <= 200',
            'atan2(3.14, pk) < 1',
            'pk < pi()',
            'log10(pk) < 0.5',
            'pk < pi() / 2',
            'pk = char(51)',
            'pk = coalesce(NULL,3,4)',
            'name = trim(\'   Apple   \')',
            'x($geometry) < -70',
            'y($geometry) > 70',
            'xmin($geometry) < -70',
            'ymin($geometry) > 70',
            'xmax($geometry) < -70',
            'ymax($geometry) > 70',
            'disjoint($geometry,geom_from_wkt( \'Polygon ((-72.2 66.1, -65.2 66.1, -65.2 72.0, -72.2 72.0, -72.2 66.1))\'))',
            'intersects($geometry,geom_from_wkt( \'Polygon ((-72.2 66.1, -65.2 66.1, -65.2 72.0, -72.2 72.0, -72.2 66.1))\'))',
            'contains(geom_from_wkt( \'Polygon ((-72.2 66.1, -65.2 66.1, -65.2 72.0, -72.2 72.0, -72.2 66.1))\'),$geometry)',
            'distance($geometry,geom_from_wkt( \'Point (-70 70)\')) > 7',
            'intersects($geometry,geom_from_gml( \'<gml:Polygon srsName="EPSG:4326"><gml:outerBoundaryIs><gml:LinearRing><gml:coordinates>-72.2,66.1 -65.2,66.1 -65.2,72.0 -72.2,72.0 -72.2,66.1</gml:coordinates></gml:LinearRing></gml:outerBoundaryIs></gml:Polygon>\'))',
            'x($geometry) < -70',
            'y($geometry) > 79',
            'xmin($geometry) < -70',
            'ymin($geometry) < 76',
            'xmax($geometry) > -68',
            'ymax($geometry) > 80',
            'area($geometry) > 10',
            'perimeter($geometry) < 12',
            'relate($geometry,geom_from_wkt( \'Polygon ((-68.2 82.1, -66.95 82.1, -66.95 79.05, -68.2 79.05, -68.2 82.1))\')) = \'FF2FF1212\'',
            'relate($geometry,geom_from_wkt( \'Polygon ((-68.2 82.1, -66.95 82.1, -66.95 79.05, -68.2 79.05, -68.2 82.1))\'), \'****F****\')',
            'crosses($geometry,geom_from_wkt( \'Linestring (-68.2 82.1, -66.95 82.1, -66.95 79.05)\'))',
            'overlaps($geometry,geom_from_wkt( \'Polygon ((-68.2 82.1, -66.95 82.1, -66.95 79.05, -68.2 79.05, -68.2 82.1))\'))',
            'within($geometry,geom_from_wkt( \'Polygon ((-75.1 76.1, -75.1 81.6, -68.8 81.6, -68.8 76.1, -75.1 76.1))\'))',
            'overlaps(translate($geometry,-1,-1),geom_from_wkt( \'Polygon ((-75.1 76.1, -75.1 81.6, -68.8 81.6, -68.8 76.1, -75.1 76.1))\'))',
            'overlaps(buffer($geometry,1),geom_from_wkt( \'Polygon ((-75.1 76.1, -75.1 81.6, -68.8 81.6, -68.8 76.1, -75.1 76.1))\'))',
            'intersects(centroid($geometry),geom_from_wkt( \'Polygon ((-74.4 78.2, -74.4 79.1, -66.8 79.1, -66.8 78.2, -74.4 78.2))\'))',
            'intersects(point_on_surface($geometry),geom_from_wkt( \'Polygon ((-74.4 78.2, -74.4 79.1, -66.8 79.1, -66.8 78.2, -74.4 78.2))\'))'
        ])
        return filters

    def testAddFeatureWrongGeomType(self):
        """
        We override this test for Oracle provider, because Oracle DBs don't care
        about geometry type constraints
        """
        pass

    # HERE GO THE PROVIDER SPECIFIC TESTS
    def testDateTimeTypes(self):
        vl = QgsVectorLayer('%s table="QGIS"."DATE_TIMES" sql=' %
                            (self.dbconn), "testdatetimes", "oracle")
        self.assertTrue(vl.isValid())

        fields = vl.dataProvider().fields()
        self.assertEqual(fields.at(fields.indexFromName(
            'date_field')).type(), QVariant.Date)
        self.assertEqual(fields.at(fields.indexFromName(
            'datetime_field')).type(), QVariant.DateTime)

        f = next(vl.getFeatures(QgsFeatureRequest()))

        date_idx = vl.fields().lookupField('date_field')
        self.assertIsInstance(f.attributes()[date_idx], QDate)
        self.assertEqual(f.attributes()[date_idx], QDate(2004, 3, 4))
        datetime_idx = vl.fields().lookupField('datetime_field')
        self.assertIsInstance(f.attributes()[datetime_idx], QDateTime)
        self.assertEqual(f.attributes()[datetime_idx], QDateTime(
            QDate(2004, 3, 4), QTime(13, 41, 52)))

    def testDefaultValue(self):
        self.assertEqual(self.source.defaultValue(1), NULL)
        self.assertEqual(self.source.defaultValue(2), "'qgis'")

    def testPoints(self):
        vl = QgsVectorLayer('%s table="QGIS"."POINT_DATA" (GEOM) srid=4326 type=POINT sql=' %
                            (self.dbconn), "testpoints", "oracle")
        self.assertTrue(vl.isValid())

        features = [f for f in vl.getFeatures()]
        self.assertEqual(features[0].geometry().asWkt(), 'Point (1 2)')
        self.assertEqual(features[1].geometry().asWkt(), 'PointZ (1 2 3)')
        self.assertEqual(features[2].geometry().asWkt(), 'MultiPointZ ((1 2 3),(4 5 6))')
        self.assertEqual(features[3].geometry().asWkt(), 'MultiPoint ((1 2),(3 4))')
        self.assertEqual(features[4].geometry().asWkt(), 'MultiPointZ ((1 2 3),(4 5 6))')
        self.assertEqual(features[5].geometry().asWkt(), 'Point (1 2)')
        self.assertEqual(features[6].geometry().asWkt(), 'Point (3 4)')
        self.assertEqual(features[7].geometry().asWkt(), 'Point (5 6)')

    def testCurves(self):
        vl = QgsVectorLayer('%s table="QGIS"."LINE_DATA" (GEOM) srid=4326 type=LINESTRING sql=' %
                            (self.dbconn), "testlines", "oracle")
        self.assertTrue(vl.isValid())

        features = {f['pk']: f for f in vl.getFeatures()}
        self.assertTrue(compareWkt(features[1].geometry().asWkt(), 'LineString (1 2, 3 4, 5 6)', 0.00001), features[1].geometry().asWkt())
        self.assertTrue(compareWkt(features[2].geometry().asWkt(), 'CircularString (1 2, 5 4, 7 2.2, 10 0.1, 13 4)', 0.00001), features[2].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[3].geometry().asWkt(), 'CompoundCurve ((-1 -5, 1 2),CircularString (1 2, 5 4, 7 2.20, 10 0.1, 13 4),(13 4, 17 -6))', 0.00001), features[3].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[4].geometry().asWkt(), 'LineStringZ (1 2 3, 4 5 6, 7 8 9)', 0.00001), features[4].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[5].geometry().asWkt(), 'MultiLineString ((1 2, 3 4),(5 6, 7 8, 9 10))', 0.00001), features[5].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[6].geometry().asWkt(), 'MultiLineStringZ ((1 2 11, 3 4 -11),(5 6 9, 7 8 1, 9 10 -3))', 0.00001), features[6].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[7].geometry().asWkt(), 'MultiCurve (CircularString (1 2, 5 4, 7 2.2, 10 0.1, 13 4),CircularString (-11 -3, 5 7, 10 -1))', 0.00001), features[7].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[8].geometry().asWkt(), 'MultiCurve (CompoundCurve ((-1 -5, 1 2),CircularString (1 2, 5 4, 7 2.2, 10 0.1, 13 4),(13 4, 17 -6)),CompoundCurve (CircularString (1 2, 5 4, 7 2.2, 10 0.1, 13 4)),CompoundCurve ((-11 -3, 5 7, 10 -1)))', 0.00001), features[8].geometry().asWkt())

    def testSurfaces(self):
        vl = QgsVectorLayer('%s table="QGIS"."POLY_DATA" (GEOM) srid=4326 type=POLYGON sql=' %
                            (self.dbconn), "testpoly", "oracle")
        self.assertTrue(vl.isValid())

        features = {f['pk']: f for f in vl.getFeatures()}
        self.assertTrue(compareWkt(features[1].geometry().asWkt(), 'Polygon ((1 2, 11 2, 11 22, 1 22, 1 2))', 0.00001), features[1].geometry().asWkt())
        self.assertTrue(compareWkt(features[2].geometry().asWkt(), 'PolygonZ ((1 2 3, 11 2 13, 11 22 15, 1 22 7, 1 2 3))', 0.00001), features[2].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[3].geometry().asWkt(), 'Polygon ((1 2, 11 2, 11 22, 1 22, 1 2),(5 6, 8 9, 8 6, 5 6),(3 4, 5 6, 3 6, 3 4))', 0.00001), features[3].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[4].geometry().asWkt(), 'PolygonZ ((1 2 3, 11 2 13, 11 22 15, 1 22 7, 1 2 3),(5 6 1, 8 9 -1, 8 6 2, 5 6 1))', 0.00001), features[4].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[5].geometry().asWkt(), 'Polygon ((1 2, 11 2, 11 22, 1 22, 1 2))', 0.00001), features[5].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[6].geometry().asWkt(), 'CurvePolygon (CircularString (6.76923076923076916 22.82875364393326834, 17.98259979777942519 11.61538461538461497, 6.76923076923076916 0.40201558683595984, -4.44413825931788598 11.61538461538461497, 6.76923076923076916 22.82875364393326834))', 0.00001), features[6].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[7].geometry().asWkt(), 'MultiPolygon (((1 2, 11 2, 11 22, 1 22, 1 2)),((1 2, 11 2, 11 22, 1 22, 1 2),(5 6, 8 9, 8 6, 5 6),(3 4, 5 6, 3 6, 3 4)))', 0.00001), features[7].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[8].geometry().asWkt(), 'MultiPolygonZ (((1 2 3, 11 2 13, 11 22 15, 1 22 7, 1 2 3)),((1 2 3, 11 2 13, 11 22 15, 1 22 7, 1 2 3),(5 6 1, 8 9 -1, 8 6 2, 5 6 1)))', 0.00001), features[8].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[9].geometry().asWkt(), 'CurvePolygon (CircularString (1 3, 3 5, 4 7, 7 3, 1 3))', 0.00001), features[9].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[10].geometry().asWkt(), 'CurvePolygon (CircularString (1 3, 3 5, 4 7, 7 3, 1 3),CircularString (3.1 3.3, 3.3 3.5, 3.4 3.7, 3.7 3.3, 3.1 3.3))', 0.00001), features[10].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[11].geometry().asWkt(), 'CurvePolygon(CompoundCurve ((-1 -5, 1 2),CircularString (1 2, 5 4, 7 2.20, 10 0.1, 13 4),(13 4, 17 -6),CircularString (17 -6, 5 -7, -1 -5)))', 0.00001), features[11].geometry().asWkt())
        self.assertTrue(
            compareWkt(features[12].geometry().asWkt(), 'MultiSurface (CurvePolygon (CircularString (1 3, 3 5, 4 7, 7 3, 1 3)),CurvePolygon (CircularString (11 3, 13 5, 14 7, 17 3, 11 3)))', 0.00001), features[12].geometry().asWkt())


if __name__ == '__main__':
    unittest.main()
