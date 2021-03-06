# -*- coding: utf-8 -*-
"""QGIS Unit tests for postgres transaction groups.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '11/06/2018'
__copyright__ = 'Copyright 2018, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

import os

from qgis.core import (
    QgsVectorLayer,
    QgsProject,
    QgsTransaction,
    QgsDataSourceUri
)

from qgis.testing import start_app, unittest

start_app()


class TestQgsPostgresTransaction(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Setup the involved layers and relations for a n:m relation
        :return:
        """
        cls.dbconn = 'dbname=\'qgis_test\''
        if 'QGIS_PGTEST_DB' in os.environ:
            cls.dbconn = os.environ['QGIS_PGTEST_DB']
        # Create test layer
        cls.vl_b = QgsVectorLayer(cls.dbconn + ' sslmode=disable key=\'pk\' table="qgis_test"."books" sql=', 'books',
                                  'postgres')
        cls.vl_a = QgsVectorLayer(cls.dbconn + ' sslmode=disable key=\'pk\' table="qgis_test"."authors" sql=',
                                  'authors', 'postgres')

        QgsProject.instance().addMapLayer(cls.vl_b)
        QgsProject.instance().addMapLayer(cls.vl_a)

        cls.relMgr = QgsProject.instance().relationManager()

        assert (cls.vl_a.isValid())
        assert (cls.vl_b.isValid())

    def startTransaction(self):
        """
        Start a new transaction and set all layers into transaction mode.

        :return: None
        """
        lyrs = [self.vl_a, self.vl_b]

        self.transaction = QgsTransaction.create(lyrs)
        self.transaction.begin()
        for l in lyrs:
            l.startEditing()

    def rollbackTransaction(self):
        """
        Rollback all changes done in this transaction.
        We always rollback and never commit to have the database in a pristine
        state at the end of each test.

        :return: None
        """
        lyrs = [self.vl_a, self.vl_b]
        for l in lyrs:
            l.commitChanges()
        self.transaction.rollback()

    def test_transactionsGroup(self):
        conn_string = QgsDataSourceUri(self.vl_b.source()).connectionInfo()

        # No transaction group.
        QgsProject.instance().setAutoTransaction(False)
        noTg = QgsProject.instance().transactionGroup("postgres", conn_string)
        self.assertIsNone(noTg)

        # start transaction - no auto transaction
        self.startTransaction()
        noTg = QgsProject.instance().transactionGroup("postgres", conn_string)
        self.assertIsNone(noTg)
        self.rollbackTransaction()

        # with auto transactions
        QgsProject.instance().setAutoTransaction(True)
        self.startTransaction()
        noTg = QgsProject.instance().transactionGroup("postgres", conn_string)
        self.assertIsNotNone(noTg)
        self.rollbackTransaction()

        # bad provider key
        self.startTransaction()
        noTg = QgsProject.instance().transactionGroup("xxxpostgres", conn_string)
        self.assertIsNone(noTg)
        self.rollbackTransaction()


if __name__ == '__main__':
    unittest.main()
