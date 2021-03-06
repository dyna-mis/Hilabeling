# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsLayoutItemPage.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = '(C) 2017 by Nyall Dawson'
__date__ = '23/10/2017'
__copyright__ = 'Copyright 2017, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

from qgis.testing import start_app, unittest
from qgis.core import QgsLayoutItemPage

from test_qgslayoutitem import LayoutItemTestCase

start_app()


class TestQgsLayoutPage(unittest.TestCase, LayoutItemTestCase):

    @classmethod
    def setUpClass(cls):
        cls.item_class = QgsLayoutItemPage


if __name__ == '__main__':
    unittest.main()
