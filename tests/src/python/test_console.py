# -*- coding: utf-8 -*-
"""QGIS Unit tests for the console

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Matthias Kuhn'
__date__ = '15.4.2016'
__copyright__ = 'Copyright 2015, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA
import os

from qgis.testing import unittest, start_app
from console import console
from qgis.core import QgsSettings
from qgis.PyQt.QtCore import QCoreApplication

start_app()


class TestConsole(unittest.TestCase):

    def setUp(self):
        QgsSettings().setValue('pythonConsole/contextHelpOnFirstLaunch', False)

    def test_show_console(self):
        if os.name == 'nt':
            QCoreApplication.setOrganizationName("QGIS")
            QCoreApplication.setOrganizationDomain("qgis.org")
            QCoreApplication.setApplicationName("QGIS-TEST")

        my_console = console.show_console()
        my_console_widget = my_console.console


if __name__ == "__main__":
    unittest.main()
