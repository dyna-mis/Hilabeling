# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsFeedback.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '12/02/2017'
__copyright__ = 'Copyright 2017, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

from qgis.core import (QgsFeedback)
from qgis.PyQt.QtTest import QSignalSpy
from qgis.testing import unittest


class TestQgsFeedback(unittest.TestCase):

    def testCancel(self):
        f = QgsFeedback()
        self.assertFalse(f.isCanceled())

        cancel_spy = QSignalSpy(f.canceled)

        f.cancel()
        self.assertTrue(f.isCanceled())
        self.assertEqual(len(cancel_spy), 1)

    def testProgress(self):
        f = QgsFeedback()
        self.assertEqual(f.progress(), 0.0)

        progress_spy = QSignalSpy(f.progressChanged)

        f.setProgress(25)
        self.assertEqual(f.progress(), 25.0)
        self.assertEqual(len(progress_spy), 1)
        self.assertEqual(progress_spy[0][0], 25.0)


if __name__ == '__main__':
    unittest.main()
