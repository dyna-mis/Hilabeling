# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsRange

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '11.04.2017'
__copyright__ = 'Copyright 2017, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import qgis  # NOQA

from qgis.testing import unittest
from qgis.core import (QgsIntRange,
                       QgsDateRange)
from qgis.PyQt.QtCore import QDate


class TestQgsIntRange(unittest.TestCase):

    def testGetters(self):
        range = QgsIntRange(1, 11)
        self.assertEqual(range.lower(), 1)
        self.assertEqual(range.upper(), 11)
        self.assertTrue(range.includeLower())
        self.assertTrue(range.includeUpper())

        range = QgsIntRange(-1, 3, False, False)
        self.assertEqual(range.lower(), -1)
        self.assertEqual(range.upper(), 3)
        self.assertFalse(range.includeLower())
        self.assertFalse(range.includeUpper())

    def testIsEmpty(self):
        range = QgsIntRange(1, 1)
        # should not be empty because 1 is included
        self.assertFalse(range.isEmpty())

        range = QgsIntRange(1, 1, False, False)
        # should be empty because 1 is NOT included
        self.assertTrue(range.isEmpty())

        # invalid range is empty
        range = QgsIntRange(1, -1)
        self.assertTrue(range.isEmpty())

    def testIsSingleton(self):
        range = QgsIntRange(1, 1)
        self.assertTrue(range.isSingleton())

        range = QgsIntRange(1, 10)
        self.assertFalse(range.isSingleton())

        range = QgsIntRange(1, 1, False, False)
        # should not be singleton because 1 is NOT included
        self.assertFalse(range.isSingleton())

        # invalid range is not singleton
        range = QgsIntRange(1, -1)
        self.assertFalse(range.isSingleton())

    def testContains(self):
        # includes both ends
        range = QgsIntRange(0, 10)
        self.assertTrue(range.contains(QgsIntRange(1, 9)))
        self.assertTrue(range.contains(QgsIntRange(1, 10)))
        self.assertTrue(range.contains(QgsIntRange(0, 9)))
        self.assertTrue(range.contains(QgsIntRange(0, 10)))
        self.assertFalse(range.contains(QgsIntRange(-1, 9)))
        self.assertFalse(range.contains(QgsIntRange(1, 11)))

        # does not include left end
        range = QgsIntRange(0, 10, False, True)
        self.assertTrue(range.contains(QgsIntRange(1, 9)))
        self.assertTrue(range.contains(QgsIntRange(1, 10)))
        self.assertFalse(range.contains(QgsIntRange(0, 9)))
        self.assertFalse(range.contains(QgsIntRange(0, 10)))
        self.assertFalse(range.contains(QgsIntRange(-1, 9)))
        self.assertFalse(range.contains(QgsIntRange(1, 11)))

        # does not include right end
        range = QgsIntRange(0, 10, True, False)
        self.assertTrue(range.contains(QgsIntRange(1, 9)))
        self.assertFalse(range.contains(QgsIntRange(1, 10)))
        self.assertTrue(range.contains(QgsIntRange(0, 9)))
        self.assertFalse(range.contains(QgsIntRange(0, 10)))
        self.assertFalse(range.contains(QgsIntRange(-1, 9)))
        self.assertFalse(range.contains(QgsIntRange(1, 11)))

    def testContainsElement(self):
        # includes both ends
        range = QgsIntRange(0, 10)
        self.assertTrue(range.contains(0))
        self.assertTrue(range.contains(5))
        self.assertTrue(range.contains(10))
        self.assertFalse(range.contains(-1))
        self.assertFalse(range.contains(11))

        # includes left end
        range = QgsIntRange(0, 10, True, False)
        self.assertTrue(range.contains(0))
        self.assertTrue(range.contains(5))
        self.assertFalse(range.contains(10))
        self.assertFalse(range.contains(-1))
        self.assertFalse(range.contains(11))

        # includes right end
        range = QgsIntRange(0, 10, False, True)
        self.assertFalse(range.contains(0))
        self.assertTrue(range.contains(5))
        self.assertTrue(range.contains(10))
        self.assertFalse(range.contains(-1))
        self.assertFalse(range.contains(11))

        # includes neither end
        range = QgsIntRange(0, 10, False, False)
        self.assertFalse(range.contains(0))
        self.assertTrue(range.contains(5))
        self.assertFalse(range.contains(10))
        self.assertFalse(range.contains(-1))
        self.assertFalse(range.contains(11))

    def testOverlaps(self):
        # includes both ends
        range = QgsIntRange(0, 10)
        self.assertTrue(range.overlaps(QgsIntRange(1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(10, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 0)))
        self.assertFalse(range.overlaps(QgsIntRange(-10, -1)))
        self.assertFalse(range.overlaps(QgsIntRange(11, 12)))

        # includes left end
        range = QgsIntRange(0, 10, True, False)
        self.assertTrue(range.overlaps(QgsIntRange(1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 11)))
        self.assertFalse(range.overlaps(QgsIntRange(10, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 0)))
        self.assertFalse(range.overlaps(QgsIntRange(-10, -1)))
        self.assertFalse(range.overlaps(QgsIntRange(11, 12)))

        # includes right end
        range = QgsIntRange(0, 10, False, True)
        self.assertTrue(range.overlaps(QgsIntRange(1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(10, 11)))
        self.assertFalse(range.overlaps(QgsIntRange(-1, 0)))
        self.assertFalse(range.overlaps(QgsIntRange(-10, -1)))
        self.assertFalse(range.overlaps(QgsIntRange(11, 12)))

        # includes neither end
        range = QgsIntRange(0, 10, False, False)
        self.assertTrue(range.overlaps(QgsIntRange(1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(0, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 10)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 9)))
        self.assertTrue(range.overlaps(QgsIntRange(1, 11)))
        self.assertTrue(range.overlaps(QgsIntRange(-1, 11)))
        self.assertFalse(range.overlaps(QgsIntRange(10, 11)))
        self.assertFalse(range.overlaps(QgsIntRange(-1, 0)))
        self.assertFalse(range.overlaps(QgsIntRange(-10, -1)))
        self.assertFalse(range.overlaps(QgsIntRange(11, 12)))


class TestQgsDateRange(unittest.TestCase):

    def testGetters(self):
        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))
        self.assertEqual(range.begin(), QDate(2010, 3, 1))
        self.assertEqual(range.end(), QDate(2010, 6, 2))
        self.assertTrue(range.includeBeginning())
        self.assertTrue(range.includeEnd())

        range = QgsDateRange(QDate(), QDate(2010, 6, 2))
        self.assertFalse(range.begin().isValid())
        self.assertEqual(range.end(), QDate(2010, 6, 2))

        range = QgsDateRange(QDate(2010, 3, 1), QDate())
        self.assertEqual(range.begin(), QDate(2010, 3, 1))
        self.assertFalse(range.end().isValid())

    def testIsEmpty(self):
        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))
        self.assertFalse(range.isEmpty())

        range = QgsDateRange(QDate(), QDate(2010, 6, 2))
        self.assertFalse(range.isEmpty())

        range = QgsDateRange(QDate(2010, 3, 1), QDate())
        self.assertFalse(range.isEmpty())

        # check QgsDateRange docs - this is treated as an infinite range, so is NOT empty
        range = QgsDateRange(QDate(), QDate())
        self.assertFalse(range.isEmpty())

        range = QgsDateRange(QDate(2017, 3, 1), QDate(2010, 6, 2))
        self.assertTrue(range.isEmpty())

        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 3, 1))
        self.assertFalse(range.isEmpty())

        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 3, 1), False, False)
        self.assertTrue(range.isEmpty())

    def testContains(self):
        # includes both ends
        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))))
        self.assertFalse(range.contains(QgsDateRange(QDate(2009, 4, 1), QDate(2010, 4, 5))))
        self.assertFalse(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2017, 4, 5))))
        self.assertFalse(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate())))
        self.assertFalse(range.contains(QgsDateRange(QDate(), QDate(2010, 4, 1))))

        # infinite left end
        range = QgsDateRange(QDate(), QDate(2010, 6, 2))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2009, 4, 1), QDate(2010, 4, 5))))
        self.assertFalse(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2017, 4, 5))))
        self.assertFalse(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate())))
        self.assertTrue(range.contains(QgsDateRange(QDate(), QDate(2010, 4, 1))))

        # infinite right end
        range = QgsDateRange(QDate(2010, 3, 1), QDate())
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))))
        self.assertFalse(range.contains(QgsDateRange(QDate(2009, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate(2017, 4, 5))))
        self.assertTrue(range.contains(QgsDateRange(QDate(2010, 4, 1), QDate())))
        self.assertFalse(range.contains(QgsDateRange(QDate(), QDate(2010, 4, 1))))

    def testContainsElement(self):
        # includes both ends
        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))
        self.assertTrue(range.contains(QDate(2010, 3, 1)))
        self.assertTrue(range.contains(QDate(2010, 5, 2)))
        self.assertTrue(range.contains(QDate(2010, 6, 2)))
        self.assertFalse(range.contains(QDate(2009, 6, 2)))
        self.assertFalse(range.contains(QDate(2017, 6, 2)))
        self.assertFalse(range.contains(QDate()))

        # infinite left end
        range = QgsDateRange(QDate(), QDate(2010, 6, 2))
        self.assertTrue(range.contains(QDate(2010, 3, 1)))
        self.assertTrue(range.contains(QDate(2010, 5, 2)))
        self.assertTrue(range.contains(QDate(2010, 6, 2)))
        self.assertTrue(range.contains(QDate(2009, 6, 2)))
        self.assertFalse(range.contains(QDate(2017, 6, 2)))
        self.assertFalse(range.contains(QDate()))

        # infinite right end
        range = QgsDateRange(QDate(2010, 3, 1), QDate())
        self.assertTrue(range.contains(QDate(2010, 3, 1)))
        self.assertTrue(range.contains(QDate(2010, 5, 2)))
        self.assertTrue(range.contains(QDate(2010, 6, 2)))
        self.assertFalse(range.contains(QDate(2009, 6, 2)))
        self.assertTrue(range.contains(QDate(2017, 6, 2)))
        self.assertFalse(range.contains(QDate()))

    def testOverlaps(self):
        # includes both ends
        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2009, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2017, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate())))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(), QDate(2010, 4, 1))))
        self.assertFalse(range.overlaps(QgsDateRange(QDate(2009, 4, 1), QDate(2009, 8, 5))))
        self.assertFalse(range.overlaps(QgsDateRange(QDate(2019, 4, 1), QDate(2019, 8, 5))))

        range = QgsDateRange(QDate(), QDate(2010, 6, 2))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2009, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2017, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate())))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(), QDate(2010, 4, 1))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2009, 4, 1), QDate(2009, 8, 5))))
        self.assertFalse(range.overlaps(QgsDateRange(QDate(2019, 4, 1), QDate(2019, 8, 5))))

        range = QgsDateRange(QDate(2010, 3, 1), QDate())
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2009, 4, 1), QDate(2010, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate(2017, 4, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2010, 4, 1), QDate())))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(), QDate(2010, 4, 1))))
        self.assertFalse(range.overlaps(QgsDateRange(QDate(2009, 4, 1), QDate(2009, 8, 5))))
        self.assertTrue(range.overlaps(QgsDateRange(QDate(2019, 4, 1), QDate(2019, 8, 5))))

    def testIsInstant(self):
        self.assertFalse(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2)).isInstant())
        self.assertTrue(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 3, 1)).isInstant())
        self.assertFalse(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 3, 1), False, False).isInstant())
        self.assertFalse(QgsDateRange(QDate(), QDate()).isInstant())

    def testIsInfinite(self):
        self.assertFalse(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2)).isInfinite())
        self.assertFalse(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 3, 1)).isInfinite())
        self.assertFalse(QgsDateRange(QDate(2010, 3, 1), QDate(2010, 3, 1), False, False).isInfinite())
        self.assertTrue(QgsDateRange(QDate(), QDate()).isInfinite())

    def testEquality(self):
        range = QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2), False, False)
        self.assertEqual(range, QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2), False, False))
        self.assertNotEqual(range, QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2), False, True))
        self.assertNotEqual(range, QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 2), True, False))
        self.assertNotEqual(range, QgsDateRange(QDate(2010, 3, 1), QDate(2010, 6, 3), False, False))
        self.assertNotEqual(range, QgsDateRange(QDate(2010, 3, 2), QDate(2010, 6, 2), False, False))


if __name__ == "__main__":
    unittest.main()
