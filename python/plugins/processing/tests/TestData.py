# -*- coding: utf-8 -*-

"""
***************************************************************************
    TestData.py
    ---------------------
    Date                 : March 2013
    Copyright            : (C) 2013 by Victor Olaya
    Email                : volayaf at gmail dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Victor Olaya'
__date__ = 'March 2013'
__copyright__ = '(C) 2013, Victor Olaya'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import os.path

testDataPath = os.path.join(os.path.dirname(__file__), 'testdata')


def table():
    return os.path.join(testDataPath, 'table.dbf')


def points():
    return os.path.join(testDataPath, 'points.gml')


def invalid_geometries():
    return os.path.join(testDataPath, 'invalidgeometries.gml')
