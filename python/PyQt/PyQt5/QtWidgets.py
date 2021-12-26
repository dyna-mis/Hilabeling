# -*- coding: utf-8 -*-

"""
***************************************************************************
    QtWidgets.py
    ---------------------
    Date                 : November 2015
    Copyright            : (C) 2015 by Matthias Kuhn
    Email                : matthias at opengis dot ch
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Matthias Kuhn'
__date__ = 'November 2015'
__copyright__ = '(C) 2015, Matthias Kuhn'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

from PyQt5.QtWidgets import *

QLayout.setMargin = lambda self, m: self.setContentsMargins(m, m, m, m)
