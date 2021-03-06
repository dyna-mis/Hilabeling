# -*- coding: utf-8 -*-

"""
***************************************************************************
    __init__.py
    ---------------------
    Date                 : August 2012
    Copyright            : (C) 2012 by Victor Olaya
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
__date__ = 'August 2012'
__copyright__ = '(C) 2012, Victor Olaya'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

from processing.tools.dataobjects import *          # NOQA
from processing.tools.general import *              # NOQA
from processing.tools.vector import *               # NOQA
from processing.tools.raster import *               # NOQA
from processing.tools.system import *               # NOQA

# monkey patch Python specific Processing API into stable qgis.processing module
import qgis.processing
qgis.processing.algorithmHelp = algorithmHelp
qgis.processing.run = run
qgis.processing.runAndLoadResults = runAndLoadResults
qgis.processing.createAlgorithmDialog = createAlgorithmDialog
qgis.processing.execAlgorithmDialog = execAlgorithmDialog
qgis.processing.createContext = createContext


def classFactory(iface):
    from processing.ProcessingPlugin import ProcessingPlugin
    return ProcessingPlugin(iface)
