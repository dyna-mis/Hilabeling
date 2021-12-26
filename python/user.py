# -*- coding: utf-8 -*-

"""
***************************************************************************
    user.py
    ---------------------
    Date                 : January 2015
    Copyright            : (C) 2015 by Nathan Woodrow
    Email                : woodrow dot nathan at gmail dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Nathan Woodrow'
__date__ = 'January 2015'
__copyright__ = '(C) 2015, Nathan Woodrow'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import os
import sys
import glob
import traceback

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import Qgis, QgsApplication, QgsMessageLog


def load_user_expressions(path):
    """
    Load all user expressions from the given paths
    """
    #Loop all py files and import them
    modules = glob.glob(path + "/*.py")
    names = [os.path.basename(f)[:-3] for f in modules]
    for name in names:
        if name == "__init__":
            continue
        # As user expression functions should be registered with qgsfunction
        # just importing the file is enough to get it to load the functions into QGIS
        try:
            __import__("expressions.{0}".format(name), locals(), globals())
        except:
            error = traceback.format_exc()
            msgtitle = QCoreApplication.translate("UserExpressions", "User expressions")
            msg = QCoreApplication.translate("UserExpressions", "The user expression {0} is not valid").format(name)
            QgsMessageLog.logMessage(msg + "\n" + error, msgtitle, Qgis.Warning)


userpythonhome = os.path.join(QgsApplication.qgisSettingsDirPath(), "python")
expressionspath = os.path.join(userpythonhome, "expressions")

sys.path.append(userpythonhome)

if not os.path.exists(expressionspath):
    os.makedirs(expressionspath)

initfile = os.path.join(expressionspath, "__init__.py")
if not os.path.exists(initfile):
    open(initfile, "w").close()

template = """from qgis.core import *
from qgis.gui import *

@qgsfunction(args='auto', group='Custom')
def my_sum(value1, value2, feature, parent):
    \"\"\"
    Calculates the sum of the two parameters value1 and value2.
    <h2>Example usage:</h2>
    <ul>
      <li>my_sum(5, 8) -> 13</li>
      <li>my_sum(\"field1\", \"field2\") -> 42</li>
    </ul>
    \"\"\"
    return value1 + value2
"""

default_expression_template = template


try:
    import expressions

    expressions.load = load_user_expressions
    expressions.load(expressionspath)
    expressions.template = template
except ImportError:
    # We get a import error and crash for some reason even if we make the expressions package
    # TODO Fix the crash on first load with no expressions folder
    # But for now it's not the end of the world if it doesn't load the first time
    pass
