# -*- coding: utf-8 -*-

"""
***************************************************************************
    CreateNewScriptAction.py
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

import os

from qgis.PyQt.QtCore import QCoreApplication

from qgis.core import QgsApplication

from processing.gui.ToolboxAction import ToolboxAction

from processing.script.ScriptEditorDialog import ScriptEditorDialog


class CreateNewScriptAction(ToolboxAction):

    def __init__(self):
        self.name = QCoreApplication.translate("CreateNewScriptAction", "Create New Script…")
        self.group = self.tr("Tools")

    def execute(self):
        dlg = ScriptEditorDialog(None)
        dlg.show()
