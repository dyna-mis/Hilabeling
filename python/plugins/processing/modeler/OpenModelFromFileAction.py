# -*- coding: utf-8 -*-

"""
***************************************************************************
    OpenModelFromFileAction.py
    ---------------------
    Date                 : February 2018
    Copyright            : (C) 2018 by Nyall Dawson
    Email                : nyall dot dawson at gmail dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Nyall Dawson'
__date__ = 'February 2018'
__copyright__ = '(C) 2018, Nyall Dawson'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import os
from qgis.PyQt.QtWidgets import QFileDialog
from qgis.PyQt.QtCore import QFileInfo, QCoreApplication

from qgis.core import QgsApplication, QgsSettings

from processing.gui.ToolboxAction import ToolboxAction
from processing.modeler.ModelerDialog import ModelerDialog

pluginPath = os.path.split(os.path.dirname(__file__))[0]


class OpenModelFromFileAction(ToolboxAction):

    def __init__(self):
        self.name = QCoreApplication.translate('OpenModelFromFileAction', 'Open Existing Model…')
        self.group = self.tr('Tools')

    def getIcon(self):
        return QgsApplication.getThemeIcon("/processingModel.svg")

    def execute(self):
        settings = QgsSettings()
        lastDir = settings.value('Processing/lastModelsDir', '')
        filename, selected_filter = QFileDialog.getOpenFileName(self.toolbox,
                                                                self.tr('Open Model', 'AddModelFromFileAction'), lastDir,
                                                                self.tr('Processing models (*.model3 *.MODEL3)', 'AddModelFromFileAction'))
        if filename:
            settings.setValue('Processing/lastModelsDir',
                              QFileInfo(filename).absoluteDir().absolutePath())

            dlg = ModelerDialog()
            dlg.loadModel(filename)
            dlg.show()
