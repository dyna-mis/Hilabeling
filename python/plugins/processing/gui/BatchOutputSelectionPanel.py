# -*- coding: utf-8 -*-

"""
***************************************************************************
    BatchOutputSelectionPanel.py
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
import re

from qgis.core import (QgsMapLayer,
                       QgsSettings,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterVectorLayer,
                       QgsProcessingParameterMultipleLayers,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterMatrix)
from qgis.PyQt.QtWidgets import QWidget, QPushButton, QLineEdit, QHBoxLayout, QSizePolicy, QFileDialog

from processing.gui.AutofillDialog import AutofillDialog
from processing.gui.ParameterGuiUtils import getFileFilter


class BatchOutputSelectionPanel(QWidget):

    def __init__(self, output, alg, row, col, panel):
        super(BatchOutputSelectionPanel, self).__init__(None)
        self.alg = alg
        self.row = row
        self.col = col
        self.output = output
        self.panel = panel
        self.table = self.panel.tblParameters
        self.horizontalLayout = QHBoxLayout(self)
        self.horizontalLayout.setSpacing(2)
        self.horizontalLayout.setMargin(0)
        self.text = QLineEdit()
        self.text.setText('')
        self.text.setMinimumWidth(300)
        self.text.setSizePolicy(QSizePolicy.Expanding,
                                QSizePolicy.Expanding)
        self.horizontalLayout.addWidget(self.text)
        self.pushButton = QPushButton()
        self.pushButton.setText('???')
        self.pushButton.clicked.connect(self.showSelectionDialog)
        self.horizontalLayout.addWidget(self.pushButton)
        self.setLayout(self.horizontalLayout)

    def showSelectionDialog(self):
        if isinstance(self.output, QgsProcessingParameterFolderDestination):
            self.selectDirectory()
            return

        filefilter = getFileFilter(self.output)
        settings = QgsSettings()
        if settings.contains('/Processing/LastBatchOutputPath'):
            path = str(settings.value('/Processing/LastBatchOutputPath'))
        else:
            path = ''
        filename, selectedFileFilter = QFileDialog.getSaveFileName(self,
                                                                   self.tr('Save File'), path, filefilter)
        if filename:
            if not filename.lower().endswith(
                    tuple(re.findall("\\*(\\.[a-z]{1,10})", filefilter))):
                ext = re.search("\\*(\\.[a-z]{1,10})", selectedFileFilter)
                if ext:
                    filename += ext.group(1)
            settings.setValue('/Processing/LastBatchOutputPath', os.path.dirname(filename))
            dlg = AutofillDialog(self.alg)
            dlg.exec_()
            if dlg.mode is not None:
                try:
                    if dlg.mode == AutofillDialog.DO_NOT_AUTOFILL:
                        self.table.cellWidget(self.row,
                                              self.col).setValue(filename)
                    elif dlg.mode == AutofillDialog.FILL_WITH_NUMBERS:
                        n = self.table.rowCount() - self.row
                        for i in range(n):
                            name = filename[:filename.rfind('.')] \
                                + str(i + 1) + filename[filename.rfind('.'):]
                            self.table.cellWidget(i + self.row,
                                                  self.col).setValue(name)
                    elif dlg.mode == AutofillDialog.FILL_WITH_PARAMETER:
                        n = self.table.rowCount() - self.row
                        for i in range(n):
                            widget = self.table.cellWidget(i + self.row,
                                                           dlg.param_index)
                            param = self.alg.parameterDefinitions()[dlg.param_index]
                            if isinstance(param, (QgsProcessingParameterRasterLayer,
                                                  QgsProcessingParameterFeatureSource,
                                                  QgsProcessingParameterVectorLayer,
                                                  QgsProcessingParameterMultipleLayers)):
                                v = widget.value()
                                if isinstance(v, QgsMapLayer):
                                    s = v.name()
                                else:
                                    s = os.path.basename(v)
                                    s = os.path.splitext(s)[0]
                            elif isinstance(param, QgsProcessingParameterBoolean):
                                s = str(widget.currentIndex() == 0)
                            elif isinstance(param, QgsProcessingParameterEnum):
                                s = str(widget.currentText())
                            elif isinstance(param, QgsProcessingParameterMatrix):
                                s = str(widget.table)
                            else:
                                s = str(widget.text())
                            name = filename[:filename.rfind('.')] + s \
                                + filename[filename.rfind('.'):]
                            self.table.cellWidget(i + self.row,
                                                  self.col).setValue(name)
                except:
                    pass

    def selectDirectory(self):

        settings = QgsSettings()
        if settings.contains('/Processing/LastBatchOutputPath'):
            lastDir = str(settings.value('/Processing/LastBatchOutputPath'))
        else:
            lastDir = ''

        dirName = QFileDialog.getExistingDirectory(self,
                                                   self.tr('Output Directory'), lastDir, QFileDialog.ShowDirsOnly)

        if dirName:
            self.table.cellWidget(self.row, self.col).setValue(dirName)
            settings.setValue('/Processing/LastBatchOutputPath', dirName)

    def setValue(self, text):
        return self.text.setText(text)

    def getValue(self):
        return str(self.text.text())
