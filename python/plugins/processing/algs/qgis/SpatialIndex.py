# -*- coding: utf-8 -*-

"""
***************************************************************************
    SpatialIndex.py
    ---------------------
    Date                 : January 2016
    Copyright            : (C) 2016 by Alexander Bruy
    Email                : alexander dot bruy at gmail dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Alexander Bruy'
__date__ = 'January 2016'
__copyright__ = '(C) 2016, Alexander Bruy'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

import os

from qgis.core import (QgsVectorDataProvider,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterVectorLayer,
                       QgsProcessingOutputVectorLayer,
                       QgsProcessing)

from processing.algs.qgis.QgisAlgorithm import QgisAlgorithm

pluginPath = os.path.split(os.path.split(os.path.dirname(__file__))[0])[0]


class SpatialIndex(QgisAlgorithm):

    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'

    def group(self):
        return self.tr('Vector general')

    def groupId(self):
        return 'vectorgeneral'

    def __init__(self):
        super().__init__()

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterVectorLayer(self.INPUT,
                                                            self.tr('Input Layer'),
                                                            [QgsProcessing.TypeVectorPolygon, QgsProcessing.TypeVectorPoint, QgsProcessing.TypeVectorLine]))
        self.addOutput(QgsProcessingOutputVectorLayer(self.OUTPUT, self.tr('Indexed layer')))

    def flags(self):
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def name(self):
        return 'createspatialindex'

    def displayName(self):
        return self.tr('Create spatial index')

    def processAlgorithm(self, parameters, context, feedback):
        layer = self.parameterAsVectorLayer(parameters, self.INPUT, context)
        provider = layer.dataProvider()

        if provider.capabilities() & QgsVectorDataProvider.CreateSpatialIndex:
            if not provider.createSpatialIndex():
                feedback.pushInfo(self.tr('Could not create spatial index'))
        else:
            feedback.pushInfo(self.tr("Layer's data provider does not support "
                                      "spatial indexes"))

        return {self.OUTPUT: layer.id()}
