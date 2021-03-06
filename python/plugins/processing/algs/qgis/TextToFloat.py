# -*- coding: utf-8 -*-

"""
***************************************************************************
    TextToFloat.py
    ---------------------
    Date                 : May 2010
    Copyright            : (C) 2010 by Michael Minn
    Email                : pyqgis at michaelminn dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Michael Minn'
__date__ = 'May 2010'
__copyright__ = '(C) 2010, Michael Minn'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

from qgis.PyQt.QtCore import QVariant
from qgis.core import (QgsField,
                       QgsProcessing,
                       QgsProcessingParameterField,
                       QgsProcessingFeatureSource)
from processing.algs.qgis.QgisAlgorithm import QgisFeatureBasedAlgorithm


class TextToFloat(QgisFeatureBasedAlgorithm):

    FIELD = 'FIELD'

    def group(self):
        return self.tr('Vector table')

    def groupId(self):
        return 'vectortable'

    def __init__(self):
        super().__init__()
        self.field_name = None
        self.field_idx = -1

    def initParameters(self, config=None):
        self.addParameter(QgsProcessingParameterField(self.FIELD,
                                                      self.tr('Text attribute to convert to float'),
                                                      parentLayerParameterName='INPUT',
                                                      type=QgsProcessingParameterField.String
                                                      ))

    def name(self):
        return 'texttofloat'

    def displayName(self):
        return self.tr('Text to float')

    def outputName(self):
        return self.tr('Float from text')

    def inputLayerTypes(self):
        return [QgsProcessing.TypeVector]

    def outputFields(self, inputFields):
        self.field_idx = inputFields.lookupField(self.field_name)
        if self.field_idx >= 0:
            inputFields[self.field_idx] = QgsField(self.field_name, QVariant.Double, '', 24, 15)
        return inputFields

    def prepareAlgorithm(self, parameters, context, feedback):
        self.field_name = self.parameterAsString(parameters, self.FIELD, context)
        return True

    def supportInPlaceEdit(self, layer):
        return False

    def sourceFlags(self):
        return QgsProcessingFeatureSource.FlagSkipGeometryValidityChecks

    def processFeature(self, feature, context, feedback):
        value = feature[self.field_idx]
        try:
            if '%' in value:
                feature[self.field_idx] = float(value.replace('%', '')) / 100.0
            else:
                feature[self.field_idx] = float(value)
        except:
            feature[self.field_idx] = None
        return [feature]
