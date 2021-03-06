# -*- coding: utf-8 -*-

"""
***************************************************************************
    PostGISExecuteAndLoadSQL.py
    ---------------------
    Date                 : May 2018
    Copyright            : (C) 2018 by Anita Graser
    Email                : anitagraser at gmx dot at
    ---------------------
    based on PostGISExecuteSQL.py by  Victor Olaya and Carterix Geomatics
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Anita Graser'
__date__ = 'May 2018'
__copyright__ = '(C) 2018, Anita Graser'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

from qgis.core import (Qgis,
                       QgsProcessingException,
                       QgsProcessingParameterString,
                       QgsApplication,
                       QgsVectorLayer,
                       QgsProject,
                       QgsProcessing,
                       QgsProcessingException,
                       QgsProcessingOutputVectorLayer,
                       QgsProcessingContext,
                       QgsProcessingFeedback)
from processing.algs.qgis.QgisAlgorithm import QgisAlgorithm
from processing.tools import postgis


class PostGISExecuteAndLoadSQL(QgisAlgorithm):

    DATABASE = 'DATABASE'
    SQL = 'SQL'
    OUTPUT = 'OUTPUT'
    ID_FIELD = 'ID_FIELD'
    GEOMETRY_FIELD = 'GEOMETRY_FIELD'

    def group(self):
        return self.tr('Database')

    def groupId(self):
        return 'database'

    def __init__(self):
        super().__init__()

    def initAlgorithm(self, config=None):
        db_param = QgsProcessingParameterString(
            self.DATABASE,
            self.tr('Database (connection name)'))
        db_param.setMetadata({
            'widget_wrapper': {
                'class': 'processing.gui.wrappers_postgis.ConnectionWidgetWrapper'}})
        self.addParameter(db_param)
        self.addParameter(QgsProcessingParameterString(
            self.SQL,
            self.tr('SQL query'),
            multiLine=True))
        self.addParameter(QgsProcessingParameterString(
            self.ID_FIELD,
            self.tr('Unique ID field name'),
            defaultValue='id'))
        self.addParameter(QgsProcessingParameterString(
            self.GEOMETRY_FIELD,
            self.tr('Geometry field name'),
            defaultValue='geom',
            optional=True))
        self.addOutput(QgsProcessingOutputVectorLayer(
            self.OUTPUT,
            self.tr("Output layer"),
            QgsProcessing.TypeVectorAnyGeometry))

    def name(self):
        return 'postgisexecuteandloadsql'

    def displayName(self):
        return self.tr('PostgreSQL execute and load SQL')

    def shortDescription(self):
        return self.tr('Executes a SQL command on a PostgreSQL database and loads the result as a table')

    def tags(self):
        return self.tr('postgis,table,database').split(',')

    def processAlgorithm(self, parameters, context, feedback):
        connection = self.parameterAsString(parameters, self.DATABASE, context)
        id_field = self.parameterAsString(parameters, self.ID_FIELD, context)
        geom_field = self.parameterAsString(
            parameters, self.GEOMETRY_FIELD, context)
        uri = postgis.uri_from_name(connection)
        sql = self.parameterAsString(parameters, self.SQL, context)
        sql = sql.replace('\n', ' ')
        uri.setDataSource("", "(" + sql + ")", geom_field, "", id_field)

        vlayer = QgsVectorLayer(uri.uri(), "layername", "postgres")

        if not vlayer.isValid():
            raise QgsProcessingException(self.tr("""This layer is invalid!
                Please check the PostGIS log for error messages."""))

        context.temporaryLayerStore().addMapLayer(vlayer)
        context.addLayerToLoadOnCompletion(
            vlayer.id(),
            QgsProcessingContext.LayerDetails('SQL layer',
                                              context.project(),
                                              self.OUTPUT))

        return {self.OUTPUT: vlayer.id()}
