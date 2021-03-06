# -*- coding: utf-8 -*-

"""
***************************************************************************
    ogr2ogrtabletopostgislist.py
    ---------------------
    Date                 : November 2012
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
__date__ = 'November 2012'
__copyright__ = '(C) 2012, Victor Olaya'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

from qgis.core import QgsSettings

from processing.algs.gdal.GdalAlgorithm import GdalAlgorithm
from processing.algs.gdal.GdalUtils import GdalUtils

from processing.tools.postgis import uri_from_name, GeoDB
from processing.tools.system import isWindows


class Ogr2OgrTableToPostGisList(GdalAlgorithm):

    DATABASE = 'DATABASE'
    INPUT_LAYER = 'INPUT_LAYER'
    SHAPE_ENCODING = 'SHAPE_ENCODING'
    HOST = 'HOST'
    PORT = 'PORT'
    USER = 'USER'
    DBNAME = 'DBNAME'
    PASSWORD = 'PASSWORD'
    SCHEMA = 'SCHEMA'
    TABLE = 'TABLE'
    PK = 'PK'
    PRIMARY_KEY = 'PRIMARY_KEY'
    WHERE = 'WHERE'
    GT = 'GT'
    OVERWRITE = 'OVERWRITE'
    APPEND = 'APPEND'
    ADDFIELDS = 'ADDFIELDS'
    LAUNDER = 'LAUNDER'
    SKIPFAILURES = 'SKIPFAILURES'
    PRECISION = 'PRECISION'
    OPTIONS = 'OPTIONS'

    def __init__(self):
        GdalAlgorithm.__init__(self)

    def dbConnectionNames(self):
        settings = QgsSettings()
        settings.beginGroup('/PostgreSQL/connections/')
        return settings.childGroups()

    def __init__(self):
        super().__init__()

    def initAlgorithm(self, config=None):
        self.DB_CONNECTIONS = self.dbConnectionNames()
        self.addParameter(ParameterSelection(self.DATABASE,
                                             self.tr('Database (connection name)'), self.DB_CONNECTIONS))
        self.addParameter(ParameterTable(self.INPUT_LAYER,
                                         self.tr('Input layer')))
        self.addParameter(ParameterString(self.SHAPE_ENCODING,
                                          self.tr('Shape encoding'), "", optional=True))
        self.addParameter(ParameterString(self.SCHEMA,
                                          self.tr('Schema name'), 'public', optional=True))
        self.addParameter(ParameterString(self.TABLE,
                                          self.tr('Table name, leave blank to use input name'),
                                          '', optional=True))
        self.addParameter(ParameterString(self.PK,
                                          self.tr('Primary key'), 'id', optional=True))
        self.addParameter(ParameterTableField(self.PRIMARY_KEY,
                                              self.tr('Primary key (existing field, used if the above option is left empty)'),
                                              self.INPUT_LAYER, optional=True))
        self.addParameter(ParameterString(self.WHERE,
                                          self.tr('Select features using a SQL "WHERE" statement (Ex: column=\'value\')'),
                                          '', optional=True))
        self.addParameter(ParameterString(self.GT,
                                          self.tr('Group N features per transaction (Default: 20000)'),
                                          '', optional=True))
        self.addParameter(ParameterBoolean(self.OVERWRITE,
                                           self.tr('Overwrite existing table'), True))
        self.addParameter(ParameterBoolean(self.APPEND,
                                           self.tr('Append to existing table'), False))
        self.addParameter(ParameterBoolean(self.ADDFIELDS,
                                           self.tr('Append and add new fields to existing table'), False))
        self.addParameter(ParameterBoolean(self.LAUNDER,
                                           self.tr('Do not launder columns/table names'), False))
        self.addParameter(ParameterBoolean(self.SKIPFAILURES,
                                           self.tr('Continue after a failure, skipping the failed record'),
                                           False))
        self.addParameter(ParameterBoolean(self.PRECISION,
                                           self.tr('Keep width and precision of input attributes'),
                                           True))
        self.addParameter(ParameterString(self.OPTIONS,
                                          self.tr('Additional creation options'), '', optional=True))

    def name(self):
        return 'importlayertableasgeometrylesstableintopostgresqldatabase'

    def displayName(self):
        return self.tr('Import layer/table as geometryless table into PostgreSQL database')

    def group(self):
        return self.tr('Vector miscellaneous')

    def groupId(self):
        return 'vectormiscellaneous'

    def getConsoleCommands(self, parameters, context, feedback, executing=True):
        connection = self.DB_CONNECTIONS[self.getParameterValue(self.DATABASE)]
        uri = uri_from_name(connection)
        if executing:
            # to get credentials input when needed
            uri = GeoDB(uri=uri).uri

        inLayer = self.getParameterValue(self.INPUT_LAYER)
        ogrLayer = GdalUtils.ogrConnectionStringFromLayer(inLayer)
        shapeEncoding = self.getParameterValue(self.SHAPE_ENCODING)
        schema = str(self.getParameterValue(self.SCHEMA))
        table = str(self.getParameterValue(self.TABLE))
        pk = str(self.getParameterValue(self.PK))
        pkstring = "-lco FID=" + pk
        primary_key = self.getParameterValue(self.PRIMARY_KEY)
        where = str(self.getParameterValue(self.WHERE))
        wherestring = '-where "' + where + '"'
        gt = str(self.getParameterValue(self.GT))
        overwrite = self.getParameterValue(self.OVERWRITE)
        append = self.getParameterValue(self.APPEND)
        addfields = self.getParameterValue(self.ADDFIELDS)
        launder = self.getParameterValue(self.LAUNDER)
        launderstring = "-lco LAUNDER=NO"
        skipfailures = self.getParameterValue(self.SKIPFAILURES)
        precision = self.getParameterValue(self.PRECISION)
        options = str(self.getParameterValue(self.OPTIONS))

        arguments = []
        arguments.append('-progress')
        arguments.append('--config PG_USE_COPY YES')
        if len(shapeEncoding) > 0:
            arguments.append('--config')
            arguments.append('SHAPE_ENCODING')
            arguments.append('"' + shapeEncoding + '"')
        arguments.append('-f')
        arguments.append('PostgreSQL')
        arguments.append('PG:"')
        for token in uri.connectionInfo(executing).split(' '):
            arguments.append(token)
        arguments.append('active_schema={}'.format(schema or 'public'))
        arguments.append('"')
        arguments.append(ogrLayer)
        arguments.append('-nlt NONE')
        arguments.append(GdalUtils.ogrLayerName(inLayer))
        if launder:
            arguments.append(launderstring)
        if append:
            arguments.append('-append')
        if addfields:
            arguments.append('-addfields')
        if overwrite:
            arguments.append('-overwrite')
        if len(pk) > 0:
            arguments.append(pkstring)
        elif primary_key is not None:
            arguments.append("-lco FID=" + primary_key)
        if len(table) == 0:
            table = GdalUtils.ogrLayerName(inLayer).lower()
        if schema:
            table = '{}.{}'.format(schema, table)
        arguments.append('-nln')
        arguments.append(table)
        if skipfailures:
            arguments.append('-skipfailures')
        if where:
            arguments.append(wherestring)
        if len(gt) > 0:
            arguments.append('-gt')
            arguments.append(gt)
        if not precision:
            arguments.append('-lco PRECISION=NO')
        if len(options) > 0:
            arguments.append(options)

        commands = []
        if isWindows():
            commands = ['cmd.exe', '/C ', 'ogr2ogr.exe',
                        GdalUtils.escapeAndJoin(arguments)]
        else:
            commands = ['ogr2ogr', GdalUtils.escapeAndJoin(arguments)]

        return commands

    def commandName(self):
        return "ogr2ogr"
