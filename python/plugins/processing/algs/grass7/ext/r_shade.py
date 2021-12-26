# -*- coding: utf-8 -*-

"""
***************************************************************************
    r_shade.py
    ----------
    Date                 : February 2016
    Copyright            : (C) 2016 by Médéric Ribreux
    Email                : medspx at medspx dot fr
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Médéric Ribreux'
__date__ = 'February 2016'
__copyright__ = '(C) 2016, Médéric Ribreux'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'


def processInputs(alg, parameters, context, feedback):
    # We need to import all the bands and color tables of the input rasters
    alg.loadRasterLayerFromParameter('shade', parameters, context,
                                     False, None)
    alg.loadRasterLayerFromParameter('color', parameters, context,
                                     False, None)


def processOutputs(alg, parameters, context, feedback):
    # Keep color table
    alg.exportRasterLayerFromParameter('output', parameters, context, True)
