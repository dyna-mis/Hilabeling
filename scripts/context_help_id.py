#!/usr/bin/env python3
"""
/***************************************************************************
                               context_help_id.py
                              -------------------
    begin                : 2009-11-16
    copyright            : (C) 2009 by Gary E.Sherman
    email                : gsherman at mrcc.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

    This script generates a unique context id based for use in the QGIS
    context sensitive help system. It uses the SHA1 hash for the class name
    and converts the first 12 characters to a unique integer.

    To create a context id, pass the name of the QGIS class on the command line.
    Example:
        ./context_help_id.py QgsAbout

    This script requires Python 2.5 or higher (hashlib was introduced at 2.5).

    NOTE: Due to a change in the way context ids are generated, ids
    generated by the old method (Java hashCode function) will be different than
    the id generated by the new method for the same class.
"""
import hashlib
import sys
# check to see if a class name was specified and if so, create the context id
if len(sys.argv) > 1:
    hash = hashlib.sha1()
    # set the hash to the name passed on the command line
    hash.update(sys.argv[1])
    # generate the context id by converting the first 12 characters of the hash
    # to decimal
    context_id = int(hash.hexdigest()[:12], 16)
    # print the result
    print context_id
else:
    # if no class name was specified, give a bit of help
    print "To generate a context sensitive help id, specify the QGIS class name on the command line"
