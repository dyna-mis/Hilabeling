# -*- coding: utf-8 -*-
"""QGIS Unit tests for QgsGeometryValidator.

.. note:: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
"""
__author__ = 'Nyall Dawson'
__date__ = '03/10/2016'
__copyright__ = 'Copyright 2016, The QGIS Project'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '176c06ceefb5f555205e72b20c962740cc0ec183'

from qgis.core import (
    QgsGeometry,
    QgsGeometryValidator
)

from qgis.testing import (
    unittest
)


class TestQgsGeometryValidator(unittest.TestCase):

    def testIssue15660(self):
        """ Test crash when validating geometry (#15660) """
        g = QgsGeometry.fromWkt(
            "Polygon ((0.44256348235389709 -47.87645625696347906, -2.88231630340906797 -47.90003919913998232,"
            "-2.88589842578005751 -48.91215450743293047, -2.8858984257800584 -48.91215450743293047,"
            "-2.88589842578005751 -48.91215450743292337, -2.88589842632776516 -48.9121545074024624,"
            "-3.24858148608664266 -48.89198543875494352, -3.27689849271356159 -49.40119850743292318,"
            "-4.37689842578006072 -49.40119850743292318, -4.3768984257800625 -49.40119850743293739,"
            "-6.11689842578005738 -49.40119850743293739, -6.11689842578006093 -49.40119850743292318,"
            "-7.03689842578006086 -49.40119850743292318, -7.02239489567359776 -48.93302358132471852,"
            "-7.02177188091450688 -48.91291272055079276, -7.02177188393206286 -48.91291272045731375,"
            "-7.02141642578006309 -48.9014385074329212, -7.7002102410998674 -48.88041051649613422,"
            "-7.70077301577442341 -48.89187793078160382, -7.70077301577442341 -48.89187793078160382,"
            "-7.70233865095334291 -48.92378019651650334, -7.72576842578006051 -49.40119850743292318,"
            "-9.26576842578005966 -49.40119850743292318, -9.26576842578006321 -49.40119850743293739,"
            "-13.28076842578006023 -49.40119850743293739, -13.04700849136197149 -44.82717853953759857,"
            "-12.22739168108193297 -44.85224022031006541, -12.22501286465108805 -44.774446133668377,"
            "-12.22288921611744783 -44.774511069530881, -12.2155540445085915 -44.53462318893357264,"
            "-13.0310217329353506 -44.50968822589503304, -12.87640859053235687 -41.29089836691012749,"
            "-7.83390711693117936 -41.74840291007100745, -7.88211379129075596 -42.99075321817508666,"
            "-7.43245210877673657 -43.00820115628285123, -7.50410812906098013 -44.67868742523263847,"
            "-7.52086717830689011 -44.67817497540850979, -7.52820234991574644 -44.91806285600581816,"
            "-7.51439432253991058 -44.91848507095944854, -7.57421591980290287 -46.3130804785479242,"
            "-8.32385639731993621 -46.28985691678211367, -8.44985043007881842 -48.85718773355701217,"
            "-6.26478736265246283 -48.92487774800262912, -6.18500945357052245 -46.35611749220839073,"
            "-6.94163842159974198 -46.33267751510010157, -6.82382190915497944 -40.77447960826383166,"
            "-5.48772831582523146 -40.77837853452808758, -5.48973219643841759 -39.92687558952010107,"
            "-2.75658441771447116 -40.04490036239724304, -3.1241861109063862 -46.6551270968877958,"
            "-2.78977790434136885 -46.6737244011090695, -2.78977790434136796 -46.6737244011090695,"
            "-2.71083842578005996 -44.83541850743291945, -2.71083842578005729 -44.83541850743291945,"
            "-0.86779302740823816 -44.89143693883772812, -0.86745855610774569 -44.87743669555854353,"
            "0.29843811058281489 -44.90401226269042922, 0.20437651721061911 -46.69301920907949466,"
            "0.50389019278376956 -46.71008040893148916, 0.44256348235389709 -47.87645625696347906))")

        self.assertTrue(g)
        # make sure validating this geometry doesn't crash QGIS
        QgsGeometryValidator.validateGeometry(g)


if __name__ == '__main__':
    unittest.main()
