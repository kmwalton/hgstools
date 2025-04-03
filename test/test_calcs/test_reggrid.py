#!/usr/bin/env python


import os
import unittest
import numpy as np

import warnings

from pyhgs.mesh import HGSGrid, Domain
from pyhgs.calcs import AvRegGrid

from pyhgs.test import sims_join

import logging
logger = logging.getLogger('pyhgs.calcs')
perf_logger = logging.getLogger('pyhgs.calcs_perf')

class Test_Av(unittest.TestCase):

    @unittest.skipIf( not os.path.isfile(sims_join('04b_very_coarse_mesh','module4bo.eco')),
            'Could not find source simulation')
    def test_coarse_regridding(self):
        _d = sims_join('04b_very_coarse_mesh','module4b')
        g = HGSGrid(_d)
        calc = AvRegGrid(g, [25., 1., 25.,])
        self.assertTrue(np.allclose(calc.gl[0], np.array([0., 25., 50.])))
        self.assertTrue(np.allclose(calc.gl[1], np.array([0., 1.])))
        self.assertTrue(np.allclose(calc.gl[2], np.array([0., 25.,])))

        _pm_map = calc.e2g[Domain.PM]
        _fx_map = calc.e2g[Domain.FRAC]

        #breakpoint()
        with self.subTest('mapping dimensions'):
            self.assertEqual(_pm_map.shape, (20,2))
            self.assertEqual(_fx_map.shape, (11,2))

        with self.subTest('PM, spotcheck'):
            self.assertTrue(_pm_map[7,0] == 1.)
            self.assertTrue(_pm_map[7,1] == 0.)
            self.assertTrue(_pm_map[8,0] == 0.)
            self.assertTrue(_pm_map[8,1] == 1.)

        with self.subTest('Fx, spotcheck'):
            self.assertTrue(_fx_map[8,0] == 1.)
            self.assertTrue(_fx_map[8,1] == 0.)
            self.assertTrue(_fx_map[9,0] == 0.)
            self.assertTrue(_fx_map[9,1] == 1.)

        with self.subTest('boundary, is grouped'):
            self.assertTrue(_fx_map[6,0] > 0.)
            self.assertTrue(_fx_map[6,1] > 0.)

        with self.subTest('boundary, is proportioned'):
            self.assertTrue(_fx_map[6,0] == 0.5)
            self.assertTrue(_fx_map[6,1] == 0.5)
        
    @unittest.skipIf( not os.path.isfile(sims_join('04b_very_coarse_mesh','module4bo.coordinates_pm')),
            'Could not find source simulation')
    def test_overflow_regridding(self):

        _d = sims_join('04b_very_coarse_mesh','module4b')

        g = HGSGrid(_d)
        calc = AvRegGrid(g, [12.5, 1., 10.,])

        self.assertTrue(np.allclose(calc.gl[0], np.linspace(0.,50.,5)))
        self.assertTrue(np.allclose(calc.gl[1], np.array([0., 1.])))
        self.assertTrue(np.allclose(calc.gl[2], np.linspace(0.,30.,4)))

        _pm_map = calc.e2g[Domain.PM]
        _fx_map = calc.e2g[Domain.FRAC]
        with self.subTest('mapping dimensions'):
            self.assertEqual(_pm_map.shape, (20,12))
            self.assertEqual(_fx_map.shape, (11,12))

        
    @unittest.skipIf( not os.path.isfile(sims_join('04b_very_coarse_mesh','module4bo.coordinates_pm')),
            'Could not find source simulation')
    def test_medium_regridding(self):

        _d = sims_join('04b_very_coarse_mesh','module4b')
        g = HGSGrid(_d)
        calc = AvRegGrid(g, [5., 1., 5.,])

        self.assertTrue(np.allclose(calc.gl[0], np.linspace(0., 50., 11)))
        self.assertTrue(np.allclose(calc.gl[1], np.array([0., 1.])))
        self.assertTrue(np.allclose(calc.gl[2], np.linspace(0., 25., 6)))

        _pm_map = calc.e2g[Domain.PM]
        _fx_map = calc.e2g[Domain.FRAC]
        with self.subTest('mapping dimensions'):
            self.assertEqual(_pm_map.shape, (20,50))
            self.assertEqual(_fx_map.shape, (11,50))

        with self.subTest('rows spotcheck'):
            a = _pm_map[0,:].todense()
            b = np.zeros((1,50))
            b[0,0:2] = 25./60.
            b[0,10:12] = 5./60.
            np.testing.assert_allclose(a,b)

            a = _pm_map[11,:].todense()
            b = np.zeros((1,50), dtype=np.single)
            b[0,22:24] = 15./80.
            b[0,32:34] = 25./80.
            np.testing.assert_allclose(a,b)

        with self.subTest('columns spotcheck'):
            a = _pm_map[:,21].todense()
            b = np.zeros((20,1))
            b[5,0] = 10./60.
            b[10,0] = 15./80.
            np.testing.assert_allclose(a,b)

if __name__=='__main__':

    import sys
    import argparse
    argp = argparse.ArgumentParser()
    argp.add_argument('-v', '--verbose', action='count', default=0)
    args, unmatched_args = argp.parse_known_args(sys.argv)

    if args.verbose > 1:
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.DEBUG)
        perf_logger.addHandler(logging.StreamHandler())
        perf_logger.setLevel(logging.INFO)
    elif args.verbose > 0:
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)

    unittest.main()

