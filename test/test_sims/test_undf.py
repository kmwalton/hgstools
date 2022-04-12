#!/usr/bin/env python
"""Python unittests involving 04b_Saturated_Fracture_Transport"""

import unittest
import os
from itertools import product
from math import prod
import numpy as np
import numpy.testing as nptest

import pyhgs
from pyhgs.parser import parse
from pyhgs._test import skip_if_no_sim_output 
from pyhgs.mesh import (HGSGrid, Domain)

TESTP = os.path.dirname(__file__)
'Path to this testing directory'

SIM_PREFIX = os.path.join(TESTP,'undulating_coarse','undf')
'Simulation directory+prefix'

@unittest.skipIf(skip_if_no_sim_output(SIM_PREFIX,[
        'o.eco', 'o.q_pm.0001', 'o.v_pm.0001', 'o.v_frac.0001',
        'o.conc_pm.salt.0010',],),
        'HGS output missing')
class Test_Module04bCoarse(unittest.TestCase):

    def setUp(self):
        # Loading and re-loading this each test case could get expensive.
        # Consider refactoring to a "static" module level variable that can be
        # treated as read-only as a balance between data loading overhead and
        # test independence; test cases that modify the HGSGrid object can
        # reload the data.
        self.g = HGSGrid(SIM_PREFIX)

        # test HGS' pm 2 fracture node mapping
        # Note: pyhgs' fracture nodes are indexed starting at zero
        self.pmn2fxn = self.g.hgs_fx_nodes['link_pm2frac'] - self.g.nn

    def test_grid_sizes(self):
        self.assertEqual(self.g.shape,(6,5,2))
        self.assertEqual(self.g.elshape,(5,4,1))

        self.assertEqual(self.g.nn, 60)
        self.assertEqual(self.g.ne, 20)
        self.assertEqual(self.g.nfn, 24)
        self.assertEqual(self.g.nfe, 11)


    def test_fx_elem_volume(self):

        ap = 100e-6
        exp = ap * np.array(
             [ 10., 12.5, 7.5, 6.25, 
                9.,  9., 12.,7.5, 
               6.25, 20., 7.5, ]
               # how interesting: the tallest element comes last in the sequence
            )

        fx_vol = self.g.get_element_volumes(Domain.FRAC)
        nptest.assert_allclose(fx_vol, exp)

    def test_pm_elem_volume(self):
        exp = np.array(
              [ [ 60.,  75.,  45.,  37.5, 120., ], # y- bottom row
                [ 60.,  75.,  45.,  37.5, 120., ],
                [ 80., 100.,  60.,   50., 160., ],
                [ 50., 62.5, 37.5, 31.25, 100., ], ]  # top row
            ).T
            # Transpose to match HGS; Easier to see x-rows horizontal on
            #when creating the test problem looking at plan view

        pm_vol = self.g.get_element_volumes()

        nptest.assert_allclose(pm_vol[:,:,0], exp)


if __name__ == '__main__':
    unittest.main()
