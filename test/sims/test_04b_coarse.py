#!/usr/bin/env python
"""Python unittests involving 04b_Saturated_Fracture_Transport"""

import unittest
import os
from pyhgs._test import skip_if_no_sim_output 
from pyhgs.mesh import (HGSGrid, Domain)

TESTP = os.path.dirname(__file__)
'Path to this testing directory'

SIM_PREFIX = os.path.join(TESTP,'04b_very_coarse_mesh','module4b')
'Simulation directory+prefix'

@unittest.skipIf(skip_if_no_sim_output(SIM_PREFIX,['o.eco',]), 'HGS output missing')
class Test_Module04bCoarse(unittest.TestCase):

    def setUp(self):
        # Loading and re-loading this each test case could get expensive.
        # Consider refactoring to a "static" module level variable that can be
        # treated as read-only as a balance between data loading overhead and
        # test independence; test cases that modify the HGSGrid object can
        # reload the data.
        self.g = HGSGrid(SIM_PREFIX)

    def test_grid_sizes(self):
        self.assertEqual(self.g.shape,(6,2,5))
        self.assertEqual(self.g.elshape,(5,1,4))

    def test_adjacency(self):

        actual = self.g.make_pm_to_fx_adjacency()
        self.assertEqual(14, len(actual), 'Number of PM with adjacent fx')

        # spot-check some nodes with fracture neighbours
        some_desired = {
            2:[4,8,],
            3:[4,9,],
            10:[0,],
            15:[0,],
        }
        for ipm,desired in some_desired.items():
            with self.subTest(pm_element=ipm):
                self.assertEqual(desired, actual[ipm])

        # spot-check some nodes with no neighbours
        for ipm in [0, 1, 5, 6, 14, 19,]:
            with self.subTest(pm_element=ipm):
                self.assertEqual([], actual[ipm])

        # try invalid
        with self.assertRaises(ValueError):
            actual[-1]
        with self.assertRaises(ValueError):
            actual[20]

if __name__ == '__main__':
    unittest.main()
