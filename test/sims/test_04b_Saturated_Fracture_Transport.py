#!/usr/bin/env python
"""Python unittests involving 04b_Saturated_Fracture_Transport
"""

import unittest
import os

from pyhgs.mesh import (HGSGrid, Domain)

TESTP = os.path.dirname(__file__)+os.sep
'Path to this testing directory'

class Test_Module04b(unittest.TestCase):

    def setUp(self):
        # Loading and re-loading this each test case could get expensive.
        # Consider refactoring to a "static" module level variable that can be
        # treated as read-only as a balance between data loading overhead and
        # test independence; test cases that modify the HGSGrid object can
        # reload the data.
        self.g = HGSGrid(f'{TESTP}04b_Saturated_Fracture_Transport/module4b')

    def test_grid_size(self):
        gl = self.g.get_grid_lines()
        self.assertEqual(tuple(map(len,gl)),(101,2,51))

    def test_grid_indexing(self):
        # Note: In the fracture domain, the 3DNode# is off by one versus the
        # values below, as the HGS indices (Fortran) are 1-based and the HGSGrid
        # indices (Python) are 0-based
        self.assertEqual(self.g.find_node_index(0.,0.,20.), 8080,)
        self.assertEqual(self.g.find_node_index(0.,1.,20.), 8181,)

        self.assertEqual(tuple(self.g.find_grid_index(8080)), (0,0,40))
        self.assertEqual(tuple(self.g.find_grid_index(0.,0.,20.)), (0,0,40))

        self.assertEqual(tuple(self.g.find_grid_index(8181)), (0,1,40))
        self.assertEqual(tuple(self.g.find_grid_index(0.,1.,20.)), (0,1,40))

if __name__ == '__main__':
    unittest.main()
