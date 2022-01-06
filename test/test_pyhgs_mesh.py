#!/usr/bin/env python
"""Test miscellaneous parts of `pyhgs.mesh`"""

import unittest

from collections import defaultdict
from decimal import Decimal
import numpy as np

from pyhgs.mesh import make_supersample_distance_groups, HGSGrid

PLACES = Decimal('0.001')
def D_CO(v, pl=PLACES):
    """Convert value(s) to Decimal objects"""
    if type(v) in [ str, float, int ]:
        return Decimal(v).quantize(pl)
    return list(Decimal(dd).quantize(pl) for dd in v)

class TestPYHGSMesh(unittest.TestCase):


    def test_supersample_groupings_grid(self):

        dx_floats = [
            0.5, #0
            0.5,
            0.4,
            0.05,
            0.05,
            0.05, #5
            0.05,
            0.4,
            0.5,
            0.5,
            0.5, #10
            0.4,
            0.05,
            0.05] # n=14

        dx = D_CO(dx_floats)

        self.assertEqual(
            make_supersample_distance_groups(dx,D_CO(0.05)),
            [[0,0], [1,1], [2,2], [3,4], [4,5], [5,6], [6,7], [7,7], [8,8],
            [9,9], [10,10], [11,11], [12,13], [13,14],],
        )

        self.assertEqual(
            make_supersample_distance_groups(dx,D_CO(0.5)), [[0,1], [1,2],
            [2,5], [3,7], [5,8], [8,9], [9,10], [10,11], [11,14],],
        )

        self.assertEqual(
            make_supersample_distance_groups(dx,D_CO(1.0)),
            [[0,2], [1,5], [2,8], [5,9], [8,10], [9,11], [10,14],],
        )


    def test_yield_fx_in_ssgrp(self):

        g = object.__new__(HGSGrid)   # uninitialized object

        g.elshape = (2,1,1,)
        ssranges = (
                [(0,1), (1,2), ],
                [(0,1,), ],
                [(0,1,), ],
            )
        adj = defaultdict(list, [(0,[0,]), (1,[0,]),] )

        act = g._yield_fx_in_ssgrp(ssranges, adj)
        des = 2*[[0,],]
        self.assertEqual([*act], des)



        g.elshape = (4,1,1)
        ssranges = (
                [(0,2), (2,3),],
                [(0,1,),],
                [(0,1,),],
            )

        act = g._yield_fx_in_ssgrp(ssranges, adj)
        des = [[0,], [],]
        self.assertEqual([*act],des)


        g.elshape = (4,2,2)
        ssranges = (
                [(0,2), (2,4),],
                [(0,2),],
                [(0,2),],
            )
        adj = defaultdict(list, [
            (0,[0,]),
            (1,[0,6,]),
            (2,[2,]),
            (3,[3,]),
            (4,[1,]),
            (5,[1,7,]),
            (6,[2,]),
            (7,[3,]),
            (9,[6,]),
            (10,[4,]),
            (11,[5,]),
            (13,[7,]),
            (14,[4,]),
            (15,[5,]),
        ] )

        des = [ [0,1,6,7,], [2,3,4,5,], ]

        act = g._yield_fx_in_ssgrp(ssranges, adj)
        self.assertEqual([*act],des)

if __name__ == '__main__':
    unittest.main()
