#!/usr/bin/env python
"""Test miscellaneous parts of `pyhgs.mesh`"""

import unittest

import os
import warnings
from collections import defaultdict
from decimal import Decimal
import numpy as np

from pyhgs._test import skip_if_no_sim_output
from pyhgs.mesh import make_supersample_distance_groups, HGSGrid, Domain

TESTP = os.path.dirname(__file__)+os.sep
'Path to this testing file'

PLACES = Decimal('0.001')
def D_CO(v, pl=PLACES):
    """Convert value(s) to Decimal objects"""
    if type(v) in [ str, float, int ]:
        return Decimal(v).quantize(pl)
    return list(Decimal(dd).quantize(pl) for dd in v)

class TestPYHGSMesh(unittest.TestCase):

    def test_Domain(self):
        Domain.FRAC
        Domain.PM
        self.assertEqual(Domain.a2D('pm'), Domain.PM)
        self.assertEqual(Domain.a2D('PM'), Domain.PM)
        self.assertEqual(Domain.a2D('frac'), Domain.FRAC)
        self.assertEqual(Domain.a2D('FRAC'), Domain.FRAC)

        with self.assertRaises(ValueError):
            Domain.a2D('junk')

        with self.assertRaises(AttributeError):
            Domain.JUNK

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



class Test_HGSGrid(unittest.TestCase):

    @unittest.skipIf(
        skip_if_no_sim_output(
            TESTP+'test_sims/04b_very_coarse_mesh/module4b'),
        'HGS output missing')
    def test_choose_nodes_block_pm(self):

        g = HGSGrid(
                TESTP+'test_sims/04b_very_coarse_mesh/module4b')

        self.assertEqual(
            g.choose_nodes_block('0,0,0,0,0,0'),
            [0,]
            )

        self.assertEqual(
            g.choose_nodes_block('0,10,0,0,0,0'),
            [0,1]
            )

        self.assertEqual(
            g.choose_nodes_block('0,10,0,1,0,0'),
            [0,1,6,7]
            )

        self.assertEqual(
            g.choose_nodes_block('0,10,0,0,0,13'),
            [0,1,12,13,24,25]
            )

        self.assertEqual(
            g.choose_nodes_block('0,10,0,1,0,6'),
            [0,1,6,7,12,13,18,19]
            )

        self.assertEqual(
            g.choose_nodes_block('1,9,0,1,0,6'),
            []
            )

    @unittest.skipIf(
        skip_if_no_sim_output(
            TESTP+'test_sims/04b_very_coarse_mesh/module4b'),
        'HGS output missing')
    def test_choose_nodes_block_fx(self):

        g = HGSGrid(
                TESTP+'test_sims/04b_very_coarse_mesh/module4b')


        fx_nodes = g.choose_nodes_block('0,10,0,1,19.9,20.1', dom='FRAC')
        pm_nodes = g.hgs_fx_nodes['link_frac2pm'][fx_nodes]
        self.assertEqual(list(pm_nodes), [36,37,42,43])


    @unittest.skipIf(
        skip_if_no_sim_output(
            TESTP+'test_sims/04b_very_coarse_mesh/module4b'),
        'HGS output missing')
    def test_node2el_pm(self):

        g = HGSGrid(TESTP+'test_sims/04b_very_coarse_mesh/module4b')

        with self.subTest('PM basic check'):
            _r =g._node2el()

            self.assertEqual(len(_r), 60, 'incidence for all nodes')

            self.assertEqual(_r[0], {0,}) 
            self.assertEqual(_r[13], {0,1,5,6}) 

        with self.subTest('PM basic check with set'):
            self.assertEqual(g._node2el(0)[0], {0,}) 
            self.assertEqual(g._node2el([0,])[0], {0,})

            _r = g._node2el([0,6,])
            self.assertEqual(len(_r), 2, 'incidence for only two nodes')
            self.assertEqual(_r[0], {0,})
            self.assertEqual(_r[6], {0,})
            self.assertEqual(_r[6], {0,})
            with self.assertRaises(KeyError):
                _r[2]

            self.assertEqual(g._node2el(59)[59], {19,}) 

        with self.subTest('Public method'):
            n2e = g.ni2eli('PM')
            self.assertEqual(n2e[53], {19,})
            self.assertEqual(n2e[59], {19,})
            self.assertEqual(n2e[13], {0,1,6,5})
            self.assertEqual(n2e[17], {4,9})

    @unittest.skipIf(
        skip_if_no_sim_output(
            TESTP+'test_sims/04b_very_coarse_mesh/module4b'),
        'HGS output missing')
    def test_node2el_fx(self):

        g = HGSGrid(TESTP+'test_sims/04b_very_coarse_mesh/module4b')

        with warnings.catch_warnings():
            #targeting the warning of the 0-based fracture node indices
            warnings.simplefilter('ignore')

            with self.subTest('Frac basic check'):
                _r =g._node2el(dom='FRAC')
                self.assertEqual(_r[0], {4,})
                self.assertEqual(_r[3], {4,5,8,9})

            with self.subTest('Frac basic check with set'):
                _r =g._node2el([0,], dom='FRAC')
                self.assertEqual(_r[0], {4,})

            with self.subTest('Frac basic check with set'):
                _pmi = [3,0,] #list(map(pm2fx, [15,0,]))
                _r =g._node2el(_pmi, dom='FRAC')

            with self.subTest('Public method'):
                n2e = g.ni2eli('FRAC')
                self.assertEqual(n2e[0], {4,})
                self.assertEqual(n2e[12], {0,})
                self.assertEqual(n2e[13], {0,1,})
                self.assertEqual(n2e[3], {4,5,8,9})
                self.assertEqual(n2e[7], {4,5,8,9})


    @unittest.skipIf(
        skip_if_no_sim_output(
            TESTP+'test_sims/04b_very_coarse_mesh/module4b'),
        'HGS output missing')
    def test_choose_elements_block_pm(self):

        g = HGSGrid(TESTP+'test_sims/04b_very_coarse_mesh/module4b')

        with self.subTest('full capture'):
            got = g.choose_elements_block('0,0,0,0,0,0')
            self.assertEqual(got, [], 'capture nothing')

            got = g.choose_elements_block('9,21,0,0,11,21')
            self.assertEqual(got, [], 'capture nothing')

            got = g.choose_elements_block('9,21,0,1,11,21')
            self.assertEqual(got, [11,], 'strict inside capture')

            got = g.choose_elements_block('10,20,0,1,12,20', False, True)
            self.assertEqual(
                    got,
                    ([11,],
                     [25,26,31,32,37,38,43,44,]))

        with self.subTest('partial capture'):
            got = g.choose_elements_block('10,20,0,1,12,20', True)
            self.assertEqual(got, [5,6,7,10,11,12,15,16,17,]) 

            got = g.choose_elements_block('10,20,0,1,12,20', True)
            self.assertEqual(got, [5,6,7,10,11,12,15,16,17,])

            got = g.choose_elements_block('10,20,0,1,12,20', True, True)
            self.assertEqual(
                    got,
                    ([5,6,7,10,11,12,15,16,17,],
                     [25,26,31,32,37,38,43,44,]))

    @unittest.skipIf(
        skip_if_no_sim_output(
            TESTP+'test_sims/04b_very_coarse_mesh/module4b'),
        'HGS output missing')
    def test_choose_elements_block_fx(self):

        g = HGSGrid(TESTP+'test_sims/04b_very_coarse_mesh/module4b')

        def pm2fx(i):
            """Warning-free wrapper of link_pm2fx"""
            with warnings.catch_warnings():
                #targeting the warning of the 0-based fracture node indices
                warnings.simplefilter('ignore')
                return g.hgs_fx_nodes['link_pm2fx'][i]
    
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

        with self.subTest('full capture'):
            got = g.choose_elements_block('25,25,0,1,0,6', dom='FRAC')
            self.assertEqual(got, [4,])

            got = g.choose_elements_block('10,20,0,1,20,20',False,True,'FRAC')
            self.assertEqual(got, ([1,], [13,14,18,19,]))

        with self.subTest('partial capture'):
            got = g.choose_elements_block('25,25,0,1,0,6', True, dom='FRAC')
            self.assertEqual(got, [4,5,8,9,])

            got = g.choose_elements_block('10,20,0,1,20,20',True,True,'FRAC')
            self.assertEqual(got, ([0,1,2,], [13,14,18,19,]))

        with self.assertRaises(ValueError):
            g.choose_elements_block('10,20,0,1,20,20', True, 'FRAC')

if __name__ == '__main__':
    unittest.main()
