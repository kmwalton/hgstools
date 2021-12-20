#!/usr/bin/env python
"""Python unittests involving 04b_Saturated_Fracture_Transport"""

import unittest
import os
from itertools import product
import numpy as np
import numpy.testing as nptest

import pyhgs
from pyhgs.parser import parse
from pyhgs._test import skip_if_no_sim_output 
from pyhgs.mesh import (HGSGrid, Domain)

TESTP = os.path.dirname(__file__)
'Path to this testing directory'

SIM_PREFIX = os.path.join(TESTP,'04b_very_coarse_mesh','module4b')
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
        self.assertEqual(self.g.shape,(6,2,5))
        self.assertEqual(self.g.elshape,(5,1,4))

        self.assertEqual(self.g.nn, 60)
        self.assertEqual(self.g.ne, 20)
        self.assertEqual(self.g.nfn, 24)
        self.assertEqual(self.g.nfe, 11)

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

            ipmgrid = self.g.eli2elg(ipm)
            with self.subTest(pm_element=ipmgrid):
                self.assertEqual(desired, actual[ipmgrid])

        # spot-check some nodes with no neighbours
        for ipm in [0, 1, 5, 6, 14, 19,]:
            with self.subTest(pm_element=ipm):
                self.assertEqual([], actual[ipm])

        # try invalid < 0 index
        with self.assertRaises(ValueError):
            actual[-1]

        # try invalid > index
        with self.assertRaises(ValueError):
            actual[20]

        self.assertTrue(self.pmn2fxn[0] < 0)


    def test_read_pm_conc(self):

        desired_t = 31536000.
        some_desired = {
            0:7.07029379487522647e-10,
            3:9.86041470696363831e-08,
            9:9.86041470696363831e-08,
            14:1.96902050220160163e-07,
            19:2.25981438006783719e-07,
            34:2.82556868569372455e-07,
            46:1.07707273855339736e-05,
            59:3.59244994863061606e-11,
        }
        some_more_desired = {
           (3,0,3):0.00341675779782235622,
           (3,0,2):0.000466027995571494102,
           (2,0,2):3.52357856172602624e-05,
           (2,0,3):0.0211286600679159164,
           (2,1,2):3.43717256328091025e-05,
           (2,1,3):0.0209547057747840881,
           (3,1,3):0.00340801011770963669,
           (3,1,2):0.000465662626083940268,
        }

        FN = f'{SIM_PREFIX}o.conc_pm.salt.0010'
        d = parse(FN)

        self.assertTrue(abs(float(d['ts'])-desired_t)<1e-6,
                'Correct file/time')

        with self.subTest(read='PM node concentration'):
            cnpm = self.g.get_nodal_vals(FN)

            act = cnpm.flatten(order='F')[list(some_desired.keys())]
            des = list(some_desired.values())

            nptest.assert_allclose(act, des, rtol=0.01,)

            for ind,des in some_more_desired.items():
                nptest.assert_allclose(cnpm[ind], des, rtol=0.01,)


        with self.subTest(read='PM element concentration'):
            cepm = self.g.get_element_vals(FN)

            act = cepm[2,0,2]
            des = sum(some_more_desired.values())/8.
            nptest.assert_allclose(act, des, rtol=0.01)


    def _assertEqual_fx_nodal_0010(self, actual):
            some_desired = {
                3:9.86041470696363831e-08,
                9:9.85499681860346755e-08,
                14:1.96902050220160163e-07,
                21:8.28335832920856774e-05,
            }
            
            act = actual[list( self.pmn2fxn[k] for k in some_desired.keys())]
            des = list(some_desired.values())
            nptest.assert_allclose( act, des, rtol=0.01,)

    def _assertEqual_fx_elem_0010(self, actual):
        some_desired = {
            0:0.546559039503335953,
            5:0.000274344090939848684,
        }

        act = actual[list(k for k in some_desired.keys())]
        des = list(some_desired.values())
        nptest.assert_allclose( act, des, rtol=0.01,)

    def test_read_frac_conc(self):
        """Read an example PM and fracture-domain concentration file"""
        with self.subTest(read='Frac node concentration'):
            cnfx = self.g.get_nodal_vals(f'{SIM_PREFIX}o.conc_pm.salt.0010',
                Domain.FRAC)
            self.assertEqual(cnfx.size,24)
            self._assertEqual_fx_nodal_0010(cnfx)

            with self.assertRaises(IndexError):
                cnfx[ self.pmn2fxn[0] ] # test where there is no fracture

        with self.subTest(read='Frac element concentration'):
            cefx = self.g.get_element_vals(f'{SIM_PREFIX}o.conc_pm.salt.0010',
                Domain.FRAC)
            self._assertEqual_fx_elem_0010(cefx)

    def test_read_frac_conc_from_pm(self):
        """Test fracture concentrations from PM nodal concentration inputs"""
        cnpm = self.g.get_nodal_vals(f'{SIM_PREFIX}o.conc_pm.salt.0010')
        self.assertEqual(cnpm.shape,(6,2,5))
        self.assertEqual(cnpm.size,60)

        with self.subTest(read='nodal concentration'):
            cnfx = self.g.get_nodal_vals(cnpm, Domain.FRAC)
            self.assertEqual(cnfx.size,24)

        with self.subTest(read='elemental concentration'):
            cefx = self.g.get_element_vals(cnpm, Domain.FRAC)
            self.assertEqual(cefx.size,11)
            self._assertEqual_fx_elem_0010(cefx)

    def test_read_flux(self):
        """Read example PM and fracture-domain velocity and flux files"""

        pmv = self.g.get_element_vals(f'{SIM_PREFIX}o.v_pm.0001')
        pmq = self.g.get_element_vals(f'{SIM_PREFIX}o.q_pm.0001')
        fxv = self.g.get_element_vals(f'{SIM_PREFIX}o.v_frac.0001',
                Domain.FRAC)


    def test_calc_flux_mag(self):
        """Read example PM and fracture-domain velocity and flux files"""

        pmq = self.g.get_element_vals(f'{SIM_PREFIX}o.q_pm.0001')
        fxv = self.g.get_element_vals(f'{SIM_PREFIX}o.v_frac.0001',
                Domain.FRAC)

        pmqmag = np.sqrt(np.sum(pmq**2, axis=3))
        fxvmag = np.sqrt(np.sum(fxv**2, axis=1))

        pmqmag_des = np.array([
          [[2.0518945e-10, 1.8794906e-10, 1.8819642e-10, 1.5892040e-10],],
          [[2.6000405e-10, 2.6866750e-10, 2.2190369e-10, 1.5802624e-10],],
          [[5.2409233e-11, 1.8065756e-10, 2.1618764e-10, 1.5657167e-10],],
          [[1.2459916e-10, 1.7248690e-10, 1.8916432e-10, 3.8519826e-11],],
          [[1.2648921e-10, 1.6842430e-10, 2.7908487e-10, 2.9564320e-10],]])
        fxvmag_des = np.array([
          1.1558607e-03, 1.1493608e-03, 1.1380099e-03, 5.1349805e-05,
          1.3408824e-06, 1.0958962e-03, 1.0914268e-03, 2.5401107e-06,
          3.7644098e-05, 1.1310030e-03, 1.1501204e-03])

        nptest.assert_allclose(pmqmag, pmqmag_des, rtol=1e-5)
        nptest.assert_allclose(fxvmag, fxvmag_des, rtol=1e-5)


    def test_supersample_distance_groups(self):

        des_pm10xgrps = [[0,1], [1,2], [2,4], [4,4],]
        des_pm10ygrps = [[0,1],]
        des_pm10zgrps = [[0,1], [1,2], [2,3], [3,4],]

        # test individual dimensions
        # PM only
        ssdg = pyhgs.mesh.make_supersample_distance_groups
        self.assertEqual(ssdg([10.,10.,5.,5.,20.], 10.), des_pm10xgrps)
        self.assertEqual(ssdg([1.,], 10.), des_pm10ygrps)
        self.assertEqual(ssdg([6.,6.,8.,5.], 10.), des_pm10zgrps)

        # test sequence of PM ss groups
        des_pm10 = list(product(des_pm10xgrps, des_pm10ygrps, des_pm10zgrps))
        act_pm10 = list(self.g.iter_supersample_distance_groups(10.))
        self.assertEqual(act_pm10, des_pm10)

        # test sequence
        des_fx10 = [
            [], [], [0,], [0,], # z-index is fastest!
            [], [], [1,], [1,],
            [4,8,9,], [5,8,9,], [2,3,6,], [2,3,7,],
            [], [], [], [], # note 0-length PM groups give NO fracture 10
        ]
        act_fx10 = list( t[0] for t in 
                self.g.iter_supersample_distance_groups(10., [Domain.FRAC,]))

        self.assertEqual(act_fx10, des_fx10)


if __name__ == '__main__':
    unittest.main()
