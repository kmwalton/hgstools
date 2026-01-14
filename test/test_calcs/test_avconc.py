#!/usr/bin/env python

import os
import unittest
import numpy as np

import warnings

from hgstools.pyhgs.test import sims_join
from hgstools.pyhgs.calcs import AvCalc
from hgstools.pyhgs.mesh import HGSGrid, Domain
from hgstools.pyhgs.parser.eco import EcoFile

import logging
calc_logger = logging.getLogger('pyhgs.calcs')
calc_logger.addHandler(logging.StreamHandler())

class Test_Av(unittest.TestCase):

    @unittest.skipIf( not os.path.isfile(sims_join('04b_very_coarse_mesh','module4bo.coordinates_pm')),
            'Could not find source simulation')
    def test_small_selection(self):
        """

        First, do a manual calculation of the average concentration for this
        "small selection". Then, use avcalc algorithm and compare results.
        """

        # load grid
        _d = sims_join('04b_very_coarse_mesh','module4b')
        g = HGSGrid(_d)
        blockspec = '0 5 0 .5 19 21'


        # determine static data
        # node & element indices, element volumes, porosity
        el_nd_ind = []
        el_V = []
        el_phi = []
        el_w = []
        for d in g.domains():
            # determine elements
            (el, nd)= g.choose_elements_block(blockspec, True, True, d)
            el_nd_ind.append( (el, nd,) )

            # record zones for selected elements
            _zn = g.get_elements_data(d)['porosity']
            el_phi.append(_zn.ravel('F')[el])

            # record volumes
            _V = g.get_element_volumes(d)
            el_V.append(_V.ravel('F')[el])

            with warnings.catch_warnings():
                #targeting the warning of the 0-based fracture node indices
                warnings.simplefilter('ignore')
                _inc = g.get_elements_data(d)['inc']

            # get the weight of each element (the count of times it intersects
            # blocspec)
            el_w.append(np.zeros(len(el), dtype=int))
            _ndset = set(nd)
            for i,iel in enumerate(el):
                _incset = set(_inc[iel])
                el_w[-1][i] = len(_ndset & _incset)
                    

        # handy aliases from static data
        el_pm = el_nd_ind[0][0]
        el_fx = el_nd_ind[1][0]
        nd_pm = el_nd_ind[0][1]
        nd_fx = el_nd_ind[1][1]
        nd_fx_pm = [g.hgs_fx_nodes['link_frac2pm'][i] for i in nd_fx]


        # determine possibly transient data
        el_q = []
        _el_q_pm = None
        _el_q_fx = None

        # pm
        _q = g.get_element_vals(g.prefix+'o.q_pm.0001','pm')
        _qq = np.zeros((len(el_pm), 3))
        for i in range(3):
            _qq[:,i] = _q[:,:,:,i].ravel('F')[el_pm]

        _el_q_pm = np.linalg.norm(_q,axis=3)
        el_q.append(np.linalg.norm(_qq,axis=1))

        # fx
        _q = g.get_element_vals(g.prefix+'o.v_frac.0001','frac')
        _qq = np.zeros((len(el_fx), 3))
        for i in range(3):
            _qq[:,i] = _q[:,i].ravel('F')[el_fx]

        _el_q_fx = np.linalg.norm(_q,axis=1)
        el_q.append(np.linalg.norm(_qq,axis=1))

        

        # intermediate checks
        self.assertEqual(len(el_pm), 2)
        self.assertEqual(len(el_fx), 1)

        self.assertEqual(el_V[0][0], 80.)
        self.assertEqual(el_V[0][1], 50.)
        self.assertEqual(el_V[1][0], 1e-3)

        self.assertEqual(el_phi[0][0], 0.3)
        self.assertEqual(el_phi[1][0], 1.0)

        self.assertTrue(
            abs(el_q[0][0] - 1.8819642007114837e-10) < 1e-17)
        self.assertTrue(
            abs(el_q[1][0] - 0.00115586) < 1e-9)

        self.assertTrue(np.allclose(el_w[0], [1,1,]))
        self.assertTrue(np.allclose(el_w[1], [1,]))

        # get conc data:
        # concentration (nodal)
        c = g.get_nodal_vals(_d+'o.conc_pm.salt.0010')

        # put in array
        el_C = []
        el_C.append(g.get_element_vals(c,'pm').ravel('F')[el_pm])
        el_C.append(g.get_element_vals(c,'frac')[el_fx])


        # flatten arrays
        _w = np.hstack(el_w)
        _V = np.hstack(el_V)
        _phi = np.hstack(el_phi)
        _q = np.hstack(el_q)
        _C = np.hstack(el_C)

        # compute averages
        Cbar_naive = np.average(_C)
        Cbar_nodal_all_naive = np.average(
                np.hstack((c.ravel('F')[nd_pm], c.ravel('F')[nd_fx_pm])))
        Cbar_nodal_fx_naive = np.average(c.ravel('F')[nd_fx_pm])
        Cbar_Vweight = np.average(_C, weights=_w*_V*_phi)
        Cbar_qweight = np.average(_C, weights=_w*_q/_phi)


        calc = AvCalc(_d)

        c_pm = c
        c_fx = g.get_nodal_vals(c,'frac')

        with self.subTest('elemental'):
          self.assertAlmostEqual(
            Cbar_naive,
            calc.average(blockspec, [c_pm,c_fx,], 'arith_el'))

        with self.subTest('nodal'):
          self.assertAlmostEqual(
            Cbar_nodal_all_naive,
            calc.average(blockspec, [c_pm,c_fx,], 'arith_nd'))

        with self.subTest('nodal, fx only'):
          self.assertAlmostEqual(
            Cbar_nodal_fx_naive,
            calc.average(blockspec, [c_fx,], 'arith_nd', 'frac'))

        with self.subTest('volume weight'):
          self.assertAlmostEqual(
            Cbar_Vweight,
            calc.average(blockspec, [c_pm,c_fx,], 'porvol_el'),
            delta=1e-4)

        with self.subTest('flux weight'):
          #_levsav = calc_logger.level
          #calc_logger.setLevel(logging.DEBUG)
          self.assertAlmostEqual(
            Cbar_qweight,
            calc.average(blockspec, [c_pm,c_fx,], [_el_q_pm.flatten('F'),
                _el_q_fx,]),
            delta=1e-4)
          #calc_logger.setLevel(_levsav)

    @unittest.skipIf( not os.path.isfile(sims_join('04b_very_coarse_mesh','module4bo.eco')),
            'Could not find source simulation')
    def test_doms_list(self):
        # load grid
        _d = sims_join('04b_very_coarse_mesh','module4b')
        calc = AvCalc(_d)
        self.assertEqual(calc._doms_list('pm'), [Domain.PM,])
        self.assertEqual(calc._doms_list(Domain.PM), [Domain.PM,])
        self.assertEqual(calc._doms_list('frac'), [Domain.FRAC,])
        self.assertEqual(calc._doms_list(Domain.FRAC), [Domain.FRAC,])
        self.assertEqual(calc._doms_list('all'), [Domain.PM, Domain.FRAC,])
        self.assertEqual(calc._doms_list(['pm','frac']), [Domain.PM, Domain.FRAC,])
        self.assertEqual(calc._doms_list(['frac','pm']), [Domain.PM, Domain.FRAC,])

    @unittest.skipIf( not os.path.isfile(sims_join('04b_very_coarse_mesh','module4bo.eco')),
            'Could not find source simulation')
    def test_get_block(self):
        # load grid
        _d = sims_join('04b_very_coarse_mesh','module4b')
        calc = AvCalc(_d)

        got_bl = calc._get_block('0 5 0 .5 19 21')
        self.assertEqual(len(got_bl), 4)
        self.assertEqual(len(got_bl[Domain.PM][0]), 2)
        self.assertEqual(len(got_bl[Domain.PM][1]), 1)
        self.assertEqual(got_bl['ne'], 3)
        self.assertEqual(got_bl['nn'], 2)

        got_bl = calc._get_block('0 10 0 1 12 20')
        self.assertEqual(len(got_bl), 4)
        self.assertEqual(len(got_bl[Domain.PM][0]), 6)
        self.assertEqual(len(got_bl[Domain.PM][1]), 8)
        self.assertEqual(len(got_bl[Domain.FRAC][0]), 2)
        self.assertEqual(len(got_bl[Domain.FRAC][1]), 4)

if __name__=='__main__':
    unittest.main()

