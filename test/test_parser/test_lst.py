#!/usr/bin/env python
"""Test functions in `pyhgs.parser.lst`"""

import unittest

import os
from itertools import islice
from bisect import bisect_left
from operator import itemgetter

from pyhgs.parser.lst import LSTFileParser

import logging
lstlog = logging.getLogger('pyhgs.parser.lst')

DATDIR=os.path.abspath(os.path.dirname(__file__))+os.path.sep

class Test_LST_Parser(unittest.TestCase):


    def test_name_equiv(self):

        p1 = LSTFileParser(DATDIR+'steady_flow_transient_transport')
        p0 = LSTFileParser(DATDIR+'steady_flow_transient_transporto.lst')

        self.assertEqual(p0._fnin, p1._fnin)

    def test_n_ts_read(self):

        p = LSTFileParser(DATDIR+'steady_flow_transient_transporto.lst')
        self.assertEqual(p.get_n_ts(),10)

        p = LSTFileParser(DATDIR+'sflow_ttransport_w0solv_erroro.lst')
        self.assertEqual(p.get_n_ts(),24)

        p = LSTFileParser(DATDIR+'fail')
        self.assertEqual(p.get_n_ts(), 0)

    def test_ss_flow(self):

        p = LSTFileParser(DATDIR+'steady_flow_transient_transporto.lst')
        self.assertEqual(p.ss_flow(), True)

        p = LSTFileParser(DATDIR+'fail')
        self.assertEqual(p.ss_flow(), True)

    def test_ec(self):

        p = LSTFileParser(DATDIR+'steady_flow_transient_transporto.lst')
        with self.subTest('good exit'):
            self.assertEqual(p.get_ec(), 0)
            self.assertTrue(p._has_sim_report)

        p = LSTFileParser(DATDIR+'fail')
        with self.subTest('bad exit'):
            self.assertEqual(p.get_ec(), 1)
            self.assertFalse(p._has_sim_report)

    def test_fluid_balance(self):

        p = LSTFileParser(DATDIR+'good')

        fbdata = p.get_fluid_balance(1)

        self.assertEqual(len(fbdata), 5)

        self.assertTrue(abs(fbdata['Hinit_1'][0]-0.1171435366 < 1e-6))
        self.assertEqual(fbdata['Hinit_1'][3],'porous_media')
        self.assertTrue(abs(fbdata['TOTAL'][2]-(0.0000002241) < 1e-6))

    def test_fluid_balance_trans(self):
        """Test parsing of fluid balance in transient flow and transport"""
        p = LSTFileParser(DATDIR+'tflow_ttransport_pm_frac_well')
        fbdata = p.get_fluid_balance(12)
        self.assertTrue(abs(fbdata['Hinit_6'][1]+1.1862877408) < 1e-6)

    def test_mass_balance_trans(self):
        """Test parsing of mass balance in transient flow and transport

        Simulation has only one tracer.
        """
        p = LSTFileParser(DATDIR+'tflow_ttransport_pm_frac_well')
        mbdata = p.get_mass_balance(13)

        self.assertAlmostEqual(
                mbdata['Fixed concentration nodes'][0], 6.4463364331)

    def test_iter_errors(self):

        with self.subTest('no errors, goodo.lst'):
            p = LSTFileParser(DATDIR+'good')
            self.assertEqual(len(list(p.iter_errors())), 0)

        with self.subTest('no errors, failo.lst'):
            p = LSTFileParser(DATDIR+'fail')
            self.assertEqual(len(list(p.iter_errors())), 0)

        with self.subTest('errors, sflow_...o.lst'):
            p = LSTFileParser(DATDIR+'sflow_ttransport_w0solv_erroro.lst')
            self.assertEqual(len(list(p.iter_errors())), 1)
            
            self.assertEqual(len(list(p.iter_errors(1))), 1)
            self.assertEqual(len(list(p.iter_errors(2))), 0)

    def test_n_solver_iter(self):

        p = LSTFileParser(DATDIR+'good')
        self.assertEqual(p.get_transport_solver_iterations(1082), 4)

        p = LSTFileParser(DATDIR+'fail')
        self.assertEqual(p.get_n_ts(), 0)
        self.assertEqual(p.get_flow_solver_iterations(0), 20000)

    def test_ts_solution_times(self):
        p = LSTFileParser(DATDIR+'good')
        
        self.assertAlmostEqual(p.get_ts_time(0),0.0)
        self.assertAlmostEqual(p.get_ts_time(1),0.01)
        self.assertAlmostEqual(p.get_ts_time(4),0.15)

        ts_list = [ i for i in islice(p.get_ts_time(),6) ]

        for act,des in zip(ts_list, [0., 0.01, 0.03, 0.07, 0.15, 0.31,]):
            self.assertAlmostEqual(act,des)

        self.assertEqual(p.get_n_ts(),1082)
        ts_list = list( p.get_ts_time() )
        self.assertEqual(len(ts_list),1082+1)


        self.assertAlmostEqual(p.get_ts_dtime(1), 9.9999997e-3)
        self.assertAlmostEqual(p.get_ts_dtime(1072), 229.357142421641)

    def test_ts_solution_times_trans(self):
        p = LSTFileParser(DATDIR+'tflow_ttransport_pm_frac_well')

        self.assertAlmostEqual(p.get_ts_time(12), 10.00078125)
        self.assertAlmostEqual(p.get_ts_dtime(12), 0.00078125)

        with self.assertRaises(ValueError):
            self.assertAlmostEqual(p.get_ts_dtime(0), 10.00078125)


    def test_gtt(self):
        """Test global target times harvesting"""
        p = LSTFileParser(DATDIR+'good')
        
        # total number of times
        gtt_actual = p.get_global_target_times()
        gtt_actual_times = list( map(itemgetter(0), gtt_actual) )
        self.assertEqual(len(gtt_actual),88)

        # specific time
        tgt_desired = 37412.5
        tgt_idesired = 251
        
        i = bisect_left(gtt_actual_times, tgt_desired)
        self.assertAlmostEqual(gtt_actual[i][0],tgt_desired)
        self.assertEqual(gtt_actual[i][1],tgt_idesired)


    def test_mass_storage(self):

        p = LSTFileParser(DATDIR+'tflow_ttransport_pm_frac_wello.lst')
        ms = p.get_mass_storage(18)
        self.assertEqual(len(ms),8)
        self.assertAlmostEqual(p.get_ts_time(18), 10.0040609323282)
        self.assertEqual(len(ms['Porous medium']),4)
        self.assertAlmostEqual(ms['Porous medium'][0], 0.5389514566)
        self.assertAlmostEqual(ms['Porous medium'][1], 0.0)
        self.assertAlmostEqual(ms['Porous medium'][2], 0.0)
        self.assertAlmostEqual(ms['Porous medium'][3], 0.0)
        self.assertAlmostEqual(ms['Discrete fractures'][0], 0.0038146)
        self.assertAlmostEqual(ms['Discrete fractures'][2], 0.0)
        self.assertAlmostEqual(ms['Wells'][0], 0.0)
        self.assertEqual(len(ms['Wells']),1)

        # find the key instead of using whole string
        k = next( k for k in ms.keys() if k.startswith('NET2'))
        self.assertAlmostEqual( ms[k], 6.4313503119)


        
        p = LSTFileParser(DATDIR+'steady_flow_transient_transporto.lst')
        ms = p.get_mass_storage(10)
        self.assertEqual(len(ms),7)
        self.assertAlmostEqual(ms['Porous medium'][0], 0.0)
        self.assertAlmostEqual(ms['Porous medium'][2], 0.0)
        with self.assertRaises(KeyError):
            ms['Wells']

if __name__ == '__main__':
    unittest.main()
