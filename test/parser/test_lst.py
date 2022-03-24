#!/usr/bin/env python


import os
import unittest
from itertools import islice

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
        self.assertEqual(p.get_ec(), 0)

        p = LSTFileParser(DATDIR+'fail')
        self.assertEqual(p.get_ec(), 1)

    def test_fluid_balance(self):

        p = LSTFileParser(DATDIR+'good')

        fbdata = p.get_fluid_balance(1)

        self.assertEqual(len(fbdata), 5)

        self.assertTrue(abs(fbdata['Hinit_1'][0]-0.1171435366 < 1e-6))
        self.assertEqual(fbdata['Hinit_1'][3],'porous_media')
        self.assertTrue(abs(fbdata['TOTAL'][2]-(0.0000002241) < 1e-6))

    def test_iter_errors(self):

        p = LSTFileParser(DATDIR+'good')
        self.assertEqual(len(list(p.iter_errors())), 0)

        p = LSTFileParser(DATDIR+'fail')
        self.assertEqual(len(list(p.iter_errors())), 0)

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


if __name__ == '__main__':
    unittest.main()
