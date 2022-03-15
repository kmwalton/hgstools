#!/usr/bin/env python


import unittest

from pyhgs.parser.lst import LSTFileParser


import logging
lstlog = logging.getLogger('pyhgs.parser.lst')

class Test_LST_Parser(unittest.TestCase):


    def test_name_equiv(self):

        p1 = LSTFileParser('steady_flow_transient_transport')
        p0 = LSTFileParser('steady_flow_transient_transporto.lst')

        self.assertEqual(p0._fnin, p1._fnin)

    def test_n_ts_read(self):

        p = LSTFileParser('steady_flow_transient_transporto.lst')
        self.assertEqual(p.get_n_ts(),10)

        p = LSTFileParser('sflow_ttransport_w0solv_erroro.lst')
        self.assertEqual(p.get_n_ts(),30)

    def test_ss_flow(self):

        p = LSTFileParser('steady_flow_transient_transporto.lst')
        self.assertEqual(p.ss_flow(), 0)

    def test_ec(self):

        p = LSTFileParser('steady_flow_transient_transporto.lst')
        self.assertEqual(p.get_ec(), 0)

if __name__ == '__main__':
    unittest.main()
