#!/usr/bin/env python
"""Test functions in `pyhgs.parser.hgseco`"""

import unittest

import os
from pyhgs.parser.hgseco import HGSEcoFileParser

DATDIR=os.path.abspath(os.path.dirname(__file__)+'/..')+os.path.sep

class Test_ECO_Parser(unittest.TestCase):

    def test_ext_detect(self):
        p = HGSEcoFileParser(DATDIR+'goodo.eco')
        self.assertEqual(p.get_n_zones(), 1)
        p = HGSEcoFileParser(DATDIR+'good')
        self.assertEqual(p.get_n_zones(), 1)

    def test_n_zone_reads(self):
        p = HGSEcoFileParser(DATDIR+'goodo.eco')
        self.assertEqual(p.get_n_zones(), 1)

        p = HGSEcoFileParser(DATDIR+os.path.join(
            'test_sims', 'undulating_coarse', 'undf'))

        self.assertEqual(p.get_n_zones(), 2)


    def test_zone_prop_reads(self):
        with self.subTest('goodo.eco'):
            p = HGSEcoFileParser(DATDIR+'goodo.eco')
            self.assertEqual(p.get_pm_zone_properties(), [(1,'sandstone',),])

        with self.subTest('undfo.eco'):
          p = HGSEcoFileParser(DATDIR+os.path.join(
            'test_sims', 'undulating_coarse', 'undf',))

          self.assertEqual(
            p.get_pm_zone_properties(),
            [(1,'porous medium',),
             (2,'porous medium 2',),])

    def test_get_output_times(self):
        p = HGSEcoFileParser(DATDIR+'goodo.eco')
        self.assertEqual(len(p.get_output_times()), 88)


if __name__ == '__main__':
    unittest.main()
