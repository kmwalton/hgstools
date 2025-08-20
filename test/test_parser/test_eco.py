#!/usr/bin/env python
"""Test functions in `pyhgs.parser.hgseco`"""

import unittest

import os
from hgstools.pyhgs.parser.eco import EcoFile
from hgstools.pyhgs.test import sims_join

DATDIR=os.path.abspath(os.path.dirname(__file__)+'/..')+os.path.sep

class Test_ECO_Parser(unittest.TestCase):

    @unittest.skipIf( not os.path.isfile(sims_join('goodo.eco')),
            'Could not find testing source file')
    def test_ext_detect(self):
        p = EcoFile(DATDIR+'goodo.eco')
        self.assertEqual(p.get_n_zones(), 1)
        p = EcoFile(DATDIR+'good')
        self.assertEqual(p.get_n_zones(), 1)

    @unittest.skipIf( not os.path.isfile(sims_join('undulating_coarse','undfo.eco')),
            'Could not find testing source file')
    def test_n_zone_reads(self):
        p = EcoFile(DATDIR+'goodo.eco')
        self.assertEqual(p.get_n_zones(), 1)

        p = EcoFile(DATDIR+os.path.join(
            'test_sims', 'undulating_coarse', 'undf'))

        self.assertEqual(p.get_n_zones(), 2)


    @unittest.skipIf( not os.path.isfile(sims_join('undulating_coarse','undfo.eco')),
            'Could not find testing source file')
    def test_zone_prop_reads(self):
        with self.subTest('goodo.eco'):
            p = EcoFile(DATDIR+'goodo.eco')
            self.assertEqual(p.get_pm_zone_properties(), [(1,'sandstone',),])

        with self.subTest('undfo.eco'):
          p = EcoFile(DATDIR+os.path.join(
            'test_sims', 'undulating_coarse', 'undf',))

          self.assertEqual(
            p.get_pm_zone_properties(),
            [(1,'porous medium',),
             (2,'porous medium 2',),])

    @unittest.skipIf( not os.path.isfile(DATDIR+'goodo.eco'),
            'Could not find testing source file')
    def test_get_output_times(self):
        p = EcoFile(DATDIR+'goodo.eco')
        self.assertEqual(len(p.get_output_times()), 88)


if __name__ == '__main__':
    unittest.main()
