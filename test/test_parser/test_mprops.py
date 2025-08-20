#!/usr/bin/env python
"""Test functions in `pyhgs.parser.mprops`"""


import unittest
import os

from hgstools.pyhgs.parser.mprops import parse as MProps

import logging
lstlog = logging.getLogger('pyhgs.parser.mprops')

TSTDIR=os.path.abspath(os.path.dirname(__file__))+os.path.sep
SIMDIR=os.path.abspath(os.path.join(TSTDIR,'..','test_sims'))+os.path.sep

class Test_MProps_Parser(unittest.TestCase):

    def test_file1(self):

        d = MProps(os.path.join(
                SIMDIR,'04b_very_coarse_mesh','module4b.mprops'))

        self.assertEqual(len(d), 2)

        with self.subTest('porosity values'):
            self.assertEqual(d['porous medium']['porosity'], 0.3)
            self.assertEqual(d['porous medium 2']['porosity'], 0.2)

        with self.subTest('case insensitive'):
            self.assertTrue('porous medium' in d)
            self.assertTrue('Porous Medium' in d)



if __name__=='__main__':
    unittest.main()
