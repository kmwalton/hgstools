#!/usr/bin/env python
"""Testing of ofrac2hgs module"""

import unittest

from decimal import Decimal
from hgstools.bin.ofrac2hgs import *

N_APERT_DIG = Decimal('0.000001')

def D(v,new_prec):
    """Return a decimal with the specified quantization/precision"""
    return Decimal(v).quantize(new_prec)

def D_AP(v):
    """Return a decimal with the required number of digits for apertures"""
    return D(v,N_APERT_DIG)

class TestOfrac2hgs(unittest.TestCase):
    def test_apQuanize(self):

        v = 0.123456
        self.assertEqual(apQuantize(v, 6), D_AP(0.123456))

        self.assertEqual(apQuantize(v, 5), D_AP(0.123460))

        breakpoint()
        self.assertEqual(apQuantize(0.000182, 5), D_AP(0.000180))
        self.assertEqual(apQuantize(0.000186, 5), D_AP(0.000190))
        self.assertEqual(apQuantize(0.000188, 5), D_AP(0.000190))

        self.assertEqual(apQuantize(0.000002, 4), D_AP(0.000050))
        self.assertEqual(apQuantize(0.000099, 4), D_AP(0.000050))


if __name__ == '__main__':
    unittest.main()
