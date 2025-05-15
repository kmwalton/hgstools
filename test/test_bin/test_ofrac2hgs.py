#!/usr/bin/env python
"""Testing of ofrac2hgs module"""

import os
import unittest
import tempfile

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

        self.assertEqual(apQuantize(0.000182, 5), D_AP(0.000180))
        self.assertEqual(apQuantize(0.000186, 5), D_AP(0.000190))
        self.assertEqual(apQuantize(0.000188, 5), D_AP(0.000190))

        self.assertEqual(apQuantize(0.000002, 4), D_AP(0.000050))
        self.assertEqual(apQuantize(0.000099, 4), D_AP(0.000050))

    def test_truncate_domain(self):

        from ofrac.ofracs import OFracGrid

        g = OFracGrid(domainSize=(1.,1.,1.),
              fx=[
                (0., 1., 0., 1., 0.5, 0.5, 0.001),
                (0., 1., 0.5, 0.5, 0., 1., 0.002),
                (0.5, 0.5, 0., 1., 0., 1., 0.003),
            ],
        )
        self.assertEqual(g.getFxCount(),3)

        # get a unique name
        fp = tempfile.NamedTemporaryFile(prefix='dfn_', suffix='.pkl', dir='.', delete=False)
        fp.close()
        OFracGrid.pickleTo(g, fp.name)

        # test things!
        h = RFG(fp.name, domainSize=(1,1,0.5))
        self.assertEqual(h.fxnet.getFxCount(),3)

        h = RFG(fp.name, domainSize=(1,1,0.51))
        self.assertEqual(h.fxnet.getFxCount(),3)

        h = RFG(fp.name, domainSize=(1,1,0.49))
        self.assertEqual(h.fxnet.getFxCount(),2)

        h = RFG(fp.name, translate=(-0.5, 0, 0), domainSize=(0.5,0.5,0.5),
                nudgeTo=0.1)
        self.assertEqual(h.fxnet.getFxCount(),3)

        h = RFG(fp.name, translate=(-0.5, 0, 0), domainSize=(0.5,0.49,0.49),)
        self.assertEqual(h.fxnet.getFxCount(),1)

        os.remove(fp.name)


if __name__ == '__main__':
    unittest.main()
