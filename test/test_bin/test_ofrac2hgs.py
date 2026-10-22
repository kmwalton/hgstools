#!/usr/bin/env python
"""Testing of ofrac2hgs module"""

import io
import os
import types
import argparse
import unittest
import tempfile
import contextlib
from unittest import mock

from decimal import Decimal
import hgstools.bin.ofrac2hgs as o2h
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


class TestNudge(unittest.TestCase):
    """Parsing and application of scalar and per-axis nudge increments"""

    @classmethod
    def setUpClass(cls):
        from ofrac.ofracs import OFracGrid

        g = OFracGrid(domainSize=(1.,1.,1.),
              fx=[
                (0.03, 0.97, 0., 1., 0.52, 0.52, 0.001),
                (0.5, 0.5, 0., 1., 0., 0.33, 0.002),
            ],
        )
        fp = tempfile.NamedTemporaryFile(prefix='dfn_', suffix='.pkl',
                delete=False)
        fp.close()
        OFracGrid.pickleTo(g, fp.name)
        cls.pkl = fp.name

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.pkl)

    def test_nudge_arg(self):
        self.assertEqual(nudge_arg('0.5'), Decimal('0.5'))
        self.assertEqual(nudge_arg('10, 10, 0'),
                (Decimal('10'), Decimal('10'), Decimal('0')))
        self.assertEqual(nudge_arg('(1 2 3)'),
                (Decimal('1'), Decimal('2'), Decimal('3')))

        for bad in ['', '1, 2', '1 2 3 4', 'abc', '-1', '1, -1, 1', 'inf',
                'nan']:
            with self.subTest(bad=bad):
                with self.assertRaises(argparse.ArgumentTypeError):
                    nudge_arg(bad)

    def test_nudge_arg_cli(self):
        p = make_arg_parser()
        with self.assertRaises(SystemExit), \
                contextlib.redirect_stderr(io.StringIO()):
            p.parse_args(['-n', 'abc', 'x.pkl'])

    def test_nudge_str(self):
        self.assertEqual(nudge_str(Decimal('0.1')), '0.1')
        self.assertEqual(
            nudge_str((Decimal('0.1'), Decimal('0'), Decimal('0.25'))),
            '0.1, 0, 0.25')

    def _preamble_nudge_lines(self, h):
        s = io.StringIO()
        h.spewPreamble(fout=s)
        return [l for l in s.getvalue().splitlines() if 'nudged' in l]

    def test_scalar(self):
        h = RFG(self.pkl, nudgeTo=Decimal('0.1'))
        self.assertEqual(
            [tuple(map(float, f.d[:6])) for f in h.fxnet.iterFracs()],
            [(0.0, 1.0, 0.0, 1.0, 0.5, 0.5),
             (0.5, 0.5, 0.0, 1.0, 0.0, 0.3)])
        self.assertEqual(self._preamble_nudge_lines(h),
            ['! Fractures and gridlines "nudged" to 0.1 increments.'])

    def test_zero_is_noop(self):
        for z in [Decimal('0'), (Decimal('0'),)*3]:
            with self.subTest(nudgeTo=z):
                h = RFG(self.pkl, nudgeTo=z)
                self.assertFalse(h.isNudged())
                self.assertEqual(
                    [tuple(map(float, f.d[:6])) for f in h.fxnet.iterFracs()],
                    [(0.03, 0.97, 0.0, 1.0, 0.52, 0.52),
                     (0.5, 0.5, 0.0, 1.0, 0.0, 0.33)])
                self.assertEqual(self._preamble_nudge_lines(h), [])

    @unittest.skipUnless(hasattr(o2h.ofracs_module, 'as_nudge_triple'),
            'loaded ofracs does not support per-axis nudging')
    def test_per_axis(self):
        h = RFG(self.pkl, nudgeTo=nudge_arg('0.1, 0, 0.25'))
        self.assertEqual(
            [tuple(map(float, f.d[:6])) for f in h.fxnet.iterFracs()],
            [(0.0, 1.0, 0.0, 1.0, 0.5, 0.5),
             (0.5, 0.5, 0.0, 1.0, 0.0, 0.25)])
        self.assertEqual(self._preamble_nudge_lines(h),
            ['! Fractures and gridlines "nudged" to 0.1, 0, 0.25 increments.'])

    def test_per_axis_guard(self):
        # stand-in for an ofracs that predates per-axis nudging
        old = types.SimpleNamespace(__file__='old/ofracs.py',
                parse=o2h.ofracs_module.parse)
        with mock.patch.object(o2h, 'ofracs_module', old):
            with self.assertRaisesRegex(RuntimeError, 'old/ofracs.py'):
                RFG(self.pkl, nudgeTo=nudge_arg('0.1, 0, 0.25'))

            # guard only bites when per-axis nudging is actually requested
            RFG(self.pkl, nudgeTo=Decimal('0.1'))
            RFG(self.pkl, nudgeTo=nudge_arg('0, 0, 0'))


if __name__ == '__main__':
    unittest.main()
