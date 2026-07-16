#!/usr/bin/env python
"""Unittests for `DischargeCalc.discharge_at` and `.area_at`, run against two
04b datasets."""

import unittest
from unittest import mock
import os
import numpy as np
import numpy.testing as nptest

from hgstools.pyhgs._test import skip_if_no_sim_output
from hgstools.pyhgs.mesh import (HGSGrid, Domain)
from hgstools.pyhgs.calcs import DischargeCalc
from hgstools.pyhgs.parser import peek_NNNN_time

TESTP = os.path.dirname(__file__)
'Path to this testing directory'

# flux output files required by discharge_at
_REQ = ['o.q_pm.0001', 'o.v_frac.0001']

# area_at needs the mesh alone, and must remain usable without the files in
# _REQ; gate its tests on that weaker requirement
_REQ_MESH = ['o.coordinates_pm', 'o.elements_pm',
             'o.coordinates_frac', 'o.elements_frac']


class _OutputWasRead(Exception):
    """Raised by a patched reader to prove no simulated output was touched."""


class _DischargeAtTests:
    """Dataset-agnostic checks for `DischargeCalc.discharge_at`.

    Subclasses set:

    - `SIM_PREFIX` : the simulation path+prefix;
    - `FACE` : a one-element-thick block spanning the full y-z face at the
      x=0 end of the domain (``"x0 x1 y0 y1 z0 z1"``); and
    - `A_PM_EXP` : the expected porous-medium face area (the geometric
      y-extent * z-extent of that face).
    """

    SIM_PREFIX = None
    FACE = None
    A_PM_EXP = None

    def setUp(self):
        self.g = HGSGrid(self.SIM_PREFIX)
        self.calc = DischargeCalc(self.g)

    def test_return_contract(self):
        """3-tuple by default; 4-tuple with with_counts; shared first 3."""
        r3 = self.calc.discharge_at(self.FACE, axis=0, timeidx=1)
        self.assertEqual(len(r3), 3)

        r4 = self.calc.discharge_at(self.FACE, axis=0, timeidx=1,
                                    with_counts=True)
        self.assertEqual(len(r4), 4)

        # the first three returned values must be identical
        self.assertEqual(r3[0], r4[0])
        nptest.assert_array_equal(r3[1], r4[1])
        nptest.assert_array_equal(r3[2], r4[2])

    def test_totals_reconcile(self):
        """Index 0 is the PM+FRAC total for A, Q and counts."""
        t, A, Q, counts = self.calc.discharge_at(
            self.FACE, axis=0, timeidx=1, with_counts=True)
        nptest.assert_allclose(A[0], A[1] + A[2])
        nptest.assert_allclose(Q[0], Q[1] + Q[2])
        self.assertEqual(int(counts[0]), int(counts[1]) + int(counts[2]))

    def test_time_read_from_file(self):
        """`time` is read from the q_pm flux file (not a 0.0 placeholder)."""
        t, A, Q = self.calc.discharge_at(self.FACE, axis=0, timeidx=1)
        self.assertIsNotNone(t)
        self.assertEqual(t, peek_NNNN_time(f'{self.SIM_PREFIX}o.q_pm.0001'))

    def test_pm_face_area(self):
        """PM area equals the geometric y-z area of the face."""
        t, A, Q = self.calc.discharge_at(self.FACE, axis=0, timeidx=1)
        nptest.assert_allclose(A[1], self.A_PM_EXP, rtol=1e-6)

    def test_fracture_area_small_positive(self):
        """Aperture-scale fracture area is positive but << PM area."""
        t, A, Q, counts = self.calc.discharge_at(
            self.FACE, axis=0, timeidx=1, with_counts=True)
        self.assertGreater(int(counts[2]), 0,
            'expected fracture elements intersecting this face')
        self.assertGreater(A[2], 0.0)
        self.assertLess(A[2], A[1])

    def test_axis_str_int_equivalent(self):
        """axis='x' and axis=0 give identical results."""
        ti, Ai, Qi = self.calc.discharge_at(self.FACE, axis=0, timeidx=1)
        ts, As, Qs = self.calc.discharge_at(self.FACE, axis='x', timeidx=1)
        nptest.assert_array_equal(Ai, As)
        nptest.assert_array_equal(Qi, Qs)

    def test_area_at_matches_discharge_at(self):
        """area_at reports exactly the A that discharge_at does."""
        t, A, Q = self.calc.discharge_at(self.FACE, axis=0, timeidx=1)
        nptest.assert_array_equal(self.calc.area_at(self.FACE, axis=0), A)


class _AreaAtTests:
    """Dataset-agnostic checks for `DischargeCalc.area_at`.

    `area_at` is the flux-free half of `discharge_at`: it reports the same face
    areas from the mesh alone. Subclasses set `SIM_PREFIX`, `FACE` and
    `A_PM_EXP` as for `_DischargeAtTests`.
    """

    SIM_PREFIX = None
    FACE = None
    A_PM_EXP = None

    def setUp(self):
        self.g = HGSGrid(self.SIM_PREFIX)
        self.calc = DischargeCalc(self.g)

    def test_totals_reconcile(self):
        """Index 0 is the PM+FRAC total."""
        A = self.calc.area_at(self.FACE, axis=0)
        self.assertEqual(len(A), 3)
        nptest.assert_allclose(A[0], A[1] + A[2])

    def test_pm_face_area(self):
        """PM area equals the geometric y-z area of the face."""
        A = self.calc.area_at(self.FACE, axis=0)
        nptest.assert_allclose(A[1], self.A_PM_EXP, rtol=1e-6)

    def test_fracture_area_small_positive(self):
        """Aperture-scale fracture area is positive but << PM area."""
        A = self.calc.area_at(self.FACE, axis=0)
        self.assertGreater(A[2], 0.0)
        self.assertLess(A[2], A[1])

    def test_axis_str_int_equivalent(self):
        """axis='x' and axis=0 give identical results."""
        nptest.assert_array_equal(self.calc.area_at(self.FACE, axis=0),
                                  self.calc.area_at(self.FACE, axis='x'))

    def test_reads_no_simulated_output(self):
        """area_at touches no output, where discharge_at must."""
        with mock.patch.object(self.calc.sim, 'get_element_vals',
                side_effect=_OutputWasRead('simulated output was read')):

            A = self.calc.area_at(self.FACE, axis=0)
            self.assertGreater(A[1], 0.0)

            # the same patch must still trip discharge_at, or this proves
            # nothing about which reader was avoided
            with self.assertRaises(_OutputWasRead):
                self.calc.discharge_at(self.FACE, axis=0, timeidx=1)

    def test_elem_face_areas_sum_to_area(self):
        """area_at is the sum of the per-element face areas, per domain."""
        ax, elem_length, per_dom = self.calc.elem_face_areas_at(
            self.FACE, axis=0)

        self.assertEqual(ax, 0)
        self.assertGreater(elem_length, 0.0)

        A = self.calc.area_at(self.FACE, axis=0)
        for i, dom in enumerate((Domain.PM, Domain.FRAC)):
            self.assertIn(dom, per_dom)
            elems, Ael = per_dom[dom]
            self.assertEqual(len(elems), len(Ael))
            nptest.assert_allclose(np.sum(Ael), A[i+1])

    def test_invalid_axis_rejected(self):
        """A bad axis raises rather than silently picking one."""
        with self.assertRaises(ValueError):
            self.calc.area_at(self.FACE, axis='w')


@unittest.skipIf(skip_if_no_sim_output(
        os.path.join(TESTP, '04b_Saturated_Fracture_Transport', 'module4b'),
        _REQ), 'HGS output missing')
class Test_DischargeAt_Saturated(_DischargeAtTests, unittest.TestCase):
    SIM_PREFIX = os.path.join(
        TESTP, '04b_Saturated_Fracture_Transport', 'module4b')
    FACE = '0 0.5 0 1 0 25'   # first x-layer, full y-z face
    A_PM_EXP = 1.0 * 25.0


@unittest.skipIf(skip_if_no_sim_output(
        os.path.join(TESTP, '04b_very_coarse_mesh', 'module4b'),
        _REQ), 'HGS output missing')
class Test_DischargeAt_VeryCoarse(_DischargeAtTests, unittest.TestCase):
    SIM_PREFIX = os.path.join(TESTP, '04b_very_coarse_mesh', 'module4b')
    FACE = '0 10 0 1 0 25'    # first x-layer, full y-z face
    A_PM_EXP = 1.0 * 25.0


@unittest.skipIf(skip_if_no_sim_output(
        os.path.join(TESTP, '04b_Saturated_Fracture_Transport', 'module4b'),
        _REQ_MESH), 'HGS mesh missing')
class Test_AreaAt_Saturated(_AreaAtTests, unittest.TestCase):
    SIM_PREFIX = os.path.join(
        TESTP, '04b_Saturated_Fracture_Transport', 'module4b')
    FACE = '0 0.5 0 1 0 25'   # first x-layer, full y-z face
    A_PM_EXP = 1.0 * 25.0


@unittest.skipIf(skip_if_no_sim_output(
        os.path.join(TESTP, '04b_very_coarse_mesh', 'module4b'),
        _REQ_MESH), 'HGS mesh missing')
class Test_AreaAt_VeryCoarse(_AreaAtTests, unittest.TestCase):
    SIM_PREFIX = os.path.join(TESTP, '04b_very_coarse_mesh', 'module4b')
    FACE = '0 10 0 1 0 25'    # first x-layer, full y-z face
    A_PM_EXP = 1.0 * 25.0


if __name__ == '__main__':
    unittest.main()
