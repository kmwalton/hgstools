#!/usr/bin/env python
"""Light unittests for the helpers of bin/hgs-calc-Kbulk.py.

Covers only the parts that need no simulation output: the gradient and Darcy
arithmetic, the zone-bounds check, the head/flux time reconciliation, and the
shape of the exported dict (including the ``--i-only`` case, in which the
flux-derived entries are omitted). The end-to-end behaviour of the CLI is
exercised against real datasets by test_sims/test_discharge.py and by running
the tool itself.
"""

import os
import sys
import json
import subprocess
import unittest
from unittest import mock
import importlib.util

import numpy as np
import numpy.testing as nptest

from hgstools.pyhgs._test import skip_if_no_sim_output
from hgstools.pyhgs.aabbox import AABBox

_BIN = os.path.normpath(os.path.join(
    os.path.dirname(__file__), '..', '..', 'bin', 'hgs-calc-Kbulk.py'))

# 04b_coarse_refined_at_fx holds a single, steady-state output
_SIM_STEADY = os.path.normpath(os.path.join(
    os.path.dirname(__file__), '..', 'test_sims', '04b_coarse_refined_at_fx'))
_SIM_STEADY_PFX = os.path.join(_SIM_STEADY, 'module4b')
_REQ_STEADY = ['o.head_pm.0001', 'o.q_pm.0001', 'o.v_frac.0001']

# zone bounds, on grid lines. The fracture sits at x=25, mid-domain, so LEFT
# and RIGHT split the domain either side of it.
_WHOLE = ('0', '0', '0', '50', '1', '25')
_LEFT = ('0', '0', '0', '25', '1', '25')
_RIGHT = ('25', '0', '0', '50', '1', '25')

_A_FACE_EXP = 1.0 * 25.0    # an x-face spans the full y- and z-extent


def _run_kbulk(*args):
    """Run the tool as a subprocess in the steady sim dir, as the chain does."""
    return subprocess.run([sys.executable, _BIN] + list(args),
                          cwd=_SIM_STEADY, capture_output=True, text=True)


def _load_by_path(path, name):
    """Import a module from `path`; its filename is not a legal module name."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

kb = _load_by_path(_BIN, 'hgs_calc_Kbulk')


class _StubParser:
    """Stands in for `argparse.ArgumentParser`.

    `error` raises instead of writing usage to stderr and exiting, so the
    message can be inspected without the noise of a SystemExit.
    """

    class Error(Exception):
        pass

    def error(self, message):
        raise self.Error(message)


class TestCalcI(unittest.TestCase):
    """`calc_i` divides the face-to-face head drop by the centroid distance."""

    def test_gradient_is_head_drop_over_distance(self):
        zones = [AABBox(0., 0., 0., 10., 1., 1.)]
        distances = [np.array([8., 0., 0.])]
        # per zone: (h_x0, h_x1, h_y0, h_y1, h_z0, h_z1)
        heads = [np.array([10., 2., 0., 0., 0., 0.])]

        i = kb.calc_i(zones, distances, heads, np.array([True, False, False]))

        nptest.assert_allclose(i[0][0], (2.-10.)/8.)

    def test_inactive_axes_left_zero(self):
        zones = [AABBox(0., 0., 0., 10., 1., 1.)]
        distances = [np.array([8., 4., 2.])]
        heads = [np.array([10., 2., 9., 1., 7., 3.])]

        i = kb.calc_i(zones, distances, heads, np.array([True, False, False]))

        self.assertEqual(i[0][1], 0.)
        self.assertEqual(i[0][2], 0.)


class TestCalcKbulk(unittest.TestCase):
    """`calc_Kbulk` applies Darcy's law to the mean of the two face fluxes."""

    def test_kbulk_is_darcy(self):
        zones = [AABBox(0., 0., 0., 10., 1., 1.)]
        i = [np.array([-1e-3, 0., 0.])]
        # per zone: (q_x0, q_x1, q_y0, q_y1, q_z0, q_z1)
        q = [np.array([2e-6, 4e-6, 0., 0., 0., 0.])]

        K = kb.calc_Kbulk(zones, i, q)

        nptest.assert_allclose(K[0][0], -(2e-6 + 4e-6)/(2 * -1e-3))

    def test_zero_gradient_does_not_raise(self):
        """A zero gradient yields inf/nan for the caller to interpret."""
        zones = [AABBox(0., 0., 0., 10., 1., 1.)]
        i = [np.zeros(3)]
        q = [np.array([1e-6, 1e-6, 0., 0., 0., 0.])]

        K = kb.calc_Kbulk(zones, i, q)      # must not raise

        self.assertFalse(np.isfinite(K[0][0]))


class TestAsDict(unittest.TestCase):
    """The exported dict carries the flow entries only when they exist."""

    def setUp(self):
        self.zones = [AABBox(0., 0., 0., 10., 1., 2.)]
        self.A = [np.array([2., 0., 0.])]
        self.Q = [np.array([1., 2., 0., 0., 0., 0.])]
        self.q = [np.array([0.5, 1.5, 0., 0., 0., 0.])]
        self.i = [np.array([-0.25, 0., 0.])]
        self.K = [np.array([4., 0., 0.])]
        self.mask = np.array([True, False, False])

    def _full(self):
        return kb.as_dict(self.zones, self.A, self.Q, self.q, self.i, self.K,
                          self.mask)

    def _i_only(self):
        return kb.as_dict(self.zones, self.A, None, None, self.i, None,
                          self.mask)

    def test_zone_bounds_reported(self):
        self.assertEqual(self._full()[0]['zone'], [0., 0., 0., 10., 1., 2.])

    def test_inactive_axes_omitted(self):
        self.assertEqual(set(self._full()[0]) - {'zone'}, {'x'})

    def test_full_dict_has_flow_entries(self):
        self.assertEqual(set(self._full()[0]['x']),
                         {'A', 'Q0', 'Q1', 'q', 'i', 'Kbulk'})

    def test_q_is_the_mean_of_the_two_faces(self):
        nptest.assert_allclose(self._full()[0]['x']['q'], (0.5 + 1.5)/2)

    def test_i_only_omits_exactly_the_flow_entries(self):
        self.assertEqual(set(self._i_only()[0]['x']), {'A', 'i'})

    def test_i_only_leaves_area_and_gradient_unchanged(self):
        full, ionly = self._full()[0]['x'], self._i_only()[0]['x']
        for k in ('A', 'i'):
            self.assertEqual(full[k], ionly[k])


class TestCheckZonesOnGridLines(unittest.TestCase):
    """Zone bounds must coincide with grid lines; see `iter_layer_bbox`."""

    GL = [np.array([0., 10., 20., 30.]),
          np.array([0., 1.]),
          np.array([0., 5., 10.])]

    def setUp(self):
        self.p = _StubParser()

    def test_bounds_on_grid_lines_pass(self):
        zn = AABBox(0., 0., 0., 30., 1., 10.)
        kb.check_zones_on_grid_lines(self.p, [zn], self.GL)     # must not raise

    def test_bound_inside_an_element_rejected(self):
        zn = AABBox(5., 0., 0., 30., 1., 10.)
        with self.assertRaises(_StubParser.Error) as cm:
            kb.check_zones_on_grid_lines(self.p, [zn], self.GL)

        msg = str(cm.exception)
        self.assertIn('does not lie on a x grid line', msg)
        self.assertIn('between grid lines 0 and 10', msg)
        # the pre-inset trap this check exists to catch
        self.assertIn('NOT be pre-inset', msg)

    def test_bound_outside_the_domain_rejected(self):
        zn = AABBox(0., 0., 0., 40., 1., 10.)
        with self.assertRaises(_StubParser.Error) as cm:
            kb.check_zones_on_grid_lines(self.p, [zn], self.GL)

        self.assertIn('outside the domain', str(cm.exception))

    def test_off_axis_bound_named_correctly(self):
        zn = AABBox(0., 0., 2.5, 30., 1., 10.)
        with self.assertRaises(_StubParser.Error) as cm:
            kb.check_zones_on_grid_lines(self.p, [zn], self.GL)

        self.assertIn('zone bound z0', str(cm.exception))

    def test_float_noise_within_tolerance_passes(self):
        zn = AABBox(1e-9, 0., 0., 30., 1., 10.)
        kb.check_zones_on_grid_lines(self.p, [zn], self.GL)     # must not raise


class TestReadTime(unittest.TestCase):
    """`read_time` reconciles the head and flux files, unless told not to."""

    def test_agreeing_files_return_the_time(self):
        with mock.patch.object(kb, 'peek_NNNN_time', return_value='19358.0'):
            self.assertEqual(kb.read_time('pfx', 9), '19358.0')

    def test_disagreeing_files_raise(self):
        with mock.patch.object(kb, 'peek_NNNN_time',
                               side_effect=['19358.0', '23010.0']):
            with self.assertRaises(RuntimeError) as cm:
                kb.read_time('pfx', 9)

        msg = str(cm.exception)
        self.assertIn('inconsistent across domains', msg)
        self.assertIn('19358.0', msg)
        self.assertIn('23010.0', msg)

    def test_check_flux_false_consults_only_the_head_file(self):
        """--i-only must not require o.q_pm to exist."""
        peek = mock.Mock(return_value='19358.0')
        with mock.patch.object(kb, 'peek_NNNN_time', peek):
            t = kb.read_time('pfx', 9, check_flux=False)

        self.assertEqual(t, '19358.0')
        self.assertEqual(peek.call_count, 1)
        self.assertIn('head_pm', peek.call_args[0][0])

    def test_index_is_zero_padded_to_four(self):
        peek = mock.Mock(return_value='9.0')
        with mock.patch.object(kb, 'peek_NNNN_time', peek):
            kb.read_time('pfx', 9, check_flux=False)

        self.assertIn('pfxo.head_pm.0009', peek.call_args[0][0])


@unittest.skipIf(skip_if_no_sim_output(_SIM_STEADY_PFX, _REQ_STEADY),
                 'HGS output missing')
class TestSteadyFlowEndToEnd(unittest.TestCase):
    """Run the CLI over 04b_coarse_refined_at_fx, whose one output is steady.

    Steady flow through a domain with no internal sources pins the results
    without hard-coding numbers from an earlier run: the same water must cross
    every x-face, each zone's faces must therefore balance (so the tool must
    stay quiet), and Kbulk must follow from q and i by Darcy's law.

    This dataset also reports in seconds, unlike the day-based models the tool
    is usually pointed at.
    """

    @classmethod
    def setUpClass(cls):
        cls.proc = _run_kbulk('-d', '1', '-t', '1', '--json', '-',
                              '-z', *_WHOLE, '-z', *_LEFT, '-z', *_RIGHT)
        if cls.proc.returncode != 0:
            raise RuntimeError(
                f'hgs-calc-Kbulk exited {cls.proc.returncode}:\n{cls.proc.stderr}')
        cls.out = json.loads(cls.proc.stdout)
        cls.zones = cls.out['times'][0]['zones']

    def test_one_time_reported(self):
        self.assertEqual(len(self.out['times']), 1)
        self.assertEqual(self.out['times'][0]['time_index'], 1)

    def test_steady_state_time_reported(self):
        """HGS writes its steady solution at a sentinel time."""
        self.assertEqual(self.out['times'][0]['time'], '1.e20')

    def test_units_read_from_eco(self):
        self.assertEqual(self.out['units'], ['kg', 'm', 's'])

    def test_only_the_requested_axis_reported(self):
        for z in self.zones.values():
            self.assertEqual(set(z) - {'zone'}, {'x'})

    def test_face_area_is_geometric(self):
        for z in self.zones.values():
            nptest.assert_allclose(z['x']['A'], _A_FACE_EXP, rtol=1e-9)

    def test_steady_flux_uniform_across_zones(self):
        """No sources, so the same water crosses every x-face."""
        q = [z['x']['q'] for z in self.zones.values()]
        nptest.assert_allclose(q, q[0], rtol=1e-6)

    def test_kbulk_follows_darcy_from_q_and_i(self):
        for z in self.zones.values():
            x = z['x']
            nptest.assert_allclose(x['Kbulk'], -x['q']/x['i'], rtol=1e-12)

    def test_no_net_discharge_warning_when_balanced(self):
        """Steady flow balances, so the net-discharge check must stay quiet."""
        self.assertNotIn('Net discharge', self.proc.stderr)

    def test_i_only_reproduces_i_and_A(self):
        """--i-only changes what is computed, not what i and A come out as."""
        p = _run_kbulk('-d', '1', '-t', '1', '--i-only', '--json', '-',
                       '-z', *_WHOLE, '-z', *_LEFT, '-z', *_RIGHT)
        self.assertEqual(p.returncode, 0, p.stderr)
        i_zones = json.loads(p.stdout)['times'][0]['zones']

        for key, full in self.zones.items():
            with self.subTest(zone=key):
                ionly = i_zones[key]
                self.assertEqual(set(ionly['x']), {'A', 'i'})
                self.assertEqual(ionly['x']['A'], full['x']['A'])
                self.assertEqual(ionly['x']['i'], full['x']['i'])

    def test_defaults_to_whole_domain_at_first_output(self):
        """No -z and no -t: the whole domain, at index 1."""
        p = _run_kbulk('-d', '1', '--json', '-')
        self.assertEqual(p.returncode, 0, p.stderr)
        out = json.loads(p.stdout)

        self.assertEqual(out['times'][0]['time_index'], 1)
        zones = out['times'][0]['zones']
        self.assertEqual(len(zones), 1)
        nptest.assert_allclose(zones['0']['zone'], [0., 0., 0., 50., 1., 25.])

    def test_off_grid_zone_is_rejected(self):
        """The bounds check reaches the user as a CLI error, not a traceback."""
        p = _run_kbulk('-d', '1', '-t', '1', '-z', '1.25', '0', '0',
                       '50', '1', '25')
        self.assertEqual(p.returncode, 2)
        self.assertIn('does not lie on a x grid line', p.stderr)


if __name__ == '__main__':
    unittest.main()
