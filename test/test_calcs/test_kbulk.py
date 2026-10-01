#!/usr/bin/env python
"""Tests for bin/hgs-calc-Kbulk.py, in the x- and z-directions.

Each 04b test simulation directory holds two runs on the same mesh:

``module4b``
    The original problem: head 50 at x=0 and 40 at x=50, so flow is along x.
    (Transport is also simulated, but is not used here.)

``zflow4b``
    Flow only, head 50 at z=0 and 40 at z=25, so flow is along z. Its prefix
    must not start with ``module4b``: `pyhgs.cli.parse_path_to_prefix` would
    otherwise resolve it to the ``batch.pfx`` prefix, ``module4b``.

To create the outputs, run ``grok`` and ``phgs`` for each prefix (set
``batch.pfx`` in turn, then restore it to ``module4b``). The
``04b_coarse_refined_at_fx`` mesh must first be made with
``preprocess-make-mesh.py``. Tests are skipped when outputs are missing.

The domain is 50 x 1 x 25 m with two matrix layers (K=8e-10 m/s for z<12,
1e-9 m/s above) and three fractures of aperture 1e-4 m: horizontal at z=20
(x=0..30) and z=6 (x=20..50), and vertical at x=25 spanning the full height.
"""

import os
import subprocess
import sys
import json
import importlib.util
import logging
import tempfile
import unittest
import warnings

import numpy as np

import hgstools
from hgstools.pyhgs.test import sims_join
from hgstools.pyhgs.mesh import HGSGrid
from hgstools.pyhgs.aabbox import AABBox

KBULK_SCRIPT = os.path.join(os.path.dirname(hgstools.__file__),
        'bin', 'hgs-calc-Kbulk.py')

# load the script as a module to reach its calc_* functions
_spec = importlib.util.spec_from_file_location('hgs_calc_Kbulk', KBULK_SCRIPT)
kbulk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(kbulk)

# the whole domain, x0 y0 z0 x1 y1 z1
ZONE = AABBox(0., 0., 0., 50., 1., 25.)

# fluid properties as reported in the o.eco files (kg-m-s)
RHO, MU, G = 1000., 1.124e-3, 9.80665
APERTURE = 1e-4
K_FRAC = RHO*G*APERTURE**2/(12.*MU)
"""Cubic-law fracture hydraulic conductivity (m/s)"""

K_LOWER, K_UPPER = 8e-10, 1e-9
"""Matrix hydraulic conductivity below and above z=12 (m/s)"""

def setUpModule():
    for _l in (kbulk.logger, kbulk.calcslogger):
        _l.setLevel(logging.WARNING)

def outputs_missing(simdir, pfx):
    return not all(os.path.isfile(sims_join(simdir, f'{pfx}o.{sfx}'))
        for sfx in ('lst', 'head_pm.0001', 'q_pm.0001', 'v_frac.0001',
                    'water_balance.dat'))

def water_balance_inflow(simdir, pfx):
    """Return the steady-state discharge into the ``Head_1`` boundary"""
    with open(sims_join(simdir, f'{pfx}o.water_balance.dat')) as fin:
        lines = fin.read().splitlines()
    _vars = next(l for l in lines if l.startswith('VARIABLES'))
    _vars = [v.strip('"') for v in _vars.split('=',1)[1].split(',')]
    return float(lines[-1].split()[_vars.index('Head_1')])

def calc_kbulk(simdir, pfx, ax):
    """Run the hgs-calc-Kbulk calculation along axis `ax` over `ZONE`

    Returns a dict with the face area `A`, the discharges `Q0` and `Q1` at the
    two opposing faces, the gradient `i` and `Kbulk`.
    """
    grid = HGSGrid(sims_join(simdir, pfx))
    zones = [ZONE,]
    mask = np.zeros(3, dtype=bool)
    mask[ax] = True

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dists = kbulk.calc_distances(grid, zones)
        areas, fluxes = kbulk.calc_q(grid, zones, mask, 1)
        heads = kbulk.calc_avg_heads(grid, zones, mask, 1)
        grads = kbulk.calc_i(zones, dists, heads, mask)
        K = kbulk.calc_Kbulk(zones, grads, fluxes, mask)

    A = areas[0][ax]
    return dict(A=A, Q0=fluxes[0][2*ax]*A, Q1=fluxes[0][2*ax+1]*A,
                i=grads[0][ax], Kbulk=K[0][ax])


class _KbulkTests:
    """Kbulk tests common to the 04b meshes; subclasses set the attributes"""

    SIMDIR = None
    KBULK_X = None
    KBULK_Z = None
    """Expected (regression) values, m/s"""

    RTOL_REGRESSION = 1e-3

    def setUp(self):
        for pfx in ('module4b', 'zflow4b'):
            if outputs_missing(self.SIMDIR, pfx):
                self.skipTest(f'No simulation output for {self.SIMDIR}/{pfx}')

    def check_common(self, r, Qbc):
        # discharge is conserved between the faces and equals the boundary
        # inflow reported by HGS
        np.testing.assert_allclose(r['Q0'], r['Q1'], rtol=1e-6)
        np.testing.assert_allclose(r['Q0'], Qbc, rtol=1e-6)
        # Darcy: flow down-gradient; Kbulk = -q/i
        self.assertLess(r['i'], 0.)
        self.assertGreater(r['Kbulk'], 0.)
        np.testing.assert_allclose(r['Kbulk'], -r['Q0']/r['A']/r['i'],
                rtol=1e-6)

    def test_kbulk_x(self):
        r = calc_kbulk(self.SIMDIR, 'module4b', 0)
        self.assertAlmostEqual(r['A'], 25.)
        self.check_common(r, water_balance_inflow(self.SIMDIR, 'module4b'))

        # Estimate: the connected fracture path (z=20 fracture to x=25, down
        # the vertical fracture to z=6, along to x=50; 64 m long) in parallel
        # with the layered matrix (arithmetic mean K)
        dh, L = 10., 50.
        K_arith = (12.*K_LOWER + 13.*K_UPPER)/25.
        Q_est = dh*K_FRAC*APERTURE/64. + K_arith*25.*dh/L
        np.testing.assert_allclose(r['Q0'], Q_est, rtol=0.05)

        np.testing.assert_allclose(r['Kbulk'], self.KBULK_X,
                rtol=self.RTOL_REGRESSION)

    def test_kbulk_z(self):
        r = calc_kbulk(self.SIMDIR, 'zflow4b', 2)
        self.assertAlmostEqual(r['A'], 50.)
        self.check_common(r, water_balance_inflow(self.SIMDIR, 'zflow4b'))

        # The vertical fracture spans the domain height, in parallel with the
        # matrix layers in series (harmonic mean K). Kbulk is over the PM area.
        K_harm = 25./(12./K_LOWER + 13./K_UPPER)
        K_est = K_FRAC*APERTURE*1./50. + K_harm
        np.testing.assert_allclose(r['Kbulk'], K_est, rtol=0.01)

        np.testing.assert_allclose(r['Kbulk'], self.KBULK_Z,
                rtol=self.RTOL_REGRESSION)

    def test_cli_json_x(self):
        """The script, run in the sim directory, reports the same Kbulk_x"""
        r = calc_kbulk(self.SIMDIR, 'module4b', 0)

        env = dict(os.environ)
        env['PYTHONPATH'] = os.pathsep.join(
            [os.path.dirname(os.path.dirname(hgstools.__file__)),]
            + ([env['PYTHONPATH']] if env.get('PYTHONPATH') else []))

        with tempfile.TemporaryDirectory() as tmpd:
            fn = os.path.join(tmpd, 'Kbulk_x.json')
            subprocess.run(
                [sys.executable, '-W', 'ignore', KBULK_SCRIPT,
                 '-d', '1', '-z', '0', '0', '0', '50', '1', '25',
                 '--json', fn],
                cwd=sims_join(self.SIMDIR), env=env,
                check=True, capture_output=True, timeout=300)
            with open(fn) as fin:
                d = json.load(fin)

        x = d['times'][0]['zones']['0']['x']
        np.testing.assert_allclose(x['Kbulk'], r['Kbulk'], rtol=1e-9)
        np.testing.assert_allclose(x['Q0'], r['Q0'], rtol=1e-9)


class Test_Kbulk_04b_very_coarse_mesh(_KbulkTests, unittest.TestCase):
    SIMDIR = '04b_very_coarse_mesh'
    KBULK_X = 2.5068e-8
    KBULK_Z = 1.5465e-8


class Test_Kbulk_04b_coarse_refined_at_fx(_KbulkTests, unittest.TestCase):
    SIMDIR = '04b_coarse_refined_at_fx'
    KBULK_X = 2.3900e-8
    KBULK_Z = 1.5435e-8


if __name__ == '__main__':
    unittest.main()
