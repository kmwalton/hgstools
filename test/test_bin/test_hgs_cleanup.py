#!/usr/bin/env python
"""Unittests for restart-file categorization (pyhgs.runner.BaseRunner) and
for bin/hgs-cleanup.py.

Restart files are categorized per the HGS reference manual:

- 'input.restart': initial-condition files named in the grok file, e.g., by
  'Restart file for heads' or 'Initial concentration from output file'
  (Sec. 2.6); these must predate the run.
- 'output.restart': restart_*.dat (Sec. 2.5.9), PREFIXo.hen and PREFIXo.cen,
  written by the run for a potential subsequent restart.
- 'suspect.restart': either of the above whose modification time contradicts
  its role.

Modification times are set explicitly with os.utime.
"""

import os
import sys
import subprocess
import tempfile
import unittest
import importlib.util
from pathlib import Path

from hgstools.pyhgs.runner import BaseRunner

_BIN = os.path.normpath(os.path.join(
    os.path.dirname(__file__), '..', '..', 'bin', 'hgs-cleanup.py'))

_GROK = '''\
Restart file for heads
    prev.hen
Initial concentration from output file    ! a comment
    prevo.cen.keep
'''

# run 'sim' starts (first simo.* output) at T0
T0 = 1_700_000_000
_FILES = {
    # name: mtime
    'batch.pfx': T0 - 500,
    'sim.grok': T0 - 400,
    'prev.hen': T0 - 300,
    'prevo.cen.keep': T0 - 300,
    'scratch_grok': T0 - 100,
    'simo.eco': T0,
    'simo.lst': T0 + 600,
    'simo.hen': T0 + 100,
    'simo.cen': T0 + 600,
    'restart_state.dat': T0 + 600,
    'restart_file_info.dat': T0 + 600,
    'restart_state_prev.dat': T0 + 590,
    'restart_file_info_prev.dat': T0 + 590,
}


def _flat(nested, out=None, pfx=''):
    """Ravel categorize_files' nested dict to {file: 'cat.subcat'}"""
    out = {} if out is None else out
    for k, v in nested.items():
        if isinstance(v, dict):
            _flat(v, out, pfx + k + '.')
        else:
            out.update((fn, pfx + k) for fn in v)
    return out


class _SimDirCase(unittest.TestCase):

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.d = Path(self._td.name)
        self.make_files(_FILES)

    def tearDown(self):
        self._td.cleanup()

    def make_files(self, files):
        for fn, mt in files.items():
            p = self.d/fn
            if fn == 'batch.pfx':
                p.write_text('sim\n')
            elif fn == 'sim.grok':
                p.write_text(_GROK)
            else:
                p.write_text('x')
            os.utime(p, (mt, mt))

    def set_mtime(self, fn, mt):
        os.utime(self.d/fn, (mt, mt))

    def categorize(self):
        runner = BaseRunner(None, simdir=str(self.d))
        cat, rest = runner.categorize_files(default_cat='non-hgs')
        return _flat(cat)


class Test_Restart_Categories(_SimDirCase):

    def test_fresh_run(self):
        cats = self.categorize()
        self.assertEqual(cats['prev.hen'], 'input.restart')
        self.assertEqual(cats['prevo.cen.keep'], 'input.restart')
        for fn in ['simo.hen', 'simo.cen', 'restart_state.dat',
                'restart_file_info.dat', 'restart_state_prev.dat',
                'restart_file_info_prev.dat']:
            self.assertEqual(cats[fn], 'output.restart', fn)
        self.assertEqual(cats['sim.grok'], 'input')
        self.assertEqual(cats['simo.lst'], 'output.listing')
        self.assertNotIn('suspect.restart', cats.values())

    def test_input_newer_than_run_is_suspect(self):
        self.set_mtime('prev.hen', T0 + 50)
        with self.assertLogs('hgstools.pyhgs.runner', 'WARNING') as cm:
            cats = self.categorize()
        self.assertEqual(cats['prev.hen'], 'suspect.restart')
        self.assertEqual(cats['prevo.cen.keep'], 'input.restart')
        self.assertTrue(any('prev.hen' in m for m in cm.output))

    def test_output_older_than_run_is_suspect(self):
        self.set_mtime('restart_state.dat', T0 - 50)
        with self.assertLogs('hgstools.pyhgs.runner', 'WARNING'):
            cats = self.categorize()
        self.assertEqual(cats['restart_state.dat'], 'suspect.restart')
        self.assertEqual(cats['restart_file_info.dat'], 'output.restart')

    def test_grok_names_own_final_heads(self):
        # the manual warns against restarting from prefixo.hen in place: the
        # run overwrites its own initial conditions
        (self.d/'sim.grok').write_text(
            'Restart file for heads\nsimo.hen\n')
        with self.assertLogs('hgstools.pyhgs.runner', 'WARNING'):
            cats = self.categorize()
        self.assertEqual(cats['simo.hen'], 'suspect.restart')

    def test_builtin_restart_index(self):
        (self.d/'parallelindx.dat').write_text(
            '__Number_of_CPU\n 1\n__Simulation_Restart\n 2\n')
        with self.assertLogs('hgstools.pyhgs.runner', 'INFO') as cm:
            cats = self.categorize()
        self.assertEqual(cats['restart_state.dat'], 'output.restart')
        self.assertEqual(cats['parallelindx.dat'], 'input')
        self.assertTrue(any('__Simulation_Restart is 2' in m
            for m in cm.output))

    def test_no_reference_time(self):
        for fn in ['simo.eco', 'simo.lst', 'simo.hen', 'simo.cen',
                'scratch_grok']:
            (self.d/fn).unlink()
        self.set_mtime('restart_state.dat', 1)
        cats = self.categorize()
        self.assertEqual(cats['restart_state.dat'], 'output.restart')
        self.assertEqual(cats['prev.hen'], 'input.restart')


def _load_cleanup():
    spec = importlib.util.spec_from_file_location('hgs_cleanup', _BIN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class Test_Category_Matches(unittest.TestCase):

    def test_matches(self):
        m = _load_cleanup().category_matches
        self.assertTrue(m('input', 'input'))
        self.assertTrue(m('input.restart', 'input'))
        self.assertTrue(m('output.tecplot', 'output'))
        self.assertTrue(m('output.tecplot', 'output.tec*'))
        self.assertTrue(m('suspect.restart', 'suspect'))
        self.assertFalse(m('input', 'input.restart'))
        self.assertFalse(m('output.restart', 'input'))
        self.assertFalse(m('inputs', 'input'))


class Test_Cleanup_CLI(_SimDirCase):

    def run_cleanup(self, *args):
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
        return subprocess.run([sys.executable, _BIN, *args], cwd=self.d,
            env=env, capture_output=True, text=True, check=True)

    def test_default_keeps_restart_inputs_and_suspects(self):
        self.set_mtime('restart_state.dat', T0 - 50) # suspect
        self.run_cleanup()
        left = set(os.listdir(self.d))
        self.assertTrue({'prev.hen', 'prevo.cen.keep', 'sim.grok',
            'batch.pfx', 'restart_state.dat'} <= left, left)
        self.assertFalse({'simo.hen', 'simo.cen', 'restart_file_info.dat',
            'restart_state_prev.dat'} & left, left)

    def test_no_keep_inputs_still_keeps_suspects(self):
        self.set_mtime('prev.hen', T0 + 50) # suspect
        self.run_cleanup('--no-keep-inputs')
        left = set(os.listdir(self.d))
        self.assertIn('prev.hen', left)
        self.assertNotIn('prevo.cen.keep', left)

    def test_dry_run_lists_categories(self):
        out = self.run_cleanup('--dry-run')
        txt = out.stdout + out.stderr
        self.assertIn('input.restart', txt)
        self.assertIn('output.restart', txt)
        self.assertEqual(set(os.listdir(self.d)), set(_FILES))


if __name__ == '__main__':
    unittest.main()
