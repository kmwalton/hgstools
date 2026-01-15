"""Tests reading of HGS binary spatial snapshot output files


Fortran arrays in:
- *.v_pm.0001
- *.v_frac.0001


Depends on test_sims\04b_very_coarse_mesh
"""

import sys
import unittest
from pathlib import Path
import logging
from functools import partial

import numpy as np

from hgstools.pyhgs.test import sims_join
from hgstools.pyhgs.cli import parse_path_to_prefix
from hgstools.pyhgs.parser import parse as hgs_parse
from hgstools.pyhgs.mesh import HGSGrid

TESTDIR = Path(__file__).parents[1]

def _has_very_coarse(fn):
    pfn = Path(sims_join('04b_very_coarse_mesh',fn))
    return pfn.exists()




class ArrayReads(unittest.TestCase):

    def assert_allclose(self, actual, expected, **kwargs):
        """Custom wrapper to clean up the traceback."""
        try:
            np.testing.assert_allclose(actual, expected, **kwargs)
        except AssertionError as e:
            # 'from None' suppresses the "During handling..." section
            sys.tracebacklimit = 0
            raise self.failureException(str(e)) from None

    def _check_very_coarse_q(self, pm_q, fx_q, msg):
        '''Check the velocities of the very coarse domain'''

        def _chk(arr, ind, qx, qy, qz):
            self.assert_allclose(arr[ind], np.array([qx,qy,qz,],dtype=np.float32), rtol=1e-3)

        totest_pm = [
            ( 0, 1.025711737767665e-09, 9.23705569806998262e-23, -2.19769688558635323e-11,),
            ( 1, 1.29296062745254403e-09, -7.81597010896864311e-23, 1.35295413761227223e-10,),
            (11, 6.43397224386887956e-10, 7.10542731441797336e-23, -3.64917734918535075e-10,),
            (19, 9.82371628488465376e-10, -2.36847572940007617e-22, 7.81767359181451127e-11,),
        ]

        totest_fx = [
            (0, 0.00115586072206497192, 1.03322124462373822e-16, 0,),
            (6, 0, 5.42441117031515615e-16, -0.00109142682049423456,),
            (9, 0.00113100302405655384, -4.1328849784949529e-16, 0,),
        ]

        if pm_q is not None:
            if pm_q.ndim > 2:
                pm_q = pm_q.reshape((-1,3), order='F')

            _chk_pm = partial(_chk, pm_q)
            for e in totest_pm:
                with self.subTest(f'{msg} pm flux vals (element {e[0]})'):
                    _chk_pm(*e)

        if fx_q is not None:
            _chk_fx = partial(_chk, fx_q)
            for e in totest_fx:
                with self.subTest(f'{msg} fracture velo, (element {e[0]})'):
                   _chk_fx(*e)
         
    @unittest.skipIf( not _has_very_coarse('module4bo.v_pm.0001'), 'File not found')
    def test_v_flat_element_array(self):

        simdir = Path(sims_join('04b_very_coarse_mesh'))

        _parselogger = logging.getLogger('hgstools.pyhgs.parser')
        _oldlev = _parselogger.getEffectiveLevel()
        _parselogger.setLevel(logging.CRITICAL)

        dpm = hgs_parse(simdir/'module4bo.v_pm.0001')['data']
        dfx = hgs_parse(simdir/'module4bo.v_frac.0001')['data']

        _parselogger.setLevel(_oldlev)

        self._check_very_coarse_q(dpm, dfx, 'raw from file (no-HGSGrid)')

    @unittest.skipIf( not _has_very_coarse('module4bo.v_pm.0001'), 'File not found')
    def test_v_with_hgsgrid(self):

        simdir = Path(sims_join('04b_very_coarse_mesh'))
        g = HGSGrid(simdir)

        dpm = g.get_element_vals(g.ppfx+'o.v_pm.0001','pm')
        dfx = g.get_element_vals(g.ppfx+'o.v_frac.0001','frac')

        self._check_very_coarse_q(dpm, dfx, 'with HGSGrid')


if __name__ == '__main__':
    unittest.main()
