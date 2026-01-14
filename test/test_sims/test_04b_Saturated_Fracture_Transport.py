#!/usr/bin/env python
"""Python unittests involving 04b_Saturated_Fracture_Transport
"""

import unittest
import os
from math import prod
import numpy as np

import hgstools.pyhgs as pyhgs
from hgstools.pyhgs._test import skip_if_no_sim_output
from hgstools.pyhgs.mesh import (HGSGrid, Domain)

TESTP = os.path.dirname(__file__)+os.sep
'Path to this testing directory'
SIM_PREFIX = os.path.join(TESTP,'04b_Saturated_Fracture_Transport','module4b')
'Simulation directory+prefix'

@unittest.skipIf(skip_if_no_sim_output(SIM_PREFIX), 'HGS output missing')
class Test_Module04b(unittest.TestCase):

    def setUp(self):
        # Loading and re-loading this each test case could get expensive.
        # Consider refactoring to a "static" module level variable that can be
        # treated as read-only as a balance between data loading overhead and
        # test independence; test cases that modify the HGSGrid object can
        # reload the data.
        self.g = HGSGrid(SIM_PREFIX)

    def test_grid_size(self):
        gl = self.g.get_grid_lines()
        self.assertEqual(tuple(map(len,gl)),(101,2,51))

    def test_grid_indexing(self):
        # Note: In the fracture domain, the 3DNode# is off by one versus the
        # values below, as the HGS indices (Fortran) are 1-based and the HGSGrid
        # indices (Python) are 0-based
        self.assertEqual(self.g.find_node_index(0.,0.,20.), 8080,)
        self.assertEqual(self.g.find_node_index(0.,1.,20.), 8181,)

        self.assertEqual(tuple(self.g.find_grid_index(8080)), (0,0,40))
        self.assertEqual(tuple(self.g.find_grid_index(0.,0.,20.)), (0,0,40))

        self.assertEqual(tuple(self.g.find_grid_index(8181)), (0,1,40))
        self.assertEqual(tuple(self.g.find_grid_index(0.,1.,20.)), (0,1,40))

    def test_calc_flux_weighted_conc(self):
        """Flux-weighted supersample"""

        # alias
        mesh_ssdg = pyhgs.mesh.make_supersample_distance_groups

        # get concentrations
        cnpm = self.g.get_nodal_vals(f'{SIM_PREFIX}o.conc_pm.salt.0010')
        cepm = self.g.get_element_vals(cnpm)
        cefx = self.g.get_element_vals(cnpm,Domain.FRAC)
        del cnpm

        # calculate flux magnitudes (weights)
        qepm = self.g.get_element_vals(f'{SIM_PREFIX}o.q_pm.0001')
        qepm = np.sqrt(np.sum(qepm**2, axis=3))
        qefx = self.g.get_element_vals(f'{SIM_PREFIX}o.v_frac.0001',Domain.FRAC)
        qefx = np.sqrt(np.sum(qefx**2, axis=1))

        # easy, but inefficient implementation
        maxd = (2.,1.,2.)

        # get supersample distance group sizes
        ssdg_pm = []
        for a,maxda in zip(range(3),maxd):
            gla = self.g.get_grid_lines()[a]
            da = gla[1:]-gla[:-1]
            ssdg_pm.append(mesh_ssdg(da,maxda))

        ssdata = np.zeros(prod(len(a) for a in ssdg_pm))

        _iterfunc = self.g.iter_supersample_distance_groups
        for i,grp in enumerate(_iterfunc(
                    maxd, domains=(Domain.PM,Domain.FRAC))):

            # alias
            (ixlo,ixhi),(yl,yh),(zl,zh),elfx = grp

            # guard against 0-length zone --- go with 0.0 already in array
            if any( lo == hi for (lo,hi) in grp[:-1] ): continue

            # pm element grid index slice
            sl = np.s_[ixlo:ixhi,yl:yh,zl:zh]
            c = np.concatenate((cepm[sl].flatten(),cefx[elfx]))
            w = np.concatenate((qepm[sl].flatten(),qefx[elfx]))

            # compute average
            ssdata[i] = np.average(c,weights=w)

        # reshape to ss grid
        ssdata = ssdata.reshape(tuple(len(a) for a in ssdg_pm))

        self.assertTrue( ssdata[48,0,38],
                np.average( np.concatenate(
                        (cepm[48:52,0:1,38:42].flatten('F'),
                         cefx[48:52], cefx[96:100],) ),
                    weights=np.concatenate(
                        (qepm[48:52,0:1,38:42].flatten('F'),
                         qefx[48:52], qefx[96:100],),)
                    )
                )
        

if __name__ == '__main__':
    unittest.main()
