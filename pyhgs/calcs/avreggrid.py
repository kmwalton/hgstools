"""Calculate DFM averages on a regular grid"""

import numpy as np
from scipy.sparse import lil_matrix, csc_matrix

import pyhgs
from pyhgs.mesh import HGSGrid, Domain

# logging stuff
import logging
from pyhgs.calcs._perfstack import _PerfLogStack
logger = logging.getLogger('pyhgs.calcs')
perf_logger = logging.getLogger('pyhgs.calcs_perf')
_pl = _PerfLogStack(perf_logger.info)


class AvRegGrid:
    """Calculate elemental average properties on a regular grid"""

    def __init__(self, grid, reggrid):
        """Initialize a new regular grid

        This includes the process of mapping all HGS grid elements to the new
        regular grid, which can be expensive:
            - elements intersecting each new grid block are determined using
            `choose_elements_block` with partial matches.
        """
        if not isinstance(grid, HGSGrid):
            raise ValueError('HGS mesh must be a rectilinear grid')

        self.g = grid
        reggrid = np.array(reggrid)

        logger.debug(f'Operating on {grid}')


        # determine shape of new grid
        ogl = grid.get_grid_lines()
        nshape = (np.ceil(
            np.array([(ax[-1]-ax[0]) for ax in ogl],)/reggrid)+np.ones(3)
          ).astype(np.int32)
        nsize = np.prod(nshape-[1,1,1])

        logger.debug(f'New discretization has {nshape} {reggrid}-sized blocks')

        # compute new grid cell boundaries
        self.gl = tuple(
                ogl[ax][0] + np.arange(nshape[ax])*reggrid[ax]
                for ax in range(3)
                )

        # store mappings
        _pl.push('Determining mappings...')
        self.e2g = {}
        for dom in grid.domains():
            _map = lil_matrix((grid.get_n_elements(dom), nsize),)
            for ibl, bl in enumerate(self._iter_blocks()):
                elind = grid.choose_elements_block(bl, True, dom=dom)
                for iel in elind:
                    v, pv = grid.intersect(iel, bl, dom, with_proportion=True)
                    if pv > 1e-6:
                        _map[iel,ibl] = pv

            self.e2g[dom] = csc_matrix(_map)

        _pl.pop('Mappings, element index to regular grid,')





    def _iter_blocks(self):
        """Return bounding boxes of successive blocks in the new grid"""
        for xi, xj in zip(self.gl[0][:-1], self.gl[0][1:]):
            for yi, yj in zip(self.gl[1][:-1], self.gl[1][1:]):
                for zi, zj in zip(self.gl[2][:-1], self.gl[2][1:]):
                    yield (xi, xj, yi, yj, zi, zj,)
