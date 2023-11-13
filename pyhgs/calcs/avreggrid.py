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
        self.shape = tuple(len(ax) for ax in self.gl)

        # store mappings
        _pl.push('Determining mappings...')
        self.e2g = {}
        """Mappings of original element indices to new grid cell indices.
        Mappings for PM and FRAC domains are in the form of sparse matricies,
        where rows are the original element indices and columns are the new grid
        indices. Values in the matrix indicate the proportion of volume of the
        original element falls into the new grid cell. (Typical entries are 0
        for no intesection or 1 for completely contained; fractional values mean
        that an element is only partially within the gridcell's volume.
        """
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

    def make_grid_nodes(self, i_fastest=False):
        """Return a sequence of (x,y,z) nodes for a rectilinear grid.

        This returns a 3xN array, where columns are (x,y,z)-triples of node
        locations in an "x-fastest"-ordered grid (which happens to be a
        convenient order for a Tecplot IJK-ordered dataset).  Setting
        `i_fastest` to `True` will reverse this behaviour.

        """
        _pl.push()

        _g = self.gl #alias
        _s = self.shape # alias
        retnodes = np.zeros((3,np.prod(_s),),)

        if i_fastest:
            for xchunk, xi in enumerate(_g[0]):
                ix = xchunk*_s[1]*_s[2]
                jx = ix + _s[1]*_s[2]
                retnodes[0,ix:jx] = xi

                for ychunk, yi in enumerate(_g[1]):
                    iy = ychunk*_s[2]
                    jy = iy+_s[2]
                    retnodes[1,ix+iy:ix+jy] = yi
                    retnodes[2,ix+iy:ix+jy] = _g[2]
        else:
            for zchunk, zi in enumerate(_g[2]):
                nxy = _s[0]*_s[1]
                iz = zchunk*nxy
                jz = iz + nxy
                retnodes[2,iz:jz] = zi

                for ychunk, yi in enumerate(_g[1]):
                    iy = ychunk*_s[0]
                    jy = iy+_s[0]
                    retnodes[1,iz+iy:iz+jy] = yi
                    retnodes[0,iz+iy:iz+jy] = _g[0]

        _pl.pop(f'Making nodes for regular grid')

        return retnodes
