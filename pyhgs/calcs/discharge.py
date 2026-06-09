"""Tool for calculating total discharge in a block of elements representing a face

"""
__docformat__ = 'numpy'

from pathlib import Path
import numpy as np
import tabulate

from ._basecalc import _BaseCalc
from ..mesh import (Domain, HGSGrid,)
from ..aabbox import AABBox

# logging stuff
import logging
from ._perfstack import _PerfLogStack
logger = logging.getLogger('hgstools.pyhgs.calcs')
perf_logger = logging.getLogger('hgstools.pyhgs.calcs_perf')
_pl = _PerfLogStack(perf_logger.info)


class DischargeCalc(_BaseCalc):
    """Tool for calculating discharge at a point in time

        import pyhgs.calcs
        import pyhgs.parser.parse as hgs_parse
        from pyhgs.mesh import HGSGrid

        # get simulation and concentration data
        g = HGSGrid('prefix')

        # set up `DischargeCalc` object
        calc = DischargeCalc(g)


        ..todo..

 
    This class has a memory of *blocks*, groups of nodes and elements, as
    specified in `blockspec` parameters in various methods (anything
    interpretable as an `pyhgs.aabbox.AABBox`). Specifically, these groupings
    are keyed by the exact `blockspec` value as given, so the caller should
    take care that equivalent specifications are not redundant: a `blockspec`
    location value of '0.' and '0', or the string '0 ...' and the tuple
    '(0, ...)', produce an identical block node/element set but are stored as
    distinct (duplicated) *blocks* in this object, thereby foregoing the
    efficiency of reusing a previously determined set of nodes/elements.
    """

    def __init__(self, sim):
        """
        Parameters
        ----------
        sim : `HGSGrid` or `str`
            A `HGSGrid` object or string specifying path/to/prefix of a
            simulation.
        """
        super().__init__(sim)


    def discharge_at(self, blockspec, axis, timeidx=1, with_counts=False):
        """Calculate discharge over the face given by blockspec at given time


        Parameters
        ----------
        blockspec : anything interpretable as an `pyhgs.aabbox.AABBox`
            The block of elements (inclusive of elements bordering this block)
            spanning the face, as an `AABBox`, a ``"x0 x1 y0 y1 z0 z1"``
            string, or a sequence of six bounds. See
            `pyhgs.aabbox.AABBox.from_blockspec`.
        axis : int
            An integer representing the desired axis parallel to the flux and
            discharge: x, y, z are 0, 1 or 2.
        timeidx : int
            The time index of the dataset, used for loading
            `o.flux_pm.00<timeidx>`-like files. Default 1.
        with_counts : bool
            returns the counts of pm and fracture elements used in the
            calculation. Default False

        Returns a tuple with the simulation time, list of areas (total, then
            one per domain), list of summed discharge values (total, then one
            per domain).
        """

        _pl.push('Started calcuting discharge')

        _gl= self.sim.get_grid_lines()

        # make sure that we're dealing with a plane
        bbox = AABBox.from_blockspec(blockspec)
        glidx = np.array(bbox.find_inner_grid_indices(*_gl),
            dtype=np.int32).reshape((3,2))
        #is_plane = (glidx[:,1]-glidx[:,0])<=1
        #if np.sum(is_plane) != 1:
        #    breakpoint()
        #    raise ValueError(f'Not a plane: {blockspec}')
        # normal axis
        #ax = np.argmax(is_plane) ; del is_plane
        
        ax = None
        if axis in (0,1,2):
            ax = axis
        elif axis in ('x','y','z'):
            ax = 'xyz'.index(axis)
        else:
            raise ValueError('Invalid axis')

        axs = 'xyz'[ax]
        _gla = _gl[ax]

        # grid spacing normal to the block
        elem_length = 0.
        igl = glidx[ax][0]
        if igl == len(_gla)-1:
            elem_length = _gla[-1] - _gla[-2]
        else:
            elem_length = _gla[igl+1] - _gla[igl]

        logger.info(f'Using elements at grid[{axs}][{igl:4d}]={_gla[igl]:8.3f} with d{axs}={elem_length:.3f} ')

        def log_plus_minus(dom,Qplus,Qminus,n,A):
            if logger.getEffectiveLevel() <= logging.INFO:
                Qnet = Qplus+Qminus
                logger.info(f'  {str(dom).title():4}({n:3} elems; {A:7.2f} m2) Q{axs}= +{Qplus:.4g} + {Qminus:.4g} = {Qnet:.4g}')

        _doms = list(self.sim.domains())
        _ndom = len(_doms)

        time = 0.
        A = np.zeros(2) # assume only two matter.
        nelem = 0
        Qplus = np.zeros_like(A)
        Qminus = np.zeros_like(A)

        # TODO - this seems to grab fractures on the outside faces of the volume
        # I suppose the criteion is whether all four fracture nodes are also
        # among the PM nodes, which is true ...
        # An additional check would need to be done to see if the centroid of
        # the fracture is within the volume, or on its face
        blk = self._get_block(blockspec, allow_partial=False)

        # get the size of the elements at the face

        if Domain.PM in _doms:

            # aliases
            elems,nodes = blk[Domain.PM]

            # get data
            flux_file = Path(f'{self.sim.prefix}o.q_pm.{timeidx:04d}')
            q = self.sim.get_element_vals(flux_file,'pm')

            # filter to this face
            q = q.reshape((-1,3), order='F')[elems,ax]+0.0
            Ael = self.sim.get_element_volumes('pm').ravel(order='F')[elems]/elem_length

            logger.debug(tabulate.tabulate(zip(elems,q,Ael),headers=['PM Elem',f'q{axs}','A']))

            # calculate totals of area and discharge
            mask = q>=0.0
            A[0] = np.sum(Ael)
            Qplus[0] = np.sum(q[mask]*Ael[mask])
            Qminus[0] = np.sum(q[~mask]*Ael[~mask])
            nelem += len(q)

            # get some extra insight if debugging
            log_plus_minus(Domain.PM, Qplus[0], Qminus[0], len(q), A[0])

        if Domain.FRAC in _doms:
            elems,nodes = blk[Domain.FRAC]
            flux_file = Path(f'{self.sim.prefix}o.v_frac.{timeidx:04d}')
            q = self.sim.get_element_vals(str(flux_file), 'frac')

            q = q[:,ax][elems]+0.0
            Ael = self.sim.get_element_volumes('frac')[elems]/elem_length

            logger.debug(tabulate.tabulate(zip(elems,q,Ael),headers=['Fx Elem',f'q{axs}','A']))

            # calculate totals of area and discharge
            mask = q>=0.0
            A[1] = np.sum(Ael)
            Qplus[1] = np.sum(q[mask]*Ael[mask])
            Qminus[1] = np.sum(q[~mask]*Ael[~mask])
            nelem += len(q)

            # get some extra insight if debugging
            log_plus_minus(Domain.FRAC, Qplus[1], Qminus[1], len(q), A[1])


        # provide warning for unused domains
        _unused_doms = set(_doms) - {Domain.PM, Domain.FRAC}
        if _unused_doms:
            logger.log(logging.WARNING, 'Discharge calculation does not '
                + 'include domains '
                + ', '.join(str(dd) for dd in _unused_doms)
                + f' at {blocspec!s}')


        _pl.pop('Done calcuting discharge')

        # do sums
        A = np.array([A[0], A[0], A[1],]) # A_PM, A_PM, A_Frac
        Q = np.array([sum(Qplus)+sum(Qminus), Qplus[0]+Qminus[0], Qplus[1]+Qminus[1],])

        log_plus_minus('all', sum(Qplus), sum(Qminus), nelem, A[0])

        logger.debug(f'  Calculated Q = {Q[0]:.4g} = PM:{Q[1]:.4g}+FX:{Q[2]:.4g}')

        return time, A, Q
