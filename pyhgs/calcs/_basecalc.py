"""Base class for various calculators


Contains useful common code for storing sets of nodes and elements
"""

# logging stuff
import logging
from ._perfstack import _PerfLogStack
logger = logging.getLogger('hgstools.pyhgs.calcs')
perf_logger = logging.getLogger('hgstools.pyhgs.calcs_perf')
_pl = _PerfLogStack(perf_logger.info)


from ..mesh import HGSGrid, Domain

class _BaseCalc():
    def __init__(self, sim):
        # set basic simulation data
        self.sim = sim
        """HGS simulation object. Currently, only `HGSGrid`-type allowed"""
        if not isinstance(sim, HGSGrid):
            # this will fail if not a string -- user will have to resolve
            self.sim = HGSGrid(sim)

        self.blocks = {}
        """A dictionary of blocks (sets of nodes and element indicies) and
        per-domain element and node data"""

    def _get_block(self, bs):
        """Choose the nodes and elements for this block and store"""

        if bs not in self.blocks:
            _pl.push()
            self.blocks[bs] = {}
            self.blocks[bs]['nn'] = 0
            self.blocks[bs]['ne'] = 0
            for d in self.sim.domains():
                _pl.push()
                # determine elements
                (el, nd)= self.sim.choose_elements_block(bs, True, True, d)
                self.blocks[bs][d] = (el,nd)
                self.blocks[bs]['nn'] += len(nd)
                self.blocks[bs]['ne'] += len(el)
                _pl.pop(f'Processed {d.name} block {bs}')
            _pl.pop('Calculated nodes and elements within blocks')

        return self.blocks[bs]
