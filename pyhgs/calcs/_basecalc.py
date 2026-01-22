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
    def __init__(self, sim, **kwargs):
        """
        Parameters
        ----------
        sim : object
            HGS simulation object. Currently, only `HGSGrid`-type allowed

        kwargs
        ------
        allow_partial : bool
            When searching for nodes and elements based on AABBox-like blocks,
            use this value by default (unless overriden in a method-call
            kwargs).

        """

        # set basic simulation data
        self.sim = sim
        if not isinstance(sim, HGSGrid):
            # this will fail if not a string -- user will have to resolve
            self.sim = HGSGrid(sim)

        self.blocks = {}
        """A dictionary of blocks (sets of nodes and element indicies) and
        per-domain element and node data"""


        self.defaults = dict((k,v) for k,v in kwargs.items())

    def _get_block(self, bs, **kwargs):
        """Choose the nodes and elements for this block and store

        Parameters
        ----------

        kwargs
        ------
        allow_partial : bool
            Nodes and elements on the edges of the block become part of the
            block.

        """

        _key = (bs,)+tuple(kwargs.items())

        if _key not in self.blocks:

            allow_partial = True
            if 'allow_partial' in self.defaults:
                allow_partial = self.defaults['allow_partial']
            elif 'allow_partial' in kwargs:
                allow_partial = kwargs['allow_partial']

            _pl.push()
            _val = {}
            _val['nn'] = 0
            _val['ne'] = 0
            for d in self.sim.domains():
                _pl.push()
                # determine elements
                (el, nd)= self.sim.choose_elements_block(bs, allow_partial, True, d)
                _val[d] = (el,nd)
                _val['nn'] += len(nd)
                _val['ne'] += len(el)
                _pl.pop(f'Processed {d.name} block {bs}')
            self.blocks[_key] = _val
            _pl.pop('Calculated nodes and elements within blocks')

        return self.blocks[_key]
