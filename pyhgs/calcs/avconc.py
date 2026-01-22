"""Calculate average properties/values from DFM sub-grids."""

__docformat__ = 'numpy'

import warnings
from itertools import repeat, chain
import numpy as np
from tabulate import tabulate

from ..mesh import Domain
from ._basecalc import _BaseCalc
from ..parser.eco import EcoFile

# logging stuff
import logging
from ._perfstack import _PerfLogStack
logger = logging.getLogger('hgstools.pyhgs.calcs')
perf_logger = logging.getLogger('hgstools.pyhgs.calcs_perf')
_pl = _PerfLogStack(perf_logger.info)


class AvCalc(_BaseCalc):
    """
    Calculate nodal or elemental averages using various weighting schemes.

        import pyhgs.calcs
        import pyhgs.parser.parse as hgs_parse
        from pyhgs.mesh import HGSGrid

        # get simulation and concentration data
        g = HGSGrid('prefix')
        c_fn = g.prefix+'o.c_pm.0001'
        (sim_t,_c_pm_nd) = hgs_parse(c_fn).values()

        # set up `AvCalc` object
        calc = AvCalc(g)

        # Make basic calculations
        blk = '0. 0. 0. 0. 4. 5.' # 1m vertical interval
        a1 = calc.average(blk, _c_pm_nd, weight='arith_nd')
        a2 = calc.average(blk, _c_pm_nd, weight='arith_el')
        a3 = calc.average(blk, _c_pm_nd, weight='vol_el')
        a4 = calc.average(blk, _c_pm_nd, weight='porvol_el')

        # Make calculations with custom weighting
        _c_fx = g.get_element_vals(_c_pm_nd, 'frac')
        _q_pm = g.get_element_vals(g.prefix+'o.q_pm.0001')
        _q_fx = g.get_element_vals(g.prefix+'o.q_frac.0001', 'frac')
        a5 = calc.average(blk, [_c_pm_nd, _c_fx], weight=[_q_pm, _q_fx])


    This class has a memory of *blocks*, groups of nodes and elements, as
    specified in `blockspec` parameters in various methods. Specifically,
    these groupings are keyed by the exact `blockspec` string, so the caller
    should take care with these strings such that the groups of nodes and
    elements are not redundant, like a `blockspec` location value of '0.'
    and '0' will produce an identical block node/element set, but will have
    distinct (duplicately stored) *blocks* in this object.
    """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        sim : `HGSGrid` or `str`
            A `HGSGrid` object or string specifying path/to/prefix of a
            simulation.
        """

        super().__init__(*args, **kwargs)


    def _get_weight(self, blockspec, dom, w_key, evalstr=None, execstr=None):
        """Return an array of the given weight values

        If the weight does not exist yet, create it using evalstr and add to the
        self.block library.

        evalstr is a "simple" string that will give the array of values via
        eval()

        execstr is a more complex set of python statements that assigns the
        final array of weight values to the local variable 'retarr'.
        
        """
        bl = self._get_block(blockspec)
        el, nd = bl[dom]

        if w_key not in bl:
            bl[w_key] = {}

        if dom not in bl[w_key]:

            if evalstr:
                bl[w_key][dom] = eval(evalstr)
            else:
                _l = locals()
                _l['retarr'] = None
                exec(execstr, globals(), _l)
                bl[w_key][dom] = _l['retarr']

        return bl[w_key][dom]

    def _get_volumes(self, blockspec, dom):
        """Return the volumes of the elements"""
        return self._get_weight(blockspec, dom,
            'w_vol',
            "self.sim.get_element_volumes(dom).ravel('F')[el]",
            )

    def _get_phis(self, blockspec, dom):
        """Return the volumes of the elements"""
        return self._get_weight(blockspec, dom,
                'w_phi',
                "self.sim.get_elements_data(dom)['porosity'].ravel('F')[el]",)

    def _get_el_weights(self, blockspec, dom):
        return self._get_weight(blockspec, dom,
                'w_el',
                execstr="""
with warnings.catch_warnings():
    #targeting the warning of the 0-based fracture node indices
    warnings.simplefilter('ignore')
    _inc = self.sim.get_elements_data(dom)['inc']

# get the weight of each element (the count of times it intersects
# blocspec)
retarr = np.zeros(len(el), dtype=np.single)
_ndset = set(nd)
for i,iel in enumerate(el):
    _incset = set(_inc[iel])
    retarr[i] = len(_ndset & _incset)
# divide out the number of nodes per element
retarr /= self.sim.get_elements_data(dom)['nln']
""",
            )

    def _get_data(self, dataspec, doms):
        """Return a list of data arrays corresponding to `doms`"""

        if isinstance(dataspec, (tuple, list,)) and len(dataspec) == len(doms):
            return dataspec

        elif isinstance(dataspec,str):
            v_nodal = self.sim.get_nodal_vals(dataspec)
            return [ self.sim.get_nodal_vals(v_nodal, dom=d) for d in doms ]

        elif isinstance(dataspec,np.ndarray):
            return [ self.sim.get_nodal_vals(dataspec, dom=d) for d in doms ]

        else:
            return ValueError(f'Unexpected type of data {type(dataspec)}')

    WEIGHTS=['arith_nd', 'arith_el', 'vol_el', 'porvol_el', ]
    """Weighting methods understood by `average`"""

    def average(self, blockspec, data, weight='arith_el', doms='all'):
        """Calculate the average concentration of the block within data

        Parameters
        ----------
        blockspec : str
            A string representing the block of elements (inclusive of elements
            bordering this block.
        data : `list`, array, or `str`
            The 1D or 3D data array, string representing a HGS binary
            concentration file name, or a floating point value indicating the
            desired simulation time (the closest time will be used).<br>
            __What if there is different data for different domains?__
        weight : `list`, array, or `str`
            The weighting scheme, by keyword string or an array|`str` that is
            interpreted in the same way as the `data` parameter:<br>
            - `arith_nd`: arithmetic mean of the nodes within blockspec
            - `arith_el': arithmetic mean of the elements within blockspec
            - 'vol_el': weight by element volume 
            - 'porvol_el': weight by element pore volume (V * phi). In the case
            of concentration, this assumes fully saturated conditions
            list
            - `list`, array, or `str`: interpret arrays or strings in the same way
            as the `data` parameter is iterpreted, or pass a `list` where entries
            in the list are 1) in the same sequence as the simulation's
            `domain()`s, and 2) arrays of weight values for every element in the
            domain. This allows a very flexible specificaion of weights, so the
            caller can pre-calculate flux magnitude, horizontal flux magnitude,
            etc., to be used as weight values. The data in passed weight arrays
            will be multiplied with the "inherent weight" of an element, which
            is the number of incident nodes it has within `blockspec`.
        doms : `str`, `list`, `pyhgs.mesh.Domain`
            A list of `pyhgs.mesh.Domain` objects, or a single `Domain`, or a
            string representing a single `Domain`, or the keyword 'all' (for all
            domains)

        """

        _doms = self._doms_list(doms)

        _bldata = self._get_block(blockspec)

        _data = self._get_data(data, _doms) # interperet

        # prepare for caller's weighting quantities
        _has_weight_in = False
        _weight_in = [None, None] # dummy, will be inored
        if type(weight) == list:
            _has_weight_in = True
            # interperet
            _weight_in = self._get_data(weight, _doms)
            # re-interperet to shape for element-wise lookup
            for i,ww in enumerate(_weight_in):
                _weight_in[i] = ww.flatten('F')
            
        # work arrays, allocate
        v = None
        w = None
        if _has_weight_in or 'el' in weight:
            v = np.zeros(sum(len(_bldata[d][0]) for d in _doms))
        elif 'nd' in weight:
            v = np.zeros(sum(len(_bldata[d][1]) for d in _doms))
        else:
            raise ValueError(f'nd or el-type not specified in weight {weight}')

        if _has_weight_in or 'arith' not in weight:
            w = np.ones(v.size)

        _logging_intermediates = dict()

        # get data from domains
        _iwork = 0 # index into work arrays
        for idom,(d,ddata,wdata) in enumerate(zip(_doms, _data, _weight_in)):
            _nv = 0
            _logging_intermediates[d]=dict()

            (el, nd) = _bldata[d]

            # keyword weighting types
            if weight == 'arith_el':
                _nv = len(el)
                v[_iwork:_iwork+_nv] = \
                    self.sim.get_element_vals(ddata, dom=d).ravel('F')[el]

            elif weight == 'arith_nd':
                _nv = len(nd)
                v[_iwork:_iwork+_nv] = \
                    self.sim.get_nodal_vals(ddata, dom=d).ravel('F')[nd]

            elif weight in ('vol_el','porvol_el'):
                _nv = len(el)
                v[_iwork:_iwork+_nv] = \
                    self.sim.get_element_vals(ddata, dom=d).ravel('F')[el]

                _w = self._get_el_weights(blockspec, d)
                _V = self._get_volumes(blockspec, d)
                w[_iwork:_iwork+_nv] = _w*_V

                _logging_intermediates[d]['w']=_w
                _logging_intermediates[d]['Vol']=_V

                if weight == 'porvol_el':
                    _phi = self._get_phis(blockspec, d)
                    _logging_intermediates[d]['phi']=_phi
                    w[_iwork:_iwork+_nv] *= _phi

            elif wdata is not None:
                try:
                    _nv = len(el)
                    v[_iwork:_iwork+_nv] = \
                        self.sim.get_element_vals(ddata, dom=d).ravel('F')[el]

                    _w = self._get_el_weights(blockspec, d)
                    _logging_intermediates[d]['w']=_w
                    _logging_intermediates[d]['user']=wdata[el]
                    w[_iwork:_iwork+_nv] = _w*wdata[el]

                except Exception as e:
                    #breakpoint()
                    print('found an error, try p str(e)')

            else:
                raise ValueError(f'Unknown weighting {weight}')

            _iwork += _nv

        if logger.level <= logging.DEBUG:

            # determine node/element array size and output
            _ien = 0
            _iname = 'element'
            if _has_weight_in or 'el' in weight:
                pass
            elif 'nd' in weight:
                _ien = 1
                _iname = 'node'

            tabulate_kwargs = {}
            tabulate_kwargs['headers'] = ['Dom', _iname, 'Value', 'Weight']
            _hxtra = list(_logging_intermediates[_doms[0]].keys())
            tabulate_kwargs['headers'].extend(_hxtra)

            _iter_dom_names = chain.from_iterable(
                    repeat(d.name, len(_bldata[d][_ien])) for d in _doms)
            _iter_indices = chain.from_iterable(
                    _bldata[d][_ien] for d in _doms)
            _iter_intermed_vals = lambda k: chain.from_iterable(
                    _logging_intermediates[d][k] for d in _doms)
            _iter_w = repeat(1.0)
            if w is not None:
                _iter_w = iter(w)

            #breakpoint()
            s = f'Raw data for averaging "{blockspec}":\n'
            s += tabulate(zip(
                _iter_dom_names,
                _iter_indices,
                v,
                _iter_w,
                *( _iter_intermed_vals(k) for k in _hxtra),
                ),
                **tabulate_kwargs,
            )

     
            logger.debug(s)

        # perform weighting (or not) and return
        if w is None:
            return np.average(v)
        return np.average(v,weights=w)

    def _doms_list(self, d):
        if isinstance(d, str):
            if d.lower() == 'all':
                return sorted(list(self.sim.domains()))
            return [Domain.a2D(d),]
        elif isinstance(d, Domain):
            # Domain
            return [d,]
        else:
            # listlike
            return sorted(list(Domain.a2D(dd) for dd in d))

def avconc(sim, blockspec, data, weights=None):
    """*Work in progress*

    A method to get a "quick" average calculation without creating an `AvCalc`
    object.

    Parameters
    ----------
    sim : `pyhgs.mesh.HGSGrid`
        The mesh. (Currently, limited to HGSGrid.)

    data : str or numpy.ndarray
        The concentration data, as either a filename or array of data.


    """

    raise NotImplementedError()

    el_nd_ind = []
    for d in g.domains():
        (el, nd)= g.choose_elements_block(blockspec, True, True, d)
        el_nd_ind.append( (el, nd,) )

    return 0.

