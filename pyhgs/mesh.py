"""Inspect Hydrogeosphere meshes; interperet binary files.


*Example 1:* Get data from your HGS simulations in as easy as three lines...

    # Imports
    from pyhgs.mesh import HGSGrid, Domain

    # get the required grid information
    grid = HGSGrid('pfx')

    # Read data!
    head_n = grid.get_nodal_vals('pfxo.head_pm.0001')
    conc_fx_el = grid.get_element_vals('pfxo.conc_pm.salt.0009', Domain.FRAC)

    # Find out which fracture node is incident with a pm node
    # Fracture nodes are indexed starting at zero (whereas the orginal datafile
    # o.coordinates_frac starts indexing at the number of PM nodes + 1).

Notes
-----

- Log messages (for a smattering of functions herein) are available using the
    standard `logging` interface. For example, to see `DEBUG`-level messages,
    in the console  use:
    
        meshlogger = logging.getLogger('pyhgs.mesh')
        meshlogger.addHandler(logging.StreamHandler())
        meshlogger.setLevel(logging.DEBUG)

    Logging output can be distinguished from the caller's log output using
    `Formatter`:

        meshlogger = logging.getLogger('pyhgs.mesh')
        meshlogger.setLevel(logging.DEBUG)
        meshhandler = logging.StreamHandler()
        meshhandler.setFormatter(logging.Formatter('HGSGrid - %(message)s'))
        meshlogger.addHandler(meshhandler)

Future Development Notes
------------------------
Consider creating pyghs.sim.HGSSim to unite pyghs.mesh-objects with physical
properties and temporal data. (E.g., _PropFinderDict would be moved to this
HGSSim object.) This would better encapsulate mesh considerations and separate
them from static and temporal data of the system being simulated.

"""

import os
import re
import time
from itertools import count,repeat,product,chain,pairwise
from collections import (defaultdict, namedtuple, deque)
from bisect import bisect,bisect_left,bisect_right
import decimal
import glob
from decimal import Decimal
from enum import IntEnum,auto
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor

from functools import partial

import warnings

import numpy as np
from numpy.linalg import norm
from scipy.sparse import lil_array

from .parser import (
    parse,
    parse_coordinates_pm,
    parse_elements_pm,
    parse_coordinates_frac,
    parse_elements_frac,
    HGS_DATATYPE,
    get_datatype,
    )
from .parser.mprops import CaseInsensitiveDict

import logging
logger = logging.getLogger(__name__)
perf_logger = logging.getLogger(__name__+'_perf')

# performance
class _PerfLogQ:
    def __init__(self, logfunc):
        self.perfq = deque()
        self.logfunc = logfunc
    def push(self, msg=None):
        if msg is not None:
            self.logfunc(msg)
        self.perfq.append(time.perf_counter())
    def pop(self, msg):
        self.logfunc(
            msg
            +f' took {time.perf_counter()-self.perfq.pop():.2f}s')
        return self

_plq = _PerfLogQ(perf_logger.info)



__docformat__ = 'numpy'

decimal.getcontext().prec = 12
N_COORD_DIG = Decimal('0.001')
N_APERT_DIG = Decimal('0.000001')

def D(v,new_prec):
    """Return a decimal with the specified quantization/precision"""

    try:
        return Decimal(v).quantize(new_prec)
    except decimal.InvalidOperation as e:
        if not v.is_finite():
            return v
        else:
            raise ValueError(f'Cannot re-quantize {v}') from e
    except Exception as e:
        raise ValueError(f'Argument {v} of type {type(v)}.') from e

def D_CO(v):
    """Return a decimal with the required number of digits for coordinates"""
    return D(v,N_COORD_DIG)

def D_AP(v):
    """Return a decimal with the required number of digits for apertures"""
    return D(v,N_APERT_DIG)

def toDTuple(s):
    """Return a tuple of coordinate precision Decimals

    Parameters
    ----------
    s : `str` or list-like
        strs must look like a list of numbers tuple, like '(x,y,z)',
        'x y z', 'x,y z', or '[u, v, w x y z'
        list-likes must be a sequence of number-like things.
    """

    ivals = iter([0,0,0,])

    if type(s) is str:
        ivals = iter(re.sub("[(),]",' ',s).strip().split())
    elif hasattr(s,'__iter__'):
        ivals = iter(s)
    else:
        raise ValueError(f'"toDTuple" expected string or list-like, but got '
                '{type(s)}.')

    return tuple( D_CO(v) for v in ivals )

# helpful type
_BoundingBox = namedtuple('_BoundingBox', 'x1 x2 y1 y2 z1 z2'.split(), )

class Domain(IntEnum):
    """Types of domains in Hydrogeosphere"""
    PM = auto()
    FRAC = auto()
    NODES_PM = auto()

    @staticmethod
    def a2D(a):
        """Convert "any" object to a Domain

        ...where "any" can be a string, or a `Domain`

        Raises:
            ValueError : When a Domain cannot be build from `a`.
        """
        if type(a) == Domain:
            return a
        elif type(a) == str:
            try:
                return Domain[a.upper()]
            except KeyError:
                pass # flow through to final 'raise'

        raise ValueError(f'Cannot build a Domain from {a}')

    def __str__(self):
        """Return a representation consistent with HGS output file name part"""
        return self.name.lower()

class _WarningDict(dict):
    def __init__(self, *args, **kwargs):
        self.warn_on = dict()
        """
        dict of key:param_list that will produce warnings when key is retrieved

        param_list must be valid argument to the call
        `warnings.warn(*param_list)`
        """
        super().__init__(*args, **kwargs)

    def __getitem__(self, k):
        if k in self.warn_on:
            warnings.warn(*self.warn_on[k])
        return super().__getitem__(k)



class HGSGrid():
    '''Inspect a Hydrogeosphere rectilinear grid'''

    def __init__(self, prefix):

        self.prefix = prefix
        """Simulation path and prefix"""
        self.hgs_pm_nodes = None
        """Simulation PM node data from
        `hgstools.pyhgs.parser.parse_coordinates_pm`
        """
        self.hgs_pm_elems = None
        """Simulation PM node data from
        `hgstools.pyhgs.parser.parse_elements_pm`
        """
        self.hgs_fx_nodes = None
        """Simulation fracture node data from
        `hgstools.pyhgs.parser.parse_coordinates_frac`
        """
        self.hgs_fx_elems = None
        """Simulation fracture node data from 
        `hgstools.pyhgs.parser.parse_elements_frac`
        """

        self.nn = 0
        """Number of PM nodes"""
        self.ne = 0
        """Number of PM elements"""
        self.nfn = 0
        """Number of fracture nodes"""
        self.nfe = 0
        """Number of fracture elements"""

        def _my_call(*args,**kwargs):
            return args[0](*args[1:],**kwargs)
        def _my_call_2(a):
            return a[0](a[1:])

        self.hgs_pm_nodes = parse_coordinates_pm(prefix)
        self.hgs_pm_elems = parse_elements_pm(prefix)
        self.nn = self.hgs_pm_nodes['nn'] # aliases
        self.ne = self.hgs_pm_elems['ne']

        if os.path.exists(f'{prefix}o.coordinates_frac'):
            self.hgs_fx_nodes = parse_coordinates_frac(prefix)
            self.hgs_fx_elems = parse_elements_frac(prefix)
            self.nfn = self.hgs_fx_nodes['nnfrac']
            self.nfe = self.hgs_fx_elems['nfe']

            # ..Notes::
            #
            # hgs_fx_nodes['link_pm2frac'] holds index values > the number of pm
            # nodes, i.e., the "global ordering" range.
            #
            self.hgs_fx_nodes['link_pm2frac'] -= self.nn

            # hgs_fx_elems['inc'] holds index values of PM nodes. Modify this to
            # hold fracture node indices
            _link = self.hgs_fx_nodes['link_pm2frac'] #alias
            _new_inc = np.fromiter(
                (_link[i] for i in self.hgs_fx_elems['inc'].flatten()),
                count=self.nfe*self.hgs_fx_elems['nln'],
                dtype=np.int32,
                )
            self.hgs_fx_elems['inc'] = \
                _new_inc.reshape(-1,self.hgs_fx_elems['nln'])

            # turn on warnings
            self.hgs_fx_nodes = _WarningDict(self.hgs_fx_nodes)
            self.hgs_fx_nodes.warn_on['link_pm2frac'] = [
                '0-based, fracture node indicies held here', UserWarning, 2,]

            self.hgs_fx_elems = _WarningDict(self.hgs_fx_elems)
            self.hgs_fx_elems.warn_on['inc'] = [
                '0-based, fracture node indicies held here', UserWarning, 2,]


        self.shape = tuple(self.hgs_pm_nodes[a] for a in ['nx','ny','nz'])
        """Shape of the PM node grid"""
        self.elshape = tuple(self.hgs_pm_nodes[a]-1 for a in ['nx','ny','nz'])
        """Shape of the PM element grid"""

        if self.hgs_pm_nodes['tetramesh'] or self.hgs_pm_nodes['nb2d']>0:
            raise ValueError(
                f'Simulation {prefix} does not have a rectilinear grid')

        self._n2e = dict()
        """dict with eapping of node index to element index for each domain"""

        logger.debug(f'Initialized HGSGrid with {self.shape} PM nodes')

    def __str__(self):
        ndoms = len(list(self.domains()))
        s = f'HGSGrid {self.prefix} with {ndoms} domain{"s" if ndoms>1 else ""}'
        s += f' and {self.nn:,} PM nodes'
        return s

    def domains(self):
        """Yield the active domains in this grid"""
        if self.hgs_pm_nodes is not None:
            yield Domain.PM
        if self.hgs_fx_nodes is not None:
            yield Domain.FRAC

    @staticmethod
    def _to_index_nocheck(t,shp):
        # Recursive call is supremely slow -- only good if 2D/3D grid is
        # unknown.
        #if len(t) == 1 : return t[0]
        #return t[0] + shp[0]*HGSGrid._to_index_nocheck(t[1:],shp[1:])

        # 3D grid calculation
        return t[0] + shp[0]*(t[1] + shp[1]*t[2])

    @staticmethod
    def _to_index(t,shp):
        """Return grid index as linear index"""
        if any(t[i] >= shp[i] for i in range(len(shp))):
            s = tuple( v-1 for v in shp )
            raise ValueError(
                f'Element index {t} out of 3D grid index bounds. Need <={s}.')

        if any(v<0 for v in t):
            raise ValueError(
                'Element index {t} out of grid index bounds. '\
                f'Need >={len(shp)*(0,)}.')

        return HGSGrid._to_index_nocheck(t,shp)

    @staticmethod
    def _to_grid_index(iin,shp):
        i = 3*[-1,]
        i[0] = iin % shp[0]
        i[1] = int((iin - i[0])/shp[0]) % shp[1]
        i[2] = int((iin - i[0] - shp[0]*i[1]) /(shp[0]*shp[1]))
        return i

    def get_grid_lines(self):
        """Return the ([x-],[y-],[z-grid lines]) in this rectilinear grid.

        If the grid is not rectilinear, one or more sets of gridlines will be
        returned as `None`.
        """
        # store
        if not hasattr(self,'gl'):
            self.gl = determine_grid_lines(self.hgs_pm_nodes)

        return self.gl


    def get_n_elements(self,dom=Domain.PM):
        """Return the number of elements in the specified domain"""
        dom = Domain.a2D(dom)

        n = -1

        if dom == Domain.PM:
            n = self.hgs_pm_elems['ne']
        elif dom == Domain.FRAC:
            n = self.hgs_fx_elems['nfe']
        else:
            raise NotImplementedError(f'Not implemented for {dom}')

        return n


    def iter_elements(self,dom=Domain.PM):
        """Iterate over element data.

            Yields tuples of
            <code>( (node indices), zone, *<other data\\>* )</code>

            where <code>*<otherdata\\>*</code>is:
                for <code>dom == Domain.PM</code>, nothing;
                for <code>dom == Domain.FRAC</code>, ap.
        """
        dom = Domain.a2D(dom)

        if dom == Domain.PM:
            raise NotImplementedError()

        elif dom == Domain.FRAC:
            with warnings.catch_warnings():
                #targeting the warning of the 0-based fracture node indices
                warnings.simplefilter('ignore')
                _inc = self.hgs_fx_elems['inc']

            for inc,ap,zn in zip(
                _inc,
                self.hgs_fx_elems['ap'],
                self.hgs_fx_elems['zone'] ):

                yield (inc, zn, ap,)

        else:
            raise NotImplementedError(f'Not implemented for {dom}')


    def get_nodes_data(self, dom):
        """Transparently retrieve the `hgs_X_nodes` dictionary, for domain X"""
        dom = Domain.a2D(dom)

        if dom == Domain.PM:
            return self.hgs_pm_nodes
        elif dom == Domain.FRAC:
            return self.hgs_fx_nodes
        else:
            raise ValueError(f'Domain {dom} not allowed')

    def get_elements_data(self, dom):
        """Transparently retrieve the `hgs_X_elems` dictionary, for domain X"""
        dom = Domain.a2D(dom)

        if dom == Domain.PM:
            if not isinstance(self.hgs_pm_elems, _PMPropFinderDict):
                self.hgs_pm_elems = _PMPropFinderDict(self, self.hgs_pm_elems)
            return self.hgs_pm_elems
        elif dom == Domain.FRAC:
            if not isinstance(self.hgs_fx_elems, _FxPropFinderDict):
                self.hgs_fx_elems = _FxPropFinderDict(self, self.hgs_fx_elems)
            return self.hgs_fx_elems
        else:
            raise ValueError(f'Domain {dom} not allowed')


    # TODO
    # Write a get_coords method that iterates over sets of nodal coordinates, in
    # contrast to the following method that gets one set at a time.

    def get_coords(self, iel, dom=Domain.PM):
        """
            Parameters
            ----------
            iel : int or tuple
                Scalar index, or
                Grid index triple, (ix, iy, iz), for appropriate
                domain types.

            Returns
            -------
            Nodal coordinates (x0,y0,z0), or
            Tuple of bounding nodal coordinates
                `( (x0,y0,z0), (x1,y1,z1), ... )` that bound the element
        """
        dom = Domain.a2D(dom)

        if dom == Domain.PM:
            if type(iel) == tuple:
                iel = self._to_index_nocheck(iel,self.elshape)
            inc_nodes = self.hgs_pm_elems['inc'][iel]
            return tuple(self.hgs_pm_nodes['ncoords'][i] for i in inc_nodes)

        elif dom == Domain.NODES_PM:
            if type(iel) == tuple:
                iel = self._to_index_nocheck(iel,self.shape)
            return self.hgs_pm_nodes['ncoords'][iel]

        elif dom == Domain.FRAC:
            inc_nodes = self.hgs_fx_elems['inc'][iel]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _f2p = self.hgs_fx_nodes['link_frac2pm']
            return tuple(self.hgs_pm_nodes['ncoords'][_f2p[inc_nodes]])

        else:
            raise NotImplementedError(f'Not implemented for {dom}')
        

    #
    # Several helper methods to determine the type of data in an array and its
    # domain
    #
    def _is_FRAC_elemental_scalar(self, dd):
        return self._is_FRAC_scalar(dd,self.nfe)
    def _is_FRAC_elemental_vector(self, dd):
        return self._is_FRAC_vector(dd,self.nfe)
    def _is_FRAC_elemental(self, dd):
        return self._is_FRAC_elemental_scalar(dd) \
            or self._is_PM_elemental_vector(dd)
    def _is_FRAC_nodal_scalar(self, dd):
        return self._is_FRAC_scalar(dd,self.nfn)
    def _is_FRAC_nodal_vector(self, dd):
        return self._is_FRAC_vector(dd,self.nfn)
    def _is_FRAC_nodal(self, dd):
        return self._is_FRAC_nodal_scalar(dd) \
            or self._is_FRAC_nodal_vector(dd)
    def _is_FRAC_scalar(self, dd, n):
        return dd.ndim == 1 and dd.size == n # 1D, scalar
    def _is_FRAC_vector(self, dd, n):
        return dd.ndim == 2 and dd.shape[-1] == 3 and dd.size == 3*n # 1D, vector

    def _is_PM_elemental(self, dd):
        return self._is_PM_elemental_scalar(dd) \
            or self._is_PM_elemental_vector(dd)
    def _is_PM_elemental_scalar(self, dd):
        return self._is_PM_scalar(dd,self.ne)
    def _is_PM_elemental_vector(self, dd):
        return self._is_PM_vector(dd,self.ne)
    def _is_PM_grid(self, dd):
        return dd.ndim > 2
    def _is_PM_nodal(self, dd):
        return self._is_PM_nodal_scalar(dd) \
            or self._is_PM_nodal_vector(dd)
    def _is_PM_nodal_scalar(self, dd):
        return self._is_PM_scalar(dd,self.nn)
    def _is_PM_nodal_vector(self, dd):
        return self._is_PM_vector(dd,self.nn)
    def _is_PM_scalar(self, dd, n):
        return any([dd.ndim == 3 and dd.size == n, # 3D, scalar
                dd.ndim == 1 and dd.size == n, # 1D, scalar
        ])
    def _is_PM_vector(self, dd, n):
        return any([
            dd.ndim == 4 and dd.shape[-1] == 3 and dd.size == 3*n, # 3D, vector
            dd.ndim == 2 and dd.shape[-1] == 3 and dd.size == 3*n, # 1D, vector
        ])



    def get_nodal_vals(self, data, dom=Domain.PM):
        """Return an array with nodal data values.

        Parameters
        ----------
        data : `str` or `numpy.ndarray`
            1) Datafile name (str) that will be read by
               `pyghs.parser.parse()`
            2) Data in numpy.ndarray indexed by node.
        dom : {`Domain.PM`, `Domain.FRAC`, "pm", "frac"}
            Domain on which to operate
        """
        dom = Domain.a2D(dom)

        ret = None

        d = data
        if type(data) == str:
            logger.debug(f'Nodal value read requested from file {data}')
            _count = self._get_count_from_filetype(data, dom)
            d = parse(data,count=_count)['data']

        if dom == Domain.PM:
            logger.debug(f'Nodal values being reshaped to {self.shape}')
            # hope we get a view instead of a copy
            if not self._is_PM_nodal(d):
                raise ValueError(
                    'data is incompatible wiht PM nodal data shape/size')
            ret = d.reshape(self.shape,order='F')

        elif dom == Domain.FRAC:
            logger.debug('Nodal values being reshaped to '
                    + f'{self.hgs_fx_nodes["link_frac2pm"].shape}')

            if self._is_PM_nodal(d):
                ret = d.flatten(order='F')[self.hgs_fx_nodes['link_frac2pm']]
            elif self._is_FRAC_nodal(d):
                ret = data
            else:
                raise ValueError(
                    'data is incompatible with fracture nodal data size/shape')

        else:
            raise NotImplementedError(f'Not implemented for domain {dom.name}')

        return ret

    def get_zoned_element_vals(self, zonedata, dom=Domain.PM):
        """Return an array with the zone data splayed out to each element

        e.g. Suppose a domain has two zones with different porosities,
        zonedata=[-1.0, 0.13, 0.10]. The result of this methos is a fully-filled
        array (of the same shape as the elements in the requested domain)
        with each entry being one of those three zonedata values.

        Parameters
        ----------
        zonedata : 1D array
            Array, where index is equal to the zone number, of the datum for
            each zone.

        .. Note::
            Zone indices typically start at index 1 (unless HGS' special
            "zone zero" is applied), so the normal use case is to pad
            index zero of the input zonedata array with a junk data value.
        """
        dom = Domain.a2D(dom)

        ret = None

        if dom==Domain.PM:
            ret = np.zeros(self.hgs_pm_elems['ne'])
            for i,izn in enumerate(self.hgs_pm_elems['zone']):
                ret[i] = zonedata[izn]
            ret = np.reshape(ret,self.elshape,order='F')
        else:
            raise NotImplementedError()

        return ret

    def _get_count_from_filetype(self, fn, dom):
        """From the file name, determine the count of expected data"""
        _count = None
        _dt = get_datatype(fn)
        if _dt is HGS_DATATYPE.NODAL:
            if dom == Domain.PM: _count = self.nn
            elif dom == Domain.FRAC: _count = self.nfn
        elif _dt is HGS_DATATYPE.ELEMENTAL:
            if dom == Domain.PM: _count = self.ne
            elif dom == Domain.FRAC: _count = self.nfe
        return _count

    def get_element_vals(self, data, dom=Domain.PM, method=None):
        """Return an array with nodal data values calculated per element.

        Parameters
        ----------
        data :
            1) Datafile name (str) that will be read by
               `pyghs.parser.parse()`
            2) Data in numpy.ndarray indexed by node.
        dom : `pyhgs.mesh.Domain`
            Domain on which to operate
        method : function, optional
            A function that accepts n values from the incident nodes of an
            element to compute the combined value.
            Default: n-point average, where n is the number of nodes per
            element.
        count : int (optional)
            The number of values expected in the Datafile (if a str is passed as
                    `data`. Otherwise, ignored.
        """
        dom = Domain.a2D(dom)

        if method:
            raise NotImplementedError('Cannot [yet] specify method')

        # data, format unspecified
        d = data
        if type(data) == str:
            logger.debug(f'Element value read requested from file {data}')
            _count = self._get_count_from_filetype(data, dom)
            d = parse(data,count=_count)['data']

        # function to calculate elemental data
        func = lambda i : i # identity method, for now

        # data, reshaped for sole argument to func, or error flag if None
        dd = None


        if dom == Domain.PM:

            def _grid_8pt_avg(a):
                """Do a non-weighted 8 point average of 3D array a"""

                if not self._is_PM_grid(a):
                    raise ValueError('Must be passed in 3D or 4D')
                if a.ndim == 3:
                  return ( a[:-1,:-1,:-1]
                    + a[ 1:,:-1,:-1]
                    + a[:-1, 1:,:-1]
                    + a[ 1:, 1:,:-1]
                    + a[:-1,:-1, 1:]
                    + a[ 1:,:-1, 1:]
                    + a[:-1, 1:, 1:]
                    + a[ 1:, 1:, 1:] ) / 8.
                elif a.ndim == 4:
                  return ( a[:-1,:-1,:-1,:]
                    + a[ 1:,:-1,:-1,:]
                    + a[:-1, 1:,:-1,:]
                    + a[ 1:, 1:,:-1,:]
                    + a[:-1,:-1, 1:,:]
                    + a[ 1:,:-1, 1:,:]
                    + a[:-1, 1:, 1:,:]
                    + a[ 1:, 1:, 1:,:] ) / 8.


            if self._is_PM_nodal_scalar(d):
                # 1D array of nodal values passed
                # Transform to 3D
                dd = d.reshape(self.shape,order='F')
                func = _grid_8pt_avg

            elif self._is_PM_elemental_scalar(d):
                # 1D array of elemental values passed
                # Transform to 3D
                # set method to Identity
                dd = d.reshape(self.elshape,order='F')

            elif self._is_PM_nodal_vector(d):
                # Assume that a flat list of triples have been passed
                # Transform this to 4-dimensional data and perform 'method' for
                # each index in the 4th dimension
                dd = d.reshape(self.shape+(3,), order='F')
                func = _grid_8pt_avg

            elif self._is_PM_elemental_vector(d):
                # Transform this to 4-dimensional data and perform 'method' for
                # each index in the 4th dimension
                dd = d.reshape((self.elshape+(3,)), order='F')

        elif dom == Domain.FRAC:

            with warnings.catch_warnings():
                #targeting the warning of the 0-based fracture node indices
                warnings.simplefilter('ignore')
                _f2p = self.hgs_fx_nodes['link_frac2pm']
                _finc = self.hgs_fx_elems['inc']

            # helper methods
            def _fx_node_avg_from_pm(a):
                """a is flattened PM nodal data"""
                ret = np.zeros(self.nfe)

                for i,fx_inc in enumerate(_finc):
                    ret[i] = np.sum(a[_f2p[fx_inc]])

                np.multiply(ret,1.0/self.hgs_fx_elems['nln'], out=ret)

                return ret

            def _fx_node_avg(a):
                """a is flattened FRAC nodal data, may have ndim = 1 or 2"""
            
                ret = None
                if a.ndim == 1:
                    ret = np.zeros(self.nfe)
                elif a.ndim == 2:
                    ret = np.zeros((self.nfe, a.shape[1]),)
                else:
                    raise NotImplementedError('Unexpected number of dimensions in "a"')

                for i,fx_inc in enumerate(_finc):
                    ret[i] = np.sum(a[fx_inc], axis=0)

                np.multiply(ret,1.0/self.hgs_fx_elems['nln'], out=ret)

                return ret

            if self._is_PM_nodal_scalar(d):
                # data seems to be PM nodal values
                dd = d.flatten(order='F')
                func = _fx_node_avg_from_pm

            elif self._is_FRAC_nodal_scalar(d) or self._is_FRAC_nodal_vector(d):
                # data seems to be FRAC nodal values
                dd = d
                func = _fx_node_avg

            elif self._is_FRAC_elemental_scalar(d):
                # data seems to be FRAC elemental values
                dd = d

            elif self._is_PM_nodal_vector(d):
                # seems to be flat array of tuples
                # data seems to be tuples of PM nodal values
                raise NotImplementedError()

            elif self._is_PM_elemental_vector(d):
                # data seems to be tuples of PM element values
                raise ValueError(
                    'Cannot interpret PM elemental data as FRAC data.')

            elif self._is_FRAC_elemental_vector(d):
                # data seems to be tuples of frac element values
                dd = d

            else:
                raise ValueError('Unhandled type/shape of data')
                    
        else:
            raise NotImplementedError(f'Not implemented for {dom}')

        if dd is None:
            raise ValueError('Unexpected size/shape found in input data, '
                f'{d.shape}, for domain {dom.name}.')

        return func(dd)

    def _is_rectilinear(self):
        """Return True if this is a rectilinear grid"""
        return 3 == sum(1 for gla in self.get_grid_lines() if gla is not None)

    def _pm_V(self):

        def _rect_V():
            (nx,ny,nz) = self.elshape
            gl = self.get_grid_lines()
            dx = gl[0][1:]-gl[0][:-1]
            dy = gl[1][1:]-gl[1][:-1]
            dz = gl[2][1:]-gl[2][:-1]

            return np.dot(
                    np.dot(dx.reshape((nx,1)),
                           dy.reshape((1,ny))).reshape((nx,ny,1)),
                    dz.reshape((1,nz))).reshape((nx,ny,nz),order='F')

        # is a regular grid
        # ... early exit
        if self._is_rectilinear():
            return _rect_V()
            
        #
        # decompose to 6 tetrahedra per brick
        #

        # some helper index slices
        #       iz
        #        ^
        #        | 7-----6
        #        |/|    /|
        #        4-+---5 |
        #        | 3-- +-2
        #        |/    |/
        #        0-----1------>ix
        #
        sl = [
            np.s_[:-1,:-1,:-1],
            np.s_[1:,:-1,:-1],
            np.s_[1:,1:,:-1],
            np.s_[:-1,1:,:-1],
            np.s_[:-1,:-1,1:],
            np.s_[1:,:-1,1:],
            np.s_[1:,1:,1:],
            np.s_[:-1,1:,1:], ]

        (nx,ny,nz) = self.shape
        N = self.hgs_fx_elems['nfe']
        p = self.hgs_pm_nodes['ncoords'].reshape(nx,ny,nz,3, order='F')

        ret = np.zeros(self.elshape)

        # for each of the six tetrahedra
        for (a,b,c,d) in [
                (0,1,3,5), (0,4,5,7), (0,3,5,7),
                (2,3,1,5), (2,3,5,7), (2,5,6,7),
            ]:

            # https://en.wikipedia.org/wiki/Tetrahedron#Volume
            ad = p[sl[d]] - p[sl[a]]
            abxac = np.cross(p[sl[b]]-p[sl[a]], p[sl[c]]-p[sl[a]])

            # V = 1/3 * B * h
            ret += 1./6 * norm(ad[:,:,:]*abxac[:,:,:], axis=3)

        return ret

    def _fx_V(self):
        """Compute fracture volumes.

        Computed as aperture times the area of each fracture element. The area
        is the sum of the area of two triangles in R3 (3D cartesian plane of
        real values).

        """
        
        N = self.hgs_fx_elems['nfe']
        c = self.hgs_pm_nodes['ncoords']
        fxap = self.hgs_fx_elems['ap'] #alias
        
        with warnings.catch_warnings():
            #targeting the warning of the 0-based fracture node indices
            warnings.simplefilter('ignore')
            _f2p = self.hgs_fx_nodes['link_frac2pm']

        # temporary arrays
        fxsz = np.zeros((N,3,3), dtype=np.float32)
        
        for i,(inc, izn, ap) in enumerate(
                self.iter_elements(Domain.FRAC)):

            # pm nodes indices must be used to get coordinates
            ipm = _f2p[inc]

            fxsz[i,0,:] = c[ipm[1]]-c[ipm[0]]
            fxsz[i,1,:] = c[ipm[2]]-c[ipm[0]]
            fxsz[i,2,:] = c[ipm[3]]-c[ipm[0]]

        # with help from...
        # https://math.stackexchange.com/questions/128991/
        #   how-to-calculate-the-area-of-a-3d-triangle
        # where this is aperture * 0.5 ( 2S_1 + 2S_2 ), where S_{1,2] are the
        # two triangle parts
        return fxap * 0.5 * \
            (_pa(fxsz[:,0,:], fxsz[:,1,:]) + _pa(fxsz[:,1,:], fxsz[:,2,:]))

    def intersect(self, iel, bb, dom=Domain.PM, with_proportion=False):
        """Return the volume of intersection of the element with bbox


        Parameters
        ==========
        iel : int
            The element index (0-based)
        bb : array
            The bounding box as an array of [x_min, x_max, .., z_min, z_max,]

        """
        dom = Domain.a2D(dom)

        if not self._is_rectilinear():
            raise NotImplementedError()


        def _pm_bb(iel):
            _nodes = self.hgs_pm_elems['inc'][iel]
            _bb = np.zeros(6)
            _bb[0::2] = self.hgs_pm_nodes['ncoords'][_nodes[0]]
            _bb[1::2] = self.hgs_pm_nodes['ncoords'][_nodes[6]]
            return _bb


        def _fx_bb(iel):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _nodes = self.hgs_fx_elems['inc'][iel]
                _n0 = self.hgs_fx_nodes['link_frac2pm'][_nodes[0]]
                _n2 = self.hgs_fx_nodes['link_frac2pm'][_nodes[2]]
            _bb = np.zeros(6)
            _bb[0::2] = self.hgs_pm_nodes['ncoords'][_n0]
            _bb[1::2] = self.hgs_pm_nodes['ncoords'][_n2]
            
            # determine orientation
            _pax=0
            if _n2-_n0 == self.shape[0]+1:
                _pax=2 # z perpendicular
            elif _n2-_n0 == (self.shape[0]*self.shape[1]+1):
                _pax=1
            # expand bounding plane to bounding box
            # NOTE due to floating point truncation, this may not perform well
            # in the calculation of intersection volume
            _hap = self.hgs_fx_elems['ap'][iel]/2.
            _bb[2*_pax] -= _hap
            _bb[2*_pax+1] += _hap

            return _bb


        elbb = None # element bounding box
        if dom == Domain.PM:
            elbb = _pm_bb(iel)
        elif dom == Domain.FRAC:
            elbb = _fx_bb(iel)
        else:
            raise NotImplementedError()

        ibb = np.zeros(6)
        np.maximum(bb[0::2],elbb[0::2],out=ibb[0::2])
        np.minimum(bb[1::2],elbb[1::2],out=ibb[1::2])

        ivol = np.prod(ibb[1::2]-ibb[0::2])
        if with_proportion:
            elvol = np.prod(elbb[1::2]-elbb[0::2])
            return ivol, ivol/elvol

        return ivol

    def get_element_volumes(self, dom=Domain.PM):
        """Return element volumes for the requested domain"""
        dom = Domain.a2D(dom)

        if dom == Domain.PM:
            return self._pm_V()
        elif dom == Domain.FRAC:
            return self._fx_V()
        else:
            raise ValueError(
                f'Cannot calculate volume for requested zone {dom}')

    def find_grid_index(self, *args):
        """Get a grid index (ix,iy,iz) given coordinate or node number

        Parameters
        ----------
        x, y, z : float
            The (x,y,z) coordinate to find, or
        inode : int
            The node index number (0-based).

        Returns
        -------
        numpy.array([ix,iy,iz],int32)
            the grid-line triple (0-based).
        """
        i = [-1,-1,-1,]
        gl = self.get_grid_lines()
        
        if len(args) == 3:
            x,y,z = map(float,args)
            for ii,gl,v in zip(count(),gl,(x,y,z,)):
                i[ii] = bisect_left(gl,v)
        elif len(args) == 1:
            inode = int(args[0])
            i[0] = inode % self.shape[0]
            i[1] = int((inode - i[0])/self.shape[0]) % self.shape[1]
            i[2] = int((inode - i[0] - self.shape[0]*i[1])
                        /(self.shape[0]*self.shape[1])
                       )
        else:
            raise ValueError('Unexpected argument types or count passed')

        return np.array(i,dtype='int32')

    def find_node_index(self, x,y,z, dom=Domain.PM):
        """Find the node index closest to coordinate (x,y,z)"""
        dom = Domain.a2D(dom)

        if dom != Domain.PM:
            raise NotImplementedError(f'Not implemented for {dom}')

        i = self.find_grid_index(x,y,z)
        return self.ng2ni(i)

    def ng2ni(self,i):
        """Node grid index (PM only, as list_like) to node sequence index"""
        return HGSGrid._to_index(i, self.shape)

    def elg2eli(self,i):
        """Element grid (PM only, as *list_like*) to element sequence index"""
        return HGSGrid._to_index(i, self.elshape)

    def eli2elg(self,i):
        """Element index (PM only, as `int`) to element grid index"""
        return HGSGrid._to_grid_index(i, self.elshape)

    def ni2eli(self, dom=Domain.PM):
        """Return map of node index to incident element indices"""
        if dom not in self._n2e:
            self._n2e[dom] = self._node2el(dom=dom)
        return self._n2e[dom]

    def make_pm_to_fx_adjacency(self):
        """Return a mapping of PM elements to adjacent fracture elements.

        Returns
        -------

        A `defaultdict`-like object
            Keys are PM element indices and
            values are `list`s of fracture element indices. The dictionary does not
            hold entries for PM elements that have no neighbouring fractures; when
            such an element is requested, an empty `list` is returned as the value.
            <br><br>
            PM element indices can be a 0-based integer index, or a `tuple` of
            0-based element grid indices.
        """

        if hasattr(self, '_pme2fxe'):
            return self._pme2fxe

        class _KeyCheckDict(defaultdict):
            def __init__(self, owner_hgsgrid=self):
                super().__init__(list)
                self._grid = owner_hgsgrid
                self._imax = owner_hgsgrid.hgs_pm_elems['ne']

            def __getitem__(self,i):
                if hasattr(i,'__len__'):
                    # assume grid index coordinates given
                    i = self._grid.elg2eli(i)
                if int(i) < 0 or i >= self._imax:
                    raise ValueError(f'{i} out of bounds of PM element indices')
                return super().__getitem__(i)

        r = _KeyCheckDict()

        # Note: the following is efficient in memory usage, but might be slow.
        # Do computation of first node of each fracture in batch
        # Do comutation of 2nd node or orientation in batch
        # Then, compute the rows of the adjacency matrix
        # TODO Vectorize computation of this function with numpy.ufunc.

        one = np.eye(3,dtype='int32')

        with warnings.catch_warnings():
            #targeting the warning of the 0-based fracture node indices
            warnings.simplefilter('ignore')
            _f2p = self.hgs_fx_nodes['link_frac2pm']

        for i,(fnodes, zone, ap) in enumerate(self.iter_elements(Domain.FRAC)):
            gi0 = self.find_grid_index(_f2p[fnodes[0]])
            gi2 = self.find_grid_index(_f2p[fnodes[2]])

            for a in range(3):
                if gi0[a] == gi2[a]:
                    if gi0[a] < self.elshape[a]:
                        r[ self.elg2eli(gi0) ].append(i)
                    if gi0[a]-1 >= 0:
                        r[ self.elg2eli(gi0-one[a]) ].append(i)
                    break

        self._pme2fxe = r
        return r

    def iter_supersample_distance_groups(self,
            maxd=None,
            groups=None,
            domains=(Domain.PM,)):
        """Generate supersample groups for the requested domains

        Parameters
        ----------
        maxd : `float` or list-like, optional
            The supersampling distance, or 3-tuple of distances to be applied to
            the x- y- and z- domains respectively
        groups : list-like, optional
            A 3-tuple, one entry for each dimension of the grid, of lists of
            (lo,hi PM-grid element index)-tuples of ready-made supersample
            groups. If groups are not specified, then they will be created
            according to the `maxd` argument.
        domains : list-like of `pyhgs.mesh.Domain`, optional
            List of domains to create supersample groups

        Returns
        -------
        A generator over tuples of 'distance groups'.


        Notes
        -----
        'Distance groups' are:
        - for `Domain.PM`, tuples of `((ixlo, ixhi), (iylo, iyhi), (izlo, izhi),)`,
            which are lo (inclusive) and hi (exclusive) indicies into the PM element
            grid in the respective directions, thus defining a block of space;
        - for `Domain.FRAC`, possibly empty lists of fracture elements within or
            on the edge of the preceding blocks, like `[]`, `[ifrac0,]`, or
            `[ifrac0, ifrac1, ... ]`.

        The order of these groups is what one would get from

            product( starmap(
                make_supersample_distance_groups,
                [ (dx, maxd[0]), (dy, maxd[1]), (dz, maxd[2]), ]
            ) )

        If multiple domains are specified, generated items are concatenated
        tuples of the preceding in the order specified in `domains`.

        e.g. `domains = [Domain.FRAC]` gives

        `( [ifx1, ifx2, ...], ),
        ...`

        e.g. `domains = [Domain.PM]` gives

        `((ixlo, ixhi), (iylo, iyhi), (izlo, izhi),),
        ...`
        
        e.g. `domains = [Domain.FRAC, Domain.PM]` gives

        `([ifx1, ...,], (ixlo, ixhi), (iylo, iyhi), (izlo, izhi),),
        ...`
        """

        # PM groups
        pm_ssgr = []

        if groups is not None:
            # assume groups are perfect as passed
            pm_ssgr = groups
        else:
            if not hasattr(maxd,'__len__'):
                maxd = tuple(map(float, 3*(maxd,)))

            gl = self.get_grid_lines()

            for a,maxda in zip(range(3),maxd):
                gla = gl[a]
                da = gla[1:]-gla[:-1]
                pm_ssgr.append(make_supersample_distance_groups(da,maxda))

        # make adjacency
        adj = self.make_pm_to_fx_adjacency()

        # dispatcher
        func = {
            Domain.PM: lambda grp: grp,
            Domain.FRAC: lambda grp: (self._find_fx_in_single_ssgrp(grp,adj),)
        }

        for pmgrp in product(*pm_ssgr):
            y = tuple()
            for d in domains:
                yy = func[d](pmgrp)
                y = y + yy

            yield y

    def _find_fx_in_single_ssgrp(self, ss_rng, pm2fxadj):
        _elg2eli = lambda t: self._to_index_nocheck(t,self.elshape)
        r = set()
        ranges = (range(*lohi) for lohi in ss_rng)
        for elg in product(*ranges):
            #r.update(pm2fxadj[self.elg2eli(elg)])
            r.update(pm2fxadj[_elg2eli(elg)])
        return list(sorted(r))

    def _yield_fx_in_ssgrp(self, ss_ranges, pm2fxadj):
        """Generates a list of incident fracture elements for each SS group

        The order of the yielded sets is the same as `product(*ss_ranges)`
        """
        # TODO
        # This might be improved by incrementally adding and deleting adjacent
        # fractrure elements from each group. e.g., for one dimension
        #
        # adjfx(lo1,hi1) := grp1
        # adjfx(lo2,hi2) = grp1 - adjfx(lo1,lo2) + adjfx(hi1,hi2) := grp2
        # adjfx(lo3,hi3) = grp2 - adjfx(lo2,lo3) + adjfx(hi2,hi3)
        # ...
        #
        # If there is a high degree of overlap between the ss regions, this
        # might be a significant improvement.
        for elindlohi in product(*ss_ranges):
            yield self._find_fx_in_single_ssgrp(elindlohi,pm2fxadj)


    def choose_nodes_block(self, blockspec, dom=Domain.PM):
        """Return a list of nodes contained in a 3D block"""

        _plq.push()

        dom = Domain.a2D(dom)

        if dom not in [Domain.PM, Domain.FRAC]:
            raise NotImplementedError()

        # zone 3D block bounds
        bb = _BoundingBox(*toDTuple(blockspec))

        gl = self.get_grid_lines()

        # grid line indices [inclusive, exclusive) for the bounding box
        ibb = -np.ones(6,dtype=int)

        icoord = iter(bb)
        for axis in range(3):
            if gl[axis] is not None:
                ibb[2*axis  ]= bisect_left(gl[axis],next(icoord))
                ibb[2*axis+1]= \
                    bisect_right(gl[axis],next(icoord),lo=ibb[2*axis])
            else:
                next(icoord)
                next(icoord)

        # helper methods
        def _choose_nodes_grid(_ibb):
          """Return a 3D grid of node indices"""
          # set of nodes to return
          _n = np.zeros(ibb[1::2]-ibb[0::2], dtype=int)
  
          # early exit!
          if _n.size == 0:
              return _n
  
          # calculate indices based on a regular grid
          _n[0,0,0] = HGSGrid._to_index_nocheck(ibb[0::2], self.shape)
          if _n.shape[0] > 1:
              _n[1:,0,0] = _n[0,0,0] + np.arange(1, _n.shape[0])
          if _n.shape[1] > 1:
              _n[:,1:,0] = _n[:,0,0][:,np.newaxis] \
                  + (self.shape[0]*np.arange(1, _n.shape[1]))[np.newaxis,:]
          if _n.shape[2] > 1:
              _n[:,:,1:] = _n[:,:,0][:,:,np.newaxis] \
                + ( 
                    np.prod(self.shape[0:2]) * np.arange(1,_n.shape[2])
                  )[np.newaxis,np.newaxis,:]
          return _n

        def _choose_nodes_nogrid(_ibb):

            # strategy generate a 3D array of *candidate* nodes in the block
            # assumes that this is still an i-j-k grid, just that the
            # coordinates are not in a regular grid

            # Then, test the coordinates of candidates (individually, or perhaps
            # by finding "bounding" values) to determine if they're in the
            # bounding box

            for i,gla in enumerate(gl):
                if gla is None:
                    _ibb[2*i] = 0     # candidate min coordinate
                    _ibb[2*i+1] = self.shape[i] # candidate max grid coordinate

            _n = _choose_nodes_grid(_ibb)


            _in_bb_funcs = []
            for a,gla in enumerate(gl):
                if gla is None:
                    _in_bb_funcs.append(
                        (a, lambda c: bb[2*a]<=c<=bb[2*a+1]),)

            def _is_in_bb(c):
                for i, predicate in _in_bb_funcs:
                    if not predicate(c[i]):
                        return False
                return True

            # TODO Not all nodes need to be checked --- Perhaps we can do a
            # binary search in the no-gridline axes where we can find a lower
            # and upper bound contained within the bounding box, then exclude
            # the others

            # OR
            # Re-build the gridlines for the "missing" axis. They _are not_
            # equally spaced domain-wide, but they _may_ be spaced evenly in one
            # plane, thys
            _all_coords = self.hgs_pm_nodes['ncoords']
            for ia,inode in enumerate(_n.flat):
                coord = _all_coords[inode]

                if not _is_in_bb(coord):
                    _n.flat[ia] = -1

            return _n


        # choose the PM nodes in the box
        _nn = None
        if any([gla is None for gla in gl]):
            _nn = _choose_nodes_nogrid(ibb).ravel(order='F')
            _nn = _nn[_nn > -1] # cull PM node indices
        else:
            _nn = _choose_nodes_grid(ibb).ravel(order='F')

        # reformat the nodes
        if dom == Domain.PM:
            _nn = list(_nn) # should already be sorted from ravel above

        elif dom == Domain.FRAC:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pm2fx = self.hgs_fx_nodes['link_pm2frac']
            _nn = pm2fx[_nn]      # do conversion for this node set
            _nn = _nn[ _nn > -1 ] # cull out nodes that don't map
            _nn = sorted(list(_nn))

        else:
            raise NotImplementedError(
                f'Cannot choose nodes for domain {dom.name}')

        _plq.pop(f'choose_nodes_block {dom.name} on {blockspec}')

        return _nn

    def _get_inc(self, dom, do_conversion=False):
        """Return the element->node incidence for the give domain"""
        if dom == Domain.PM:
            return self.hgs_pm_elems['inc']
        elif dom == Domain.FRAC:
            if do_conversion:
                _inc_fx = np.vectorize(self._ipm2ifrac)(self.hgs_fx_elems['inc'])
                return _inc_fx
            else:
                return self.hgs_fx_elems['inc']

    def choose_elements_block(self, blockspec, allow_partial=False,
            return_nodes=False,
            dom=Domain.PM):
        """Return a list of elements contained in a 3D block

        Parameters
        ----------
        blockspec : str
            Some specification of a 3D block: 1) As a string, "x1 x2 y1 y2 z1
            z2".
        allow_partial : bool
            If False, the element must be fully contained within the block. If
            True, interseciton of an element's corner, edge, face, or some
            partial volume with the block is sufficient for it to be included
            within the returned set
        dom : {`Domain.PM`, `Domain.FRAC`, 'pm', 'frac'}

        Returns
        -------
        `[element indices...]`, sorted in ascending order, or
        `([element indices...], [node indices...])`
            If `return_nodes` is `True`, then return a 2-tuple of the list of
            elements and the list of nodes, each in increasing order. The nodes
            returned here are within or on the edge of `blocspec` (i.e., _not_
            a the full set of nodes adjacent to the elements).
        """

        _plq.push()

        # guard against old code specifying a domain as a positional parameter
        # (where return_nodes is omitted) which gets interpreted as a bool for
        # return_nodes
        if type(return_nodes) != bool:
            raise ValueError(
                f'return_nodes must be a bool, not {type(return_nodes)}')

        dom = Domain.a2D(dom)

        _nodes = self.choose_nodes_block(blockspec, dom)
        _n2el = self._node2el(_nodes, dom)
        _els = set(chain.from_iterable(_n2el[_] for _ in _n2el))

        with warnings.catch_warnings():
            #targeting the warning of the 0-based fracture node indices
            warnings.simplefilter('ignore')
            _inc = self.get_elements_data(dom)['inc']

        def make_ret(n,e):
            if return_nodes:
                return (list(sorted(e)),list(sorted(n)))
            return list(sorted(e))

        _ret_els = None

        if allow_partial:
            _ret_els = _els

        else:
            # cull out elements who's node set is not fully in the block
            _ret_els = set()
            for iel in _els:
                _elinc = _inc[iel]
                if _elinc.size == np.intersect1d(_elinc,_nodes,True).size:
                    _ret_els.add(iel)

        _plq.pop(f'choose_elements_block {dom.name} {blockspec}')

        return make_ret(_nodes,_ret_els)


    def _node2el(self, nodes=None, dom=Domain.PM):
        """Return a mapping of nodes to incident elements

        Parameters
        ----------
        nodes : array-like
            Return a mapping for the given subset of nodes. Assumes that node
            indicies are given relative to the given domain (i.e., not the total
            node ordering or "3DNode#").

        Return
        ------
        A dict-like representation of the node->element mapping for the given
        domain and subset of node indices.
        """

        _plq.push('Making nodes->element incidence lists...')

        dom = Domain.a2D(dom)

        def _build_n2el_inc_pm():
            """work in progress"""
            ret = [set() for _ in range(self.nn)]

            inode = 0
            iel = 0
            # iterate in fortran-order
            for iez in range(self.elshape[2]):
                for iey in range(self.elshape[1]):
                    for iex in range(self.elshape[0]):
                        inode = self.ng2ni((iex  ,iey  ,iez  ,))
                        (nnx,nny,nnz) = self.shape
                        nnxy = nnx*nny
                        ret[inode].add(iel)
                        ret[inode+1].add(iel)
                        ret[inode+nnx].add(iel)
                        ret[inode+nnx+1].add(iel)
                        ret[inode+nnxy].add(iel)
                        ret[inode+nnxy+1].add(iel)
                        ret[inode+nnxy+nnx].add(iel)
                        ret[inode+nnxy+nnx+1].add(iel)
                        iel += 1

            return ret

        def _build_n2el_inc():

            # choose the nodes and elements sets
            if dom == Domain.PM:
                _nn = self.nn
            elif dom == Domain.FRAC:
                _nn = self.nfn
            else:
                raise NotImplementedError(f'Domain {dom} not implemented')

            with warnings.catch_warnings():
                #targeting the warning of the 0-based fracture node indices
                warnings.simplefilter('ignore')
                _el2ninc = self.get_elements_data(dom)['inc']

            # TODO reimplement with calculation by grid indices for PM
            ret = [set() for _ in range(_nn)]
            for iel,row in enumerate(_el2ninc):
                for inode in row:
                    ret[inode].add(iel)
            return ret


        # enable 'caching' of node->element incidence by storing them in the
        # hgs_fx_nodes and hgs_pm_nodes dictionaries in the key 'inc'
        if dom not in self._n2e:
            #if dom == Domain.PM:
            #    self._n2e[dom] = _build_n2el_inc_pm()
            #else:
            self._n2e[dom] = _build_n2el_inc()

        #alias
        _n2e = self._n2e[dom]

        # return object
        ret = None

        # If given, reformat the `nodes` data to a sorted array
        # Specify `ret` to hold all or just a subset of rows
        if type(nodes) is int:
            ret = {nodes:_n2e[nodes]}
        elif nodes is not None:
            ret = dict((i,_n2e[i]) for i in nodes)
        else:
            ret = _n2e

        _plq.pop('Making nodes->element incidence lists')

        return ret




def make_supersample_distance_groups(dx, maxd):
    """Return indices of unique, overlapping chunks of combined size <= maxd.

    This function will define sets of indices, "supersamples," for the purpose of
    creating a new spatial regions that are greater in size than each original
    cell by combining neigboring cells.

    Parameters
    ----------
    dx : list
        List of distances that represent "cell" widths, or increments
        between grid lines.
    maxd : float
        A value representing the maximum distance/size of the group along
        this dx-axis set.

    Returns
    -------
    [ (istart,iend], ... ]
        where no istart,end chunk is contained
        in any other. Each of the original increments is represented one or more
        times. In cases where an increment, idx, exceeds maxd, it will be
        represented in a zero-length chunk (idx,idx).
    """

    # deal with edge cases
    if len(dx) == 0:
        return []
    if len(dx) == 1:
        return [ [0, int(dx[0]<=maxd),], ]


    n = 0
    ssbl = np.empty((len(dx),2), dtype=np.int32) # super sample blocks
    ssbl[0][0] = -1
    NONE = ssbl[0][0]
    accd = 0

    # find initial entries that are > maxd
    for istart,d in enumerate(dx):
        if d > maxd:
            ssbl[n] = [istart,istart,]
            n += 1
        else:
            ssbl[n] = [istart,NONE,]
            n += 1
            accd = d
            break

    # build up supersample blocks
    istart = ssbl[n-1][0]+1

    for iend,d in enumerate(dx[istart:],start=istart):

        if accd <= 0 and d > maxd:
            ssbl[n-1][1] = iend
            ssbl[n] = [iend+1,NONE,]
            n += 1
            accd = 0
            continue

        accd = accd + d

        if accd > maxd:
            # end has been found; finish-up ssbl's last entry
            ssbl[n-1][1] = iend

            # begin a new entry
            istart = ssbl[n-1][0]
            ssbl[n] = [istart, NONE,]
            n += 1

            # increment the new entry's start value to accommodate this block
            while accd > maxd:
                accd = accd - dx[istart]
                istart += 1
                ssbl[n-1][0] = istart

            if istart > iend:
                ssbl[n-1] = [iend,iend,]
                if iend < len(dx)-1:
                    ssbl[n] = [iend+1,NONE,]
                    n += 1
                
    if ssbl[n-1][1] == NONE:
        ssbl[n-1][1] = len(dx)

    ret = list( list(r) for r in ssbl[:n] )
    return ret

def supersample( groups, d, *more_d, weights=None ):
    """Return supersampled version of **d**

    Supersampling is done according to the groups of 3D domain indicies provided
    in **groups**, where it assumed that the dimension and size of **d** is in
    accord with the number of groupings in **groups**.

    Parameters
    ----------
    groups : list_like
        Length of list must be equal to the number of axes in **d**; entries in
        the list must be list-likes of 3D element-grid start and end indices that
        constitute supersample groups.

    d : numpy.ndarray
        The data to be supersampled.

    weights : array_like, default None
        Either None (default) for equal weighting of all cells, or an
        array_like of the same dimensions of d with weighting values for
        each cell in d.

    """

    if more_d:
        raise NotImplementedError()

    if len(d.shape) != len(groups):
        raise ValueError(f"d's shape must match number of groupings lists")

    if weights is not None:
        if weights.shape != d.shape:
            raise ValueError('d and weights must have same shape')

    def _r1d(dd,gg,ww):
        """This problem, reduced to 1d
        Returns a ndarray( [ sum(data*weight), sum(weight) ], ... )
        of size (len(gg), 2)
        """
        _r = np.ndarray( (2, len(gg)) )

        wd = dd
        if ww is not None:
            np.multiply(ww, dd, out=wd)

        sumlo, sumhi = gg[0]
        wdsum = np.sum(wd[sumlo:sumhi])
        wsum = sumhi-sumlo
        if ww is not None:
            wsum = np.sum(ww[sumlo:sumhi])

        _r[:,0] = [wdsum, wsum]

        if ww is not None:
            for i,(lo,hi) in enumerate(gg[1:],start=1):
                wsum = wsum - np.sum(ww[sumlo:lo]) + np.sum(ww[sumhi:hi])
                wdsum = wdsum - np.sum(wd[sumlo:lo]) + np.sum(wd[sumhi:hi])
                _r[:,i] = [wdsum, wsum]
                sumlo = lo
                sumhi = hi

        else:
            for i,(lo,hi) in enumerate(gg[1:],start=1):
                wsum = hi-lo
                wdsum = wdsum - np.sum(wd[sumlo:lo]) + np.sum(wd[sumhi:hi])
                _r[:,i] = [wdsum, wsum]
                sumlo = lo
                sumhi = hi

        return _r

    # Note: case with zero-length group, (istart==iend)
    # Currently, nan will be returned in taht cell

    # set output array dimensions
    rshape = tuple(len(ax) for ax in groups)
    r = np.zeros(rshape)

    # TODO write a more generic calculation. What's the reduction?!

    # perform different calculation depending on number of input dimensions
    if len(groups) == 1:

# Hard implementation
#        wd = _r1d(d, groups[0], weights)
#        np.divide( wd[0], wd[1], out=r)

# Easy implementation:
        if weights is not None:
            for i,(lo,hi) in enumerate(groups[0]):
                w = weights[lo:hi]
                r[i] = np.average(d[lo:hi],weights=w)
        else:
            for i,(lo,hi) in enumerate(groups[0]):
                r[i] = np.average(d[lo:hi])

    elif len(groups) == 2:

# Easy implementation:
        if weights is not None:
          for i,(lo,hi) in enumerate(groups[0]):
            for j,(jlo,jhi) in enumerate(groups[1]):
              w = weights[lo:hi,jlo:jhi]
              r[i,j] = np.average(d[lo:hi,jlo:jhi], weights=w)
        else:
          for i,(lo,hi) in enumerate(groups[0]):
            for j,(jlo,jhi) in enumerate(groups[1]):
              r[i,j] = np.average(d[lo:hi,jlo:jhi])

    elif len(groups) == 3:

# Easy implementation:
        if weights is not None:
          for i,(lo,hi) in enumerate(groups[0]):
            for j,(jlo,jhi) in enumerate(groups[1]):
              for k,(klo,khi) in enumerate(groups[2]):
                w = weights[lo:hi,jlo:jhi,klo:khi]
                r[i,j,k] = np.average(d[lo:hi,jlo:jhi,klo:khi], weights=w)
        else:
          for i,(lo,hi) in enumerate(groups[0]):
            for j,(jlo,jhi) in enumerate(groups[1]):
              for k,(klo,khi) in enumerate(groups[2]):
                r[i,j,k] = np.average(d[lo:hi,jlo:jhi,klo:khi])
    else:
        raise NotImplementedError(f'Not implemented for # axes > 3.')

    return r

def _pa(a, b):
    """Return the parallelogram area of the R^3 vectors in a and b

    Parameters
    ----------
    a, b, : arrays of R^3 triples 
        i.e. shape is (N,3), where (x,y,z) is in axis 1

    Multiply this by 0.5 to get area of triangle.
    """
    return norm(np.cross(a,b), axis=1)


def determine_grid_lines(hgs_coordinates_pm):
    """Return the ([x-],[y-],[z-grid lines]) in this rectilinear grid.

    If the grid is not rectilinear, one or more sets of gridlines will be
    returned as `None`.
    """

    (nx,ny,nz) = tuple(hgs_coordinates_pm['n'+a] for a in 'xyz')
    coords = hgs_coordinates_pm['ncoords']

    # check whether the grid spacing consistent along the grid
    n_unique = [ np.unique(coords[:,a]).size for a in range(3) ]
    _gl = 3*[None,]

    # replace flag 'None' with grid line list
    if nx == n_unique[0]:
        _gl[0] = coords[0:nx,0]

    if ny == n_unique[1]:
        _gl[1] = coords[0:nx*ny:nx,1]

    if nz == n_unique[2]:
        _gl[2] = coords[0::nx*ny,2]


    return tuple(_gl)

class _PropFinderDict(dict):


    def __init__(self, m, d,):
        """Create a dict that "finds" properties from antoher dict"""
        self.m = m
        self.d = d
        self.od = dict()

        self._finders= {
            'porosity':_PropFinderDict._key_error_undef,
        }

    def __getitem__(self, k):
        if k in self.d:
            return self.d[k]
        elif k in self.od:
            return self.od[k]
        elif k in self._finders:
            self.od[k] = self._finders[k]()
            return self.od[k]
        else:
            raise KeyError()

    def __setitem__(self, k, v):
        if k in self.d:
            return self.d.__setitem__(k, v)
        else:
            return self.od.__setitem(k, v)

    def __contains__(self, k):
        return k in self.d or k in self.od
    def __len__(self):
        return len(self.d)+len(self.od)
    def __iter__(self):
        return chain(iter(self.d),iter(self.od))

    def __repr__(self):
        return f'<{self.__class__.__name__} with {len(self)} keys>'
    

    @staticmethod
    def _key_error_undef():
        raise KeyError(
            f'Key defined, but no implementation provided.')


class _PMPropFinderDict(_PropFinderDict):
    """
    Assumes self.d is pyhgs.parser.parse(prefix+'o.elements_pm')
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # define
        self._finders['porosity'] = partial(
                self._scrape_zone_data_from_mprops, self.m, 'porosity')

    @staticmethod
    def _get_mprops(pfx):

        mprops_dict = CaseInsensitiveDict()
        simdir = os.path.dirname(pfx)+os.path.sep
        grok = parse(pfx+'.grok')

        mprops_files = grok['files_mprops']
        for f in mprops_files:
            _f = f
            if not os.path.isabs(f):
                _f = os.path.normpath(simdir+f)
            d = parse(_f)
            mprops_dict.update(d)

        if len(mprops_files) == 0 and \
                os.path.isfile(simdir+'scratch_mprops'):
            mprops_dict = parse(simdir+'scratch_mprops')
            warnings.warn('Using scratch_mprops as the mprops datafile',
                    stacklevel=6)

        if len(mprops_dict) == 0:
            raise RuntimeError('Could not pick a .mprops file')

        return mprops_dict

    @staticmethod
    def _scrape_zone_data_from_mprops(mesh, k):
        """Return a mapping of zone to [scalar] property value"""
        mprops_dict = _PMPropFinderDict._get_mprops(mesh.prefix)
        eco_data = parse(f'{mesh.prefix}o.eco')
        zones = eco_data.get_pm_zone_properties()

        pvals = (1+len(zones))*[1.,]
        for izn,zname in zones:
            # izn is 1-based
            pvals[izn] = mprops_dict[zname][k]

        return mesh.get_zoned_element_vals(pvals,'pm')



class _FxPropFinderDict(_PropFinderDict):
    """
    Assumes self.d is pyhgs.parser.parse_elements_frac(self.prefix)
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._finders['porosity'] = partial(
                _FxPropFinderDict._ones, self.d['nfe'])

    @staticmethod
    def _ones(n):
        return np.ones(n)

