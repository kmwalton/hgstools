'''Inspect Hydrogeosphere meshes; interperet binary files'''

from itertools import count,repeat
from collections import defaultdict
from bisect import bisect,bisect_left,bisect_right
from enum import IntEnum,auto
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from pyhgs.parser import (
    parse,
    parse_coordinates_pm,
    parse_elements_pm,
    parse_coordinates_frac,
    parse_elements_frac,
    )

__docformat__ = 'numpy'

class Domain(IntEnum):
    """Types of domains in Hydrogeosphere"""
    PM = auto()
    FRAC = auto()

class HGSGrid():
    '''Inspect a Hydrogeosphere rectilinear grid'''

    def __init__(self, prefix):

        self.prefix = prefix
        """Simulation path and prefix"""
        self.hgs_pm_nodes = None
        """Simulation PM node data from `pyhgs.parser.parse_coordinates_pm`"""
        self.hgs_pm_elems = None
        """Simulation PM node data from `pyhgs.parser.parse_elements_pm`"""
        self.hgs_fx_nodes = None
        """Simulation fracture node data from
        `pyhgs.parser.parse_coordinates_frac`
        """
        self.hgs_fx_elems = None
        """Simulation fracture node data from 
        `pyhgs.parser.parse_elements_frac`
        """

        def _my_call(*args,**kwargs):
            return args[0](*args[1:],**kwargs)
        def _my_call_2(a):
            return a[0](a[1:])

        if False:
            # failed implementation
            with Pool(processes=4) as pool:
                res = pool.map(
                    _my_call_2,
                    zip( [ parse_coordinates_pm,
                        parse_elements_pm,
                        parse_coordinates_frac,
                        parse_elements_frac,
                       ],
                       4*[prefix,]
                    ) )

                self.hgs_pm_nodes = res[0]
                self.hgs_pm_elems = res[1]
                self.hgs_fx_nodes = res[2]
                self.hgs_fx_elems = res[3]
        elif False:
            # failed implementation
            with ThreadPoolExecutor(max_workers=4) as e:
                r0 = e.submit(parse_coordinates_pm,prefix)
                r1 = e.submit(parse_elements_pm,prefix)
                r2 = e.submit(parse_coordinates_frac,prefix)
                r3 = e.submit(parse_elements_frac,prefix)

                self.hgs_pm_nodes = r0.result()
                self.hgs_pm_elems = r1.result()
                self.hgs_fx_nodes = r2.result()
                self.hgs_fx_elems = r3.result()

        else:
            self.hgs_pm_nodes = parse_coordinates_pm(prefix)
            self.hgs_pm_elems = parse_elements_pm(prefix)
            self.hgs_fx_nodes = parse_coordinates_frac(prefix)
            self.hgs_fx_elems = parse_elements_frac(prefix)

        if self.hgs_pm_nodes['tetramesh'] or self.hgs_pm_nodes['nb2d']>0:
            raise ValueError(
                f'Simulation {prefix} does not have a rectilinear grid')

        self.nn = self.hgs_pm_nodes['nn'] # aliases
        """Number of PM nodes"""
        self.ne = self.hgs_pm_elems['ne']
        """Number of PM elements"""
        self.nfn = self.hgs_fx_nodes['nnfrac']
        """Number of fracture nodes"""
        self.nfe = self.hgs_fx_elems['nfe']
        """Number of fracture elements"""

        self.shape = tuple(self.hgs_pm_nodes[a] for a in ['nx','ny','nz'])
        """Shape of the PM node grid"""
        self.elshape = tuple(self.hgs_pm_nodes[a]-1 for a in ['nx','ny','nz'])
        """Shape of the PM element grid"""

    @staticmethod
    def _to_index_nocheck(t,shp):
        return t[0] + shp[0]*t[1] + shp[0]*shp[1]*t[2]

    @staticmethod
    def _to_index(t,shp):
        """Return 3D grid index as linear index"""
        if any(t[i] >= shp[i] for i in range (3)):
            s = tuple( v-1 for v in shp )
            raise ValueError(
                f'Element index {t} out of 3D grid index bounds. Need <={s}.')

        if any(v<0 for v in t):
            raise ValueError(
                'Element index {t} out of 3D grid index bounds. '\
                f'Need >={3*(0,)}.')

        return HGSGrid._to_index_nocheck(t,shp)

    def get_grid_lines(self):
        """Return the ([x-],[y-],[z-grid lines]) in this rectilinear grid"""
        # store
        if not hasattr(self,'gl'):
            (nx,ny,nz) = self.shape
            self.gl = (
                self.hgs_pm_nodes['ncoords'][0:nx,0],
                self.hgs_pm_nodes['ncoords'][0:nx*ny:nx,1],
                self.hgs_pm_nodes['ncoords'][0::nx*ny,2],
            )
        return self.gl


    def get_n_elements(self,dom=Domain.PM):
        """Return the number of elements in the specified domain"""

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
            <code>( (node indices), zone, *<other data\>* )</code>

            where <code>*<otherdata\>*</code>is:
                for <code>dom == Domain.PM</code>, nothing;
                for <code>dom == Domain.FRAC</code>, ap.
        """

        if dom == Domain.PM:
            raise NotImplementedError()

        elif dom == Domain.FRAC:
            for inc,ap,zn in zip(
                self.hgs_fx_elems['inc'],
                self.hgs_fx_elems['ap'],
                self.hgs_fx_elems['zone'] ):

                yield (inc, zn, ap,)

        else:
            raise NotImplementedError(f'Not implemented for {dom}')

    def get_coords(self, iel, dom=Domain.PM):
        """
            Arguments:

                iel : int or tuple
                    Element index, or
                    Element from grid indices (ix, iy, iz), for appropriate
                    domain types.

            Returns:
                Tuple of 3D coordinates ( c0, c1, ... )
        """

        if dom == Domain.PM:
            if type(iel) == tuple:
                iel = self._to_index_nocheck(iel,self.elshape)
            inc_nodes = self.hgs_pm_elems['inc'][iel]
            return tuple(self.hgs_pm_nodes['ncoords'][i] for i in inc_nodes)

        elif dom == Domain.FRAC:
            inc_nodes = self.hgs_fx_elems['inc'][iel]
            return tuple(self.hgs_pm_nodes['ncoords'][i] for i in inc_nodes)

        else:
            raise NotImplementedError(f'Not implemented for {dom}')
        

    def get_nodal_vals(self, data, dom=Domain.PM):
        """Return an array with nodal data values.

        Arguments:
            data : str or `numpy.ndarray`
                1) Datafile name (str) that will be read by
                   `pyghs.parser.parse()`
                2) Data in numpy.ndarray indexed by node.

            dom : pyhgs.mesh.Domain
                Domain on which to operate

            method : function (optional)
                A function that accepts n values from the incident nodes of an
                element to compute the combined value.
                Default: n-point average, where n is the number of nodes per
                element.
        """
        ret = None

        d = data
        if type(data) == str:
            d = parse(data)['data']

        if dom == Domain.PM:
            # hope we get a view instead of a copy
            ret = d.reshape(self.shape,order='F')

        elif dom == Domain.FRAC:
            ret = d.flatten(order='F')[self.hgs_fx_nodes['link_frac2pm']]

        else:
            raise NotImplementedError(f'Not implemented for {dom!s}')

        return ret

    def get_zoned_element_vals(self, zonedata, dom=Domain.PM):
        """Return an array with the zone data splayed out to each element

        e.g. Suppose a domain has two zones with different porosities,
        zonedata=[-1.0, 0.13, 0.10]. The result of this methos is a fully-filled
        array (of the same shape as the elements in the requested domain)
        with each entry being one of those three zonedata values.

        Arguments:
            zonedata : 1D array
                Array, where index is equal to the zone number, of the datum for
                each zone.

            .. Note::
                Zone indices typically start at index 1 (unless HGS' special
                "zone zero" is applied), so the normal use case is to pad
                index zero of the input zonedata array with a junk data value.
        """
        ret = None

        if dom==Domain.PM:
            ret = np.zeros(self.hgs_pm_elems['ne'])
            for i,izn in enumerate(self.hgs_pm_elems['zone']):
                ret[i] = zonedata[izn]
            ret = np.reshape(ret,self.elshape,order='F')
        else:
            raise NotImplementedError()

        return ret

    def get_element_vals(self, data, dom=Domain.PM, method=None):
        """Return an array with nodal data values calculated per element.

        Arguments:
            data :
                1) Datafile name (str) that will be read by
                   `pyghs.parser.parse()`
                2) Data in numpy.ndarray indexed by node.

            dom : pyhgs.mesh.Domain
                Domain on which to operate

            method : function (optional)
                A function that accepts n values from the incident nodes of an
                element to compute the combined value.
                Default: n-point average, where n is the number of nodes per
                element.
        """

        if method:
            raise NotImplementedError('Cannot [yet] specify method')

        # data, format unspecified
        d = data
        if type(data) == str:
            d = parse(data)['data']

        # function to calculate elemental data
        func = lambda i : i # identity method, for now

        # data, reshaped for sole argument to func, or error flag if None
        dd = None

        if dom == Domain.PM:

            def _grid_8pt_avg(a):
                """Do a non-weighted 8 point average of 3D array a"""
                return ( a[:-1,:-1,:-1]
                    + a[ 1:,:-1,:-1]
                    + a[:-1, 1:,:-1]
                    + a[ 1:, 1:,:-1]
                    + a[:-1,:-1, 1:]
                    + a[ 1:,:-1, 1:]
                    + a[:-1, 1:, 1:]
                    + a[ 1:, 1:, 1:] ) / 8.


            if len(d.shape) == 1:
                # flat ordering

                if d.size == self.nn:
                    # 1D array of nodal values passed
                    # Transform to 3D
                    dd = d.reshape(self.shape,order='F')
                    func = _grid_8pt_avg

                elif d.size == self.ne:
                    # 1D array of elemental values passed
                    # Transform to 3D
                    # set method to Identity
                    dd = d.reshape(self.elshape,order='F')

            elif len(d.shape) == 2:
                # flat ordering of tuples

                if d.size == d.shape[1]*self.nn:
                    # Assume that a flat list of triples have been passed
                    # Transform this to 4-dimensional data and perform 'method' for
                    # each index in the 4th dimension
                    dd = np.ndarray(self.shape+(3,))
                    for idim in range(d.shape[1]):
                        dd[:,:,:,idim] = d[:,idim].reshape(self.shape, order='F')
                    func = _grid_8pt_avg

                elif d.size == d.shape[1]*self.ne:
                    # Transform this to 4-dimensional data and perform 'method' for
                    # each index in the 4th dimension
                    dd = np.ndarray(self.elshape+(3,),dtype=np.float32)
                    for idim in range(d.shape[1]):
                        dd[:,:,:,idim] = d[:,idim].reshape(self.elshape, order='F')

            elif len(d.shape) == 3:
                # 3D grid ordering

                if d.shape == self.shape:
                    # 3D array of nodal values passed in
                    # Assume this is already in the correct shape
                    dd = d
                    func = _grid_8pt_avg

                elif d.shape == self.elshape:
                    # Assume already in elemental format
                    # set method to Identity
                    dd = d

        elif dom == Domain.FRAC:

            # helper methods
            def _fx_node_avg(a):
                """a is flattened PM nodal data"""
                ret = np.zeros(self.nfe)

                for i,fx_inc in enumerate(self.hgs_fx_elems['inc']):
                    ret[i] = np.sum(a[fx_inc])

                np.multiply(ret,1.0/self.hgs_fx_elems['nln'], out=ret)

                return ret


            if len(d.shape) == 1:
                # flat array
                if d.size == self.nn:
                    # data seems to be PM nodal values
                    dd = d
                    func = _fx_node_avg

            elif len(d.shape) == 2:
                # seems to be flat array of tuples
                if d.size == d.shape[1]*self.nn:
                    # data seems to be tuples of PM nodal values
                    raise NotImplementedError()
                elif d.size == d.shape[1]*self.nfn:
                    # data seems to be tuples of frac nodal values
                    raise NotImplementedError()
                if d.size == d.shape[1]*self.ne:
                    # data seems to be tuples of PM element values
                    raise NotImplementedError()
                elif d.size == d.shape[1]*self.nfe:
                    # data seems to be tuples of frac element values
                    dd = d
                    
            elif len(d.shape) == 3:
                # seems to be grid array of scalars
                if d.shape == self.shape:
                    # data seems to be PM nodal values, as grid
                    dd = np.flatten(d, order='F')
                    func = _fx_node_avg

            elif len(d.shape) == 4:
                # seems to be grid array of tuples
                raise NotImplementedError()

            # TODO speed this up?
            # Determine the orientation of each element so that array
            # manipulations like above can be used(??)


        else:
            raise NotImplementedError(f'Not implemented for {dom}')

        if dd is None:
            raise ValueError('Unexpected size/shape found in input data, '
                f'{d.shape}, for {dom!s}.')

        return func(dd)


    def find_grid_index(self, *args):
        """Get a grid index (ix,iy,iz) given coordinate or node number

        Arguments:
            x, y, z : float
                The (x,y,z) coordinate to find, or
            inode : int
                The node index number (0-based).

        Returns:
            numpy.array([ix,iy,iz],int32), the grid-line triple (0-based).
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

        if dom != Domain.PM:
            raise NotImplementedError(f'Not implemented for {dom}')

        i = self.find_grid_index(x,y,z)
        return self.ng2ni(i)

    def ng2ni(self,i):
        """Node grid index (PM only, as list_like) to node sequence index"""
        return HGSGrid._to_index(i, self.shape)

    def elg2eli(self,i):
        """Element grid (PM only, as list_like) to element sequence index"""
        return HGSGrid._to_index(i, self.elshape)

    def make_pm_to_fx_adjacency(self):
        """Return a dict of pm element indices to adjacent fracture elements"""

        class _KeyCheckDict(defaultdict):
            def __init__(self, elszpm=self.hgs_pm_elems['ne']):
                super().__init__(list)
                self._imax = elszpm

            def __getitem__(self,i):
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

        for i,(fnodes, zone, ap) in enumerate(self.iter_elements(Domain.FRAC)):
            gi0 = self.find_grid_index(fnodes[0])
            gi2 = self.find_grid_index(fnodes[2])

            for a in range(3):
                if gi0[a] == gi2[a]:
                    if gi0[a] < self.elshape[a]:
                        r[ self.elg2eli(gi0) ].append(i)
                    if gi0[a]-1 >= 0:
                        r[ self.elg2eli(gi0-one[a]) ].append(i)
                    break

        return r

def make_supersample_distance_groups(dx, maxd):
    """Return indices of unique, overlapping chunks of combined size <= maxd.

    This function will define sets of indices, "supersamples," for the purpose of
    creating a new spatial regions that are greater in size than each original
    cell by combining neigboring cells.

    Arguments:
        dx : list
            List of distances that represent "cell" widths, or increments
            between grid lines.

    Returns:
        The list [ (istart,iend], ... ] where no istart,end chunk is contained
        in any other. Each of the original increments is represented one or more
        times. In cases where an increment, idx, exceeds maxd, it will be
        represented in a zero-length chunk (idx,idx).
    """

    # deal with edge cases
    if len(dx) == 0:
        return []
    if len(dx) == 1:
        return [ [0, int(dx[0]<=maxd),], ]


    ssbl = [] # super sample blocks
    accd = 0

    # find initial entries that are > maxd
    for istart,d in enumerate(dx):
        if d > maxd:
            ssbl.append([istart,istart,])
        else:
            ssbl.append([istart,None,])
            accd = d
            break

    # build up supersample blocks
    istart = ssbl[-1][0]+1

    for iend,d in enumerate(dx[istart:],start=istart):

        if accd <= 0 and d > maxd:
            ssbl[-1][1] = iend
            ssbl.append([iend+1,None,])
            accd = 0
            continue

        accd = accd + d

        if accd > maxd:
            # end has been found; finish-up ssbl's last entry
            ssbl[-1][1] = iend

            # begin a new entry
            istart = ssbl[-1][0]
            ssbl.append([istart, None,])

            # increment the new entry's start value to accommodate this block
            while accd > maxd:
                accd = accd - dx[istart]
                istart += 1
                ssbl[-1][0] = istart

            if istart > iend:
                ssbl[-1] = [iend,iend,]
                if iend < len(dx)-1:
                    ssbl.append([iend+1,None,])
                
    if ssbl[-1][1] is None:
        ssbl[-1][1] = len(dx)

    return ssbl

def supersample( groups, d, *more_d, weights=None ):
    """Return supersampled version of **d**

    Supersampling is done according to the groups of 3D domain indicies provided
    in **groups**, where it assumed that the dimension and size of **d** is in
    accord with the number of groupings in **groups**.

    Arguments
    ---------
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
