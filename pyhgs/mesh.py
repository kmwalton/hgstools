'''Inspect Hydrogeosphere meshes; interperet binary files'''

from itertools import count,repeat
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

class Domain(IntEnum):
    """Types of domains in Hydrogeosphere"""
    PM = auto()
    FRAC = auto()

class HGSGrid():
    '''Inspect a Hydrogeosphere rectilinear grid'''

    def __init__(self, prefix):

        self.prefix = prefix
        self.hgs_pm_nodes = None
        self.hgs_pm_elems = None
        self.hgs_fx_nodes = None
        self.hgs_fx_elems = None

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

        self.nn = self.hgs_pm_nodes['nn'] # alias

        self.shape = tuple(self.hgs_pm_nodes[a] for a in ['nx','ny','nz'])
        self.elshape = tuple(self.hgs_pm_nodes[a]-1 for a in ['nx','ny','nz'])

    def _to_index(self, t):
        if any(t[i] >= self.shape[i]-1 for i in range (3)):
            s = tuple( v-2 for v in self.shape )
            raise ValueError(
                f'Element index {t} out of element grid bounds, <={s}.')

        if any(v<0 for v in t):
            raise ValueError(
                f'Element index {t} out of element grid bounds, >={3*(0,)}.')

        return t[0] + self.shape[0]*t[1] + self.shape[0]*self.shape[1]*t[2]

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
            ( (node indices), zone, <other data> )

            where <otherdata> is,
            for dom == Domain.PM:
                nothing
            for dom == Domain.FRAC:
                ap
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
                iel = self._to_index(iel)
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
        ret = None

        d = data
        if type(data) == str:
            d = parse(data)['data']

        if dom == Domain.PM:
            # hope we get a view instead of a copy
            ret = d.reshape(self.shape,order='F')

        else:
            raise NotImplementedError(f'Not implemented for {dom}')

        return ret

    def get_zoned_element_vals(self, zonedata, dom=Domain.PM):
        """Return an array with the zone data applied to each element

        Arguments:
            zonedata : 1D array
                Array, where index is equal to the zone number, of the datum for
                each zone.

                Note: zones typically start at index 1 (unless HGS' special zone
                zero is applied), so the normal use case is to pad index zero of
                the array with a junk data value.
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

        ret = None

        if method:
            raise NotImplementedError('Cannot [yet] use method')

        d = data
        if type(data) == str:
            d = parse(data)['data']

        if dom == Domain.PM:
            # hope we get a view instead of a copy
            dd = d.reshape(self.shape,order='F')
            # calculate an 8-point average
            ret = ( dd[:-1,:-1,:-1]
                    + dd[ 1:,:-1,:-1]
                    + dd[:-1, 1:,:-1]
                    + dd[ 1:, 1:,:-1]
                    + dd[:-1,:-1, 1:]
                    + dd[ 1:,:-1, 1:]
                    + dd[:-1, 1:, 1:]
                    + dd[ 1:, 1:, 1:] ) / 8.

        elif dom == Domain.FRAC:
            ret = np.zeros((self.hgs_fx_elems['nfe'],))

            dd = d.flatten(order='F')

            # check that data is PM, nodal values ... assuming this check is
            # sufficient
            if dd.shape != self.nn:
                raise ValueError('Expected data size to be {self.nn}')

            inc = self.hgs_fx_elems['inc']

            # TODO speed this up?
            # Determine the orientation of each element so that array
            # manipulations like above can be used(??)

            # do 'i' loop first because fracture nodes will likely be neighbours
            # in the nodal data
            for i,fx_inc in enumerate(inc):
                for j in fx_inc:
                    ret[i] += dd[j]

            np.multiply(ret,1.0/self.hgs_fx_elems['nln'], out=ret)

        else:
            raise NotImplementedError(f'Not implemented for {dom}')

        return ret



    def find_grid_index(self, x,y,z):
        """Returns the grid index (ix,iy,iz) closest to coordinate (x,y,z)"""
        i = [-1,-1,-1,]
        for ii,gl,v in zip(count(),self.gl,(x,y,z,)):
            i[ii] = bisect_left(gl,v)
        return tuple(i)

    def find_node_index(self, x,y,z, dom=Domain.PM):
        """Find the node index closest to coordinate (x,y,z)"""

        if dom != Domain.PM:
            raise NotImplementedError(f'Not implemented for {dom}')

        (ix,iy,iz) = self.find_grid_index(x,y,z)
        return ix + self.shape[0]*(iy + self.shape[1]*iz)

def supersample_distance(dx, maxd):
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
    ssbl = [] # super sample blocks
    accd = 0

    # find initial entries that are > maxd
    for istart,d in enumerate(dx[:-1]):
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
                ssbl.append([iend+1,None,])
                

    ssbl[-1][1] = len(dx)

    return ssbl
