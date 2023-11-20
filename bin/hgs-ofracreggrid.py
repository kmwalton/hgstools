#!/usr/bin/env python
"""Convert an discrete (orthogonal) fracture grid to a regular grid.

Current implementation will take solute concentration (as identified by solutes
in grok file) and reproduce these values in zones for total aqueous
concentration and fracture pororosity-only concentration.

The equation for each gridcell, i, is
    C_i = \sum_j ( V_fraction * V * porosity * C )
        / \sum_j ( V_fraction * V * porosity )

Where j are all the original fracture and/or porous matrix elements that
intersect with volume i (V_fraction of j in i is non-zero).
"""

import sys
import os
import argparse
import re
import types
import warnings

import itertools
import numpy as np

from g360.scriptmeta import get_imprint

import tecplot as tp
import tecplot.constant as TPConst
from tecplot.constant import ValueLocation as TPValLoc


import pyhgs
from pyhgs.cli import PathToPrefix
from pyhgs.mesh import Domain, HGSGrid
from pyhgs.calcs import AvRegGrid


import logging
mylog = logging.getLogger(__name__)
calclog = logging.getLogger('pyhgs.calcs')
calcperflog = logging.getLogger('pyhgs.calcs_perf')

def merge_dicts(*dargs):
    ret = dict()
    if isinstance(dargs[0], types.GeneratorType):
        for d in dargs[0]:
            ret.update(d)
    else:
        for d in dargs:
            ret.update(d)
    return ret

def make_tp_dataset(pfx, solutes):
    """Create a new dataset with the new variables and metadata"""

    VARS=['x', 'y', 'z', 'pv_pm', 'pv_fx',]
    VARS.extend( s+suf
        for s,suf in itertools.product(solutes, ('_tot', '_fx',)))

    mylog.debug(f'Setting up dataset with variables {VARS}')

    tp.add_page() # add a page so the dataset does not conflict with previous
    ds = tp.active_frame().create_dataset( f'Regularized grid for {pfx}', VARS,)
    ds.aux_data['ORIGINAL_DATASET'] = pfx
    ds.aux_data['RUN_METADATA'] = get_imprint(os.path.abspath(__file__))

    return ds

def populate_tp_zone1_spatial(ds, rgrid):
    """Set up/populate zone spatial data"""

    VARLOC = 3*(TPValLoc.Nodal,)+(ds.num_variables-3)*(TPValLoc.CellCentered,)

    nodes = rgrid.make_grid_nodes()
    with tp.session.suspend():
        zn1 = ds.add_ordered_zone('Average, Regridded', rgrid.shape,
                locations=VARLOC, strand_id=1)
        zn1.values('x')[:] = nodes[0,:]
        zn1.values('y')[:] = nodes[1,:]
        zn1.values('z')[:] = nodes[2,:]

    return zn1

def populate_tp_zones(ds, ogrid, rgrid, solutes, pathto):

    zn1 = ds.zone(0)
    def zone_generator():
        yield 0, zn1
        for i in itertools.count(1):
            yield i, zn1.copy(
                    share_variables=[ds.variable(v) for v in range(5)])
    iter_zones =  zone_generator()
                

    # porosity & volume
    # Note: these are (MxN) arrays, M = new grid elements, N = original elements
    # _Must_ use the csr_matrix.multiply operation, else the return object (from
    # np.multiply()) is unusable.
    v_phi_pm = rgrid.e2g[Domain.PM].T.multiply( (ogrid.get_element_volumes('pm')
                *g.get_elements_data('pm')['porosity']).ravel(order='F'))
    v_phi_fx = rgrid.e2g[Domain.FRAC].T.multiply(
            ogrid.get_element_volumes('frac'))

    # M-sized arrays
    _V_pm = np.asarray(np.sum(v_phi_pm, axis=1))
    _V_fx = np.asarray(np.sum(v_phi_fx, axis=1))
    
    # list of datafiles
    datafiles = [ fn for fn in os.listdir(pathto) if
            re.match(r'.*\.conc_.*\.(\d+)$', fn) ]
    has_frac_df = any(re.match(r'.*\.conc_frac.*\.(\d+)$', fn)
            for fn in datafiles)
    _lensuff = len(re.match(r'.*\.(\d+)$', datafiles[0]).group(1))
    dfsuff = sorted(set(fn[-_lensuff:] for fn in datafiles))

    # Load datafiles
    for isuff, suff in enumerate(dfsuff):
        mylog.info(f'processing temporal suffix {suff}')

        izn, zn = next(iter_zones)

        # bass-ackward workaround?
        _bad_size = np.array(rgrid.shape)-[0,0,1]
        _padded = np.zeros(_bad_size)

        if izn == 0:
            #breakpoint()
            _padded[:-1,:-1,:] = _V_pm.reshape(rgrid.elshape, order='F')
            zn.values('pv_pm')[:] = _padded.ravel(order='F')
            _padded[:-1,:-1,:] = _V_fx.reshape(rgrid.elshape, order='F')
            zn.values('pv_fx')[:] = _padded.ravel(order='F')

        for s in grok['solute']:

            # read PM dataset
            d = pyhgs.parser.parse(f'{ppfx}o.conc_pm.{s}.{suff}')
            t = d['ts']

            # concentration, PM
            data_pm = ogrid.get_element_vals(d['data']).ravel(order='F')
            
            # concentration, FRAC
            if has_frac_df:
                d = pyhgs.parser.parse(f'{ppfx}o.conc_frac.{s}.{suff}')
            data_fx = ogrid.get_element_vals(d['data'], Domain.FRAC)

            # process re-gridded data

            # M-sized arrays
            _M_pm = np.sum(v_phi_pm.multiply(data_pm), axis=1)
            _M_fx = np.sum(v_phi_fx.multiply(data_fx), axis=1)

            # write data to zone
            with tp.session.suspend(), warnings.catch_warnings():
                warnings.simplefilter('ignore') # ignore divde by zero in _V_fx
                zn1.solution_time = float(t)

                _d_fx = np.nan_to_num(np.squeeze(np.asarray(_M_fx/_V_fx)),
                            copy=False)
                _d_tot = np.nan_to_num(np.squeeze(
                                np.asarray((_M_pm+_M_fx)/(_V_pm+_V_fx))),
                            copy=False)

                #breakpoint()

                # fails, incompatible size issue:
                #zn.values(f'{s}_fx')[:] = _d_fx
                #zn.values(f'{s}_tot')[:] = _d_pm

                # no failure; output data is bad
                #zn.values(f'{s}_fx')[:rgrid.elsize] = _d_fx
                #zn.values(f'{s}_tot')[:rgrid.elsize] = _d_pm


                # with workaround
                _padded[:-1,:-1,:] = _d_fx.reshape(rgrid.elshape, order='F')
                zn.values(f'{s}_fx')[:] = _padded.ravel(order='F')
                _padded[:-1,:-1,:] = _d_tot.reshape(rgrid.elshape, order='F')
                zn.values(f'{s}_tot')[:] = _padded.ravel(order='F')


            #if isuff == 0: return # DEBUG


if __name__ == '__main__':

    argp = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    argp.add_argument('-v', '--verbose',
        action='count',
        default=0,
        help='Set the verbosity.',
    )

    argp.add_argument('--out-file',
        default='out_regridded.dat',
        help='The name of the tecplot ascii-format datafile',
    )

    argp.add_argument('DISCRETIZATION',
        #metavar=('dx', 'dy', 'dz',),
        type=float,
        nargs=3,
        help='The target discretization, as "dx dy dz".',
    )

    argp.add_argument('PATH_TO_PREFIX',
            metavar='path/to/prefix',
            action=PathToPrefix,
            help='path/to/prefix of the HGS simulation',
        )
    args = argp.parse_args()

    # set up output logging
    if args.verbose > 0:
        for l in (mylog, calclog, calcperflog,):
            l.addHandler(logging.StreamHandler())
    if args.verbose > 1:
        mylog.setLevel(logging.INFO)
        calclog.setLevel(logging.INFO)
    if args.verbose > 2:
        mylog.setLevel(logging.DEBUG)
        calclog.setLevel(logging.DEBUG)
        calcperflog.setLevel(logging.INFO)

    mylog.info(get_imprint(os.path.abspath(__file__), True))

    pathto, pfx, junk = PathToPrefix.split(args.PATH_TO_PREFIX)
    ppfx = pathto + os.sep + pfx

    mylog.debug(f'Examining HGS problem "{pfx}"')

    # Load grid
    g = HGSGrid(ppfx)
    grok = pyhgs.parser.grok.parse(ppfx+'.grok')
    #mprops = merge_dicts(pyhgs.parser.mprops.parse(pathto+os.sep+fmprops)
    #    for fmprops in grok['files_mprops'])
    #breakpoint()

    # make averaged-grid
    ag = AvRegGrid(g, args.DISCRETIZATION)

    # set up destination tecplot datafile
    ds = make_tp_dataset(pfx, grok['solute'])
    populate_tp_zone1_spatial(ds, ag)
    populate_tp_zones(ds, g, ag, grok['solute'], pathto)


    # write tecplot datafile
    mylog.debug(f'Writing tecplot data to {args.out_file}')
    tp.data.save_tecplot_ascii(args.out_file, dataset=ds)

    sys.exit(0)
