#!/usr/bin/env python
"""Convert an discrete (orthogonal) fracture grid to a regular grid.

Current implementation will take solute concentration (as identified by solutes
in grok file) and reproduce these values in zones for total aqueous
concentration and fracture pororosity-only concentration.
"""

import sys
import os
import argparse

import itertools
import numpy as np

from g360.scriptmeta import get_imprint

import tecplot as tp
import tecplot.constant as TPConst


import pyhgs
from pyhgs.cli import PathToPrefix
from pyhgs.mesh import HGSGrid
from pyhgs.calcs import AvRegGrid


import logging
mylog = logging.getLogger(__name__)
calclog = logging.getLogger('pyhgs.calcs')
calcperflog = logging.getLogger('pyhgs.calcs_perf')

def make_tp_dataset(pfx, solutes):
    """Create a new dataset with the new variables and metadata"""

    VARS=['x', 'y', 'z',]
    # TODO add volume of each porosity?
    VARS.extend( s+suf
        for s,suf in itertools.product(solutes, ('_tot', '_fx',)))

    mylog.debug(f'Setting up dataset with variables {VARS}')

    tp.add_page() # add a page so the dataset does not conflict with previous
    ds = tp.active_frame().create_dataset( f'Regularized grid for {pfx}', VARS,)
    ds.aux_data['ORIGINAL_DATASET'] = pfx
    ds.aux_data['RUN_METADATA'] = get_imprint(os.path.abspath(__file__))

    return ds

def make_tp_zones(ds, ogrid, solutes, rgrid):
    """Set up/populate zone spatial data"""

    nodes = rgrid.make_grid_nodes()
    
    zn1 = ds.add_ordered_zone('Average, Regridded', rgrid.shape)
    zn1.values('x')[:] = nodes[0,:]
    zn1.values('y')[:] = nodes[1,:]
    zn1.values('z')[:] = nodes[2,:]



if __name__ == '__main__':

    argp = argparse.ArgumentParser()
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
        metavar=('dx', 'dy', 'dz',),
        type=float,
        nargs=3,
        help='The target discretization.',
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

    # make averaged-grid
    ag = AvRegGrid(g, args.DISCRETIZATION)

    # set up destination tecplot datafile
    ds = make_tp_dataset(pfx, grok['solute'])
    zones = make_tp_zones(ds, g, grok['solute'], ag)

    # TODO
    # Load datafiles
    # concentration, porosity
    # process re-gridded data

    # write tecplot datafile
    mylog.debug(f'Writing tecplot data to {args.out_file}')
    tp.data.save_tecplot_ascii(args.out_file, dataset=ds)

    sys.exit(0)
