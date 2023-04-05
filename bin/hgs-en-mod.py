#!/usr/bin/env python
"""Perform modifications to an ._en-suffixed file

Hydrogeosphere produces binary datafiles for solutions to head (.hen),
concentration (.cen), etc. This program performs some basic manipulations of
these files.

For rectilinear grids:
    - truncate the input ._en to a smaller grid size (given a new, smaller array
            size)
    - interpolate the input ._en to a more refined grid (given a
            directory/prefix of an identical, but more spatially refined HGS
            domain)
"""



import sys
import os
import argparse
import re

import numpy as np
from scipy.io import FortranFile, FortranEOFError
from scipy.interpolate import RegularGridInterpolator

import matplotlib.pyplot as plt

from pyhgs.cli import PathToPrefix
from pyhgs.mesh import HGSGrid



argp = argparse.ArgumentParser()
argp.add_argument('FILE_IN',
    metavar='path/to/prefixo.Xen',
    action=PathToPrefix,
    help='''Input file path and name, where "X" is h for head, c for
    concentration, etc. HGS input simulation prefix is inferred from this
    value.''',
    )
slicegr = argp.add_mutually_exclusive_group()
slicegr.add_argument('--islice',
    type=str,
    metavar='"[i-range,j-range,k-range]"',
    help='''Retain the data specified by a numpy-style slice using array
    indices, inclusive ov square brackets, e.g., "[10:,0:-5,:]" to remove
    gridlinex 0 through 9, in x-, the last 5 gridlines in y-, and retain all
    gridlines in z-. This option requires knowledge of the spatial location of
    the gridlines you're slicing off and beware that this may cause a
    translation when it is read by the new HGS domain. Note: indices provided
    here are zero-based.''',
    )
slicegr.add_argument('--slice',
    type=str,
    metavar='"[x-range,y-range,z-range]"',
    help='''Retain the data specified by a numpy-style slice using model
    coordinates, inclusive ov square brackets, e.g., "[9:,0:-5,:]" to remove
    model distance 0 through 8.99999 units , in x-, the last 4.99999 units in
    y-, and retain all of z-.''',
    )

argp.add_argument('--interp-to',
    metavar='path/to/newdomainprefix',
    help='''Interpolate values in the input grid to the grid found in the
    simulation path/to/newdomainprefixo.coordinates_pm''',
    )

argp.add_argument('--show',
    action='store_true',
    default=False,
    help='''Show a rendering of z-layer zero of the two files''',
    )

argp.add_argument('FILE_OUT',
    metavar='path/to/outfile.Xen',
    action=PathToPrefix,
    help='''Output file path and name.'''
    )

args = argp.parse_args()


def _parse_en(fn, dtype, shape=None):
    with FortranFile(fn,'r') as fin:
        #ts = fin.read_ints(dtype=np.byte)
        d = fin.read_reals(dtype=dtype)
    if shape:
        d = d.reshape(shape,order='F')
    return d


path, pfx, fnrest = PathToPrefix.split(args.FILE_IN)
pathpfx = os.path.join(path,pfx)

grid = HGSGrid(pathpfx)

din = _parse_en(args.FILE_IN, np.double, shape=grid.shape)
dout = None
glout = 3*[None,]

if args.islice:
    try:
        dout = eval('din'+args.islice)
    except BaseException as e:
        argp.error('Error when processing --islice:\n'+str(e))

    glout = grid.get_grid_lines()
    for iax,axsl in enumerate(args.islice[1:-1].split(',')):
        if glout[iax] is None:
            continue
        glout[iax] = eval(f'glout[iax][{axsl}]')

elif args.slice:
    d = r'[+-]?\d+(?:\.\d*)?' # digit
    capr = '('+d+')?:('+d+')?' # captured range
    m = re.match(f'\[{capr},{capr},{capr}\]', args.slice)

    if not m:
        argp.error(f'Error when processing --slice {args.slice}')

    _i = 6*['',]

    gl = grid.get_grid_lines()

    for i,v in enumerate(m.groups()):

        if v is None:
            continue

        iax = int(i/2)

        if gl[iax] is None:
            argp.error(f'Gridlines on {"xyz"[iax]} axis are irregular'
                + ' and hence "{v}" cannot be used.')

        if i % 2 == 0:
            _i[i] = np.searchsorted(gl[iax],float(v),side='left')

        else:
            vv = float(v)
            if vv < gl[iax][0]:
                vv = gl[iax][-1]+vv
            _i[i] = np.searchsorted(gl[iax],vv,side='right')


    _islice = '[{}:{},{}:{},{}:{}]'.format(*_i)

    try:
        dout = eval('din'+_islice)
    except BaseException as e:
        argp.error(f'Error when processing --slice {args.slice}'
            +f', converted to indices {_islice}:\n'
            +str(e)
        )

    for iax,axsl in enumerate(_islice[1:-1].split(',')):
        if gl[iax] is None:
            continue
        glout[iax] = eval(f'gl[iax][{axsl}]')
else:
    dout = din.copy()
    glout = grid.get_grid_lines()

if args.interp_to:

    _glin = glout[:]
    _din = dout.copy()

    path, pfx, fnrest = PathToPrefix.split(args.interp_to)
    pathpfx = os.path.join(path,pfx)
    gridout = HGSGrid(pathpfx)

    if _glin[2] is not None:
        raise NotImplementedError('Only special case implemented')

    # Special case of 2-layer undulating z-dimesion
    # perform 2D interpolation in x-y and copy to all layers in z

    _glin = _glin[:2]
    _din = _din[:,:,0]
    _glout = gridout.get_grid_lines()[:2]

    ut,vt = np.meshgrid(*_glout, indexing='ij')
    ptsout = np.array([ut.ravel(), vt.ravel()]).T

    interp = RegularGridInterpolator(_glin, _din)
    dout = interp(ptsout).reshape(*(map(len,_glout)))

    dout = np.repeat(dout[:,:,np.newaxis], 2, axis=2)

if args.show:
    fig,ax = plt.subplots(1, 2, sharex=True, sharey=True)
    vmi = np.min(din[:,:,0])
    vma = np.max(din[:,:,0])
    ax[0].imshow(din[:,:,0].T, vmin=vmi, vmax=vma)
    ax[0].set_title('Input')
    ax[0].set_xlim([0,max(din.shape[0], dout.shape[0])])
    ax[0].set_ylim([0,max(din.shape[1], dout.shape[1])])

    ax[1].imshow(dout[:,:,0].T, vmin=vmi, vmax=vma)
    ax[1].set_title('Output')
    ax[1].set_xlim(ax[0].get_xlim())
    ax[1].set_ylim(ax[0].get_ylim())
    
    plt.tight_layout()
    plt.show()
    

with FortranFile(args.FILE_OUT,'w') as fout:
    fout.write_record(dout.ravel())
