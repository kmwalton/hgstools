#!/usr/bin/env python
"""Calculate bulk hydraulic conductivity (Kbulk) along axes of a HGS domain.

Run from within a HydroGeoSphere simulation directory, this command-line tool
estimates a Darcy bulk hydraulic conductivity for one or more axis-aligned
zones (bounding boxes) of the domain. For each zone and each active axis it:

  1. averages the simulated hydraulic head over the two opposing faces of the
     zone (volume-weighted, via ``AvCalc``) and divides their difference by the
     centroid-to-centroid distance to obtain the bulk gradient ``i``;
  2. integrates the Darcy flux over those faces to obtain the discharge ``Q``
     and the area-averaged flux ``q`` (via ``DischargeCalc``, which splits the
     porous-medium and fracture contributions); and
  3. combines them by Darcy's law, ``Kbulk = -(q0 + q1) / (2 i)``.

Results are printed as a table or, with ``--json``, emitted as JSON for
downstream post-processing. See ``main`` for the command-line interface.

It supports arbitrary 3D zones and any combination of the x/y/z axes.
"""

import sys
import os
import contextlib
import argparse
import json
from pathlib import Path
from datetime import datetime
from itertools import starmap, product, compress
from collections import namedtuple
from functools import partial
from math import log10,floor

import numpy as np

import hgstools
from hgstools.pyhgs.cli import parse_path_to_prefix
from hgstools.pyhgs.parser import parse as hgs_parse
from hgstools.pyhgs.parser import peek_NNNN_time
from hgstools.pyhgs.mesh import HGSGrid
from hgstools.pyhgs.calcs import AABBox, AvCalc, DischargeCalc


try:
    # Attempt the specific import
    from g360.scriptmeta import get_imprint
except (ImportError, ModuleNotFoundError):
    # Fallback function if the import fails
    def get_imprint(*args, **kwargs):
        """Fallback: returns current date and file name on two lines."""
        current_date = datetime.now().strftime('%Y-%m-%d')
        # __file__ gives the path to the script calling this function
        file_name = os.path.basename(__file__)
        return f"{current_date}\n{file_name}"

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
calcslogger = logging.getLogger('hgstools.pyhgs.calcs')
calcslogger.addHandler(logging.StreamHandler())


################################################################################
# helper functions
#

def _fmtlistarr(a, prefix='', **fmtopts):
    """Format a sequence of arrays as a prefixed, multi-line string.

    The arrays in `a` are stacked (`numpy.vstack`) and rendered with
    `numpy.array2string`; `prefix` (e.g. a tab) is prepended to every line so
    the block indents cleanly beneath a logging header. Extra `fmtopts` are
    passed through to `array2string`.
    """
    return '\n'.join( prefix + line
        for line in np.array2string(np.vstack(a),**fmtopts).splitlines())

# time conversion factors
convt2s = {
    'day':86400.,
    'd':86400.,
    'hour':3600.,
    'hr':3600.,
    'second':1.,
    's':1.,
}
"""Conversion factor from simulation time units to seconds"""

def calc_head_from_output(grid, filename, bnds):
    '''calculate average PM grid head at x-boundaries'''

    xs,xe,ys,ye,zs,ze = bnds
    (_time, _head_pm) = hgs_parse(filename).values()
    calc = AvCalc(grid)
    h0 = calc.average(
        f'{xs} {xs} {ys} {ye} {zs} {ze}',
        _head_pm,
        weight='vol_el')
    h1 = calc.average(
        f'{xe} {xe} {ys} {ye} {zs} {ze}',
        _head_pm,
        weight='vol_el')

    return h0, h1

def calc_distances(grid, zones):
    '''Return the element centroid-distances across opposing faces'''
    ret = []
    _gl = grid.get_grid_lines()
    for z in zones:
        cz = z.to_centroid_bbox(*_gl)
        ret.append(cz[3:]-cz[:3])
    return ret

def calc_avg_heads(grid, zones, dim_mask=3*(True,), timeidx=1):
    '''Return the average head for each face in each zone'''
    ret = []

    _eps = 10e-4
    _gl = grid.get_grid_lines()

    head_file = Path(f'{grid.ppfx}o.head_pm.{timeidx:04d}')
    _p = hgs_parse(head_file)
    sim_t, h_pm_nodal = _p['ts'], _p['data']

    logger.info(f'Using {head_file} with time {sim_t}')

    #heads = grid.get_element_vals(h_pm_nodal)
    heads = h_pm_nodal
    calc = AvCalc(grid, allow_partial=False)

    for zn in zones:
        # final calculations for each face stored here
        hbar = np.zeros(6)
        face_specs = []

        # iterate over opposite, active faces, and find their block
        # specifications
        #for ax,(_f0, _f1) in zip(np.arange(3)[dim_mask], zn.iter_face_bbox(dim_mask)):
        for ax,(_f0, _f1) in zip(np.arange(3)[dim_mask], zn.iter_layer_bbox(_gl, dim_mask)):
            # get the face in a format for the average calc x0 x1 y0 y1 z0 z1
            face0 = ' '.join(f'{v:.3f}' for v in _f0.reshape((2,3)).T.ravel())
            face1 = ' '.join(f'{v:.3f}' for v in _f1.reshape((2,3)).T.ravel())
            face_specs.append((ax,face0,face1,))

        logger.debug(f'For zone {zn}\n\tUsing grid layers (x0 x1 y0 y1 z0 z1):\n\t  '
            +'\n\t  '.join(f'{f0}\n\t  {f1}' for ax,f0,f1 in face_specs)
        )

        # Now, iterate over the pairs of faces and do the calculation
        # This 
        for ax,face0,face1 in face_specs:
            # What should the weighting be here?
            # Should these be weighted by discharge magnitude, perhaps?
            hbar[2*ax] = calc.average(face0, heads, weight='vol_el')
            hbar[2*ax+1] = calc.average(face1, heads, weight='vol_el')

        ret.append(hbar)

    return ret

def calc_i(zones, distances, heads, mask):
    '''Return the bulk hydraulic gradient vector for each zone.

    Along each axis the gradient is the head difference between the zone's two
    opposing faces divided by the centroid-to-centroid distance across the
    zone, ``i = (h_face1 - h_face0) / d``. Only axes selected by `mask` are
    computed; the rest are left zero.

    Parameters:
        zones : list of AABBox
        distances : per-zone centroid-to-centroid distances ``(dx, dy, dz)``,
            from `calc_distances`
        heads : per-zone face heads
            ``(h_x0, h_x1, h_y0, h_y1, h_z0, h_z1)``, from `calc_avg_heads`
        mask : boolean array selecting the active (x, y, z) axes
    '''
    ret = []
    for zn, d, h in zip(zones, distances, heads):
        i = np.zeros(3)
        i[mask] = (h[1::2]-h[::2])[mask] / d[mask]
        ret.append(i)
    return ret

def _max_f_width(v):
    '''Determine the maximum integer printing-width of the inputs v'''

    def _w(vv):
        sigfigs = floor(log10(abs(vv)))+1 if vv > 0. else 1
        sign = vv+0.0 < 0.
        return sign+sigfigs

    return max(map(_w, v))

def calc_q(grid, zones, dim_mask=3*(True,), timeidx=1):
    '''Return area and discharge over each face of each zone'''

    nudge = 0.001

    retA = []
    retq = []

    calc = DischargeCalc(grid)
    _gl = grid.get_grid_lines()

    _wid = _max_f_width(_gl[a][i] for a,i in product([0,1,2],[0,-1]))

    flux_pm_file = f'{grid.ppfx}o.q_pm.{timeidx:04}'
    flux_fx_file = f'{grid.ppfx}o.v_frac.{timeidx:04}'
    sim_t = peek_NNNN_time(flux_pm_file)
    logger.info(f'Using simulation time={sim_t} q-data from {flux_pm_file} and {flux_fx_file}')

    for zn in zones:
        # final calculations for each face stored here
        q = np.zeros(6)
        A = np.zeros(3)

        face_specs = []

        # Area index, 0=total, 1=PM only, 2=Fx only
        _Aidx = 1

        # iterate over opposite, active faces
        for ax,(_f0, _f1) in zip(np.arange(3)[dim_mask], zn.iter_layer_bbox(_gl, dim_mask)):
            # get the face in a format for the average calc x0 x1 y0 y1 z0 z1
            face0 = ' '.join(f'{v:{_wid+5}.4f}' for v in _f0.reshape((2,3)).T.ravel())
            face1 = ' '.join(f'{v:{_wid+5}.4f}' for v in _f1.reshape((2,3)).T.ravel())
            face_specs.append((ax, face0, face1,))

        logger.debug(f'For zone {zn}\n\tUsing grid layers (x0 x1 y0 y1 z0 z1):\n\t  '
            +'\n\t  '.join(f'{f0}\n\t  {f1}' for ax,f0,f1 in face_specs)
        )

        Qin = Qout = 0

        # do the calculations
        for ax,face0,face1 in face_specs:
            t, f0A, f0Q = calc.discharge_at(face0, ax, 1)
            q[2*ax] = f0Q[0] / f0A[_Aidx]

            t, f1A, f1Q = calc.discharge_at(face1, ax, 1)
            q[2*ax+1] = f1Q[0] / f1A[_Aidx]

            A[ax] = (f0A[_Aidx]+f1A[_Aidx])/2

            logger.debug(f'Average Q{"xyz"[ax]}={(f0Q[0]+f1Q[0])/2:.4g}')

            Qin += f0Q[0]
            Qout += f1Q[0]

        Qnet = Qin - Qout
        if abs(Qnet)>1e-3:
            logger.warning(f'WARNING: Net discharge in zone {zn} is\n  Q_net = {Qnet:.3g} = {Qin:.3g} inflow - {-Qout:.3g} outflow')
        else:
            logger.debug(f'Net discharge in zone {zn} is\n  Q_net = {Qnet:.3g} = {Qin:.3g} inflow - {-Qout:.3g} outflow')

        retA.append(A)
        retq.append(q)

    return retA, retq

def calc_Kbulk(zones, i, q, dim_mask=3*(True,)):
    '''Return the Darcy bulk conductivity vector for each zone.

    Applies Darcy's law per axis using the mean of the fluxes on the two
    opposing faces: ``Kbulk = -(q_face0 + q_face1) / (2 * i)``. Axes with a
    zero gradient produce inf/nan (the division is guarded and the result is
    left as-is for the caller to interpret).

    Parameters:
        zones : list of AABBox
        i : per-zone bulk gradient vectors, from `calc_i`
        q : per-zone face fluxes
            ``(q_x0, q_x1, q_y0, q_y1, q_z0, q_z1)``, from `calc_q`
        dim_mask : accepted for signature symmetry with the other ``calc_*``
            helpers; the computation is vectorized over all three axes, so
            inactive axes are simply ignored by callers rather than here.
    '''
    ret = []

    for zn,ii,qq in zip(zones, i, q):

        with np.errstate(divide='ignore', invalid='ignore'):
            Kbulk = -(qq[::2] + qq[1::2]) / ii / 2

        #Kbulk = np.zeros(3)
        # do by loop 
        #for ax in np.arange(3)[dim_mask]:
        #    Kbulk[ax] = -(qq[2*ax] + qq[2*ax+1]) / ii[ax] / 2

        ret.append(Kbulk)

    return ret

def as_dict(zones, A, Q, q, i, Kbulk, dim_mask=3*[True,]):
    '''Return the, zone-A-q-i-Kbulk data as a dictionary

    Parameters:
        Array of x, y, z values, per zone
        Q has x0 x1, y0, y1, z0, z1

    This is intended to be a convient format for exporting as json

    The 'zone' data can be converted back to a AABBox as

        from hgstools.pyhgs.calcs import AABBox
        # json_obj = ...
        a = AABBox(*json_obj['zone'])
    '''

    ret = dict()
    for iz,z in enumerate(zones):
        _d = { 'zone':z[:].tolist(), }
        _q = (q[iz][::2] + q[iz][1::2])/2
        for iax, ax in [(i, name) for i, (mask, name) in enumerate(zip(dim_mask, 'xyz')) if mask]:
            _d[ax] = dict()
            _d[ax]['A'] = A[iz][iax]
            _d[ax]['Q0'] = Q[iz][2*iax]
            _d[ax]['Q1'] = Q[iz][2*iax+1]
            _d[ax]['q'] = _q[iax]
            _d[ax]['i'] = i[iz][iax]
            _d[ax]['Kbulk'] = Kbulk[iz][iax]
        ret[iz] = _d
    return ret


#
# main script
#

def main():
    '''Parse command-line arguments and run the Kbulk calculation.

    Operates on the HGS simulation in the current working directory. Reads the
    zone bounding boxes and active dimensions from the arguments (see the
    ``argparse`` setup below), computes per-zone areas, discharges, fluxes,
    gradients and bulk conductivities, and writes a human-readable table to
    stdout or, when ``--json`` is given, a JSON document to the named file (or
    stdout).
    '''

    # Create the ArgumentParser
    parser = argparse.ArgumentParser('''Calculate a bulk gradient for arbitrary
            3D zones regions in a domain''')

    # Add Positional Arguments
    # By default, arguments added without a '-' or '--' prefix are positional
    # and required.

    parser.add_argument('-d', '--dims',
        default=7,
        type=int,
        help='''The dimensions to consider for Kbulk: 1 for x-, 2 for y-, and 4 for
        z-, where the selected dimensions are summed. E.g., for x- and z- only,
        use "-d 5". Default 7.''',
    )

    parser.add_argument(
        '-z', '--zone',
        nargs=6,
        metavar='V',
        #metavar="x0 y0 z0 x1 y1 z1".split(),
        type=float,
        dest='zones',
        action='append',
        help='A bounding box of a Kbulk zone, as six floats "x0 y0 z0 x1 y1 z1".'
    )

    parser.add_argument(
        '-v', '--verbose', 
        action='count', 
        default=0,
        help="Increase output verbosity (e.g., -v, -vv, -vvv)"
    )

    parser.add_argument(
    '--json', 
    metavar='JSON_FILE',
    nargs='?',           # Optional: 0 or 1 arguments
    const=sys.stdout,    # Value if --json is present but no file is provided
    default=None,        # Value if --json is not present at all
    help="Output as JSON to a file or stdout if '-' or no file is indicated"
    )

    args = parser.parse_args()

    # set loggers based on verbosity
    if args.verbose == 0:
        logger.setLevel(logging.WARNING)
        calcslogger.setLevel(logging.WARNING)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
        calcslogger.setLevel(logging.WARNING)
    elif args.verbose == 2:
        logger.setLevel(logging.DEBUG)
        calcslogger.setLevel(logging.INFO)
    elif args.verbose >= 3:
        logger.setLevel(logging.DEBUG)
        calcslogger.setLevel(logging.DEBUG)

    # parse the zones; overwrite the args.zones
    if args.zones:
        bboxes = []
        for box_values in args.zones:
            box = AABBox(*box_values)

            # Validation Logic
            if not (box.x1 > box.x0 and box.y1 > box.y0 and box.z1 > box.z0):
                parser.error(
                    f"Invalid AABBox {box_values}: Max coordinates must be "
                    f"greater than Min coordinates (x1>x0, y1>y0, z1>z0)."
                )
            bboxes.append(box)
        args.zones = bboxes

    # Post-processing to handle the '-' case
    if args.json == '-':
        args.json = sys.stdout

    np.set_printoptions(linewidth=np.inf, sign=' ', precision=3)

    # get simulation general data
    pth, pfx, _ = parse_path_to_prefix('.')
    grid = HGSGrid(pfx)
    _gl = grid.get_grid_lines()
    logger.debug(f'Found sim {pfx}')
    if not args.zones:
        def glgetter(j,i): return _gl[i][j]
        args.zones=[ AABBox(*starmap(glgetter, product([0,-1],[0,1,2]))), ]

    # create a mask for dimensions with bitwise operation and the dimensions of
    # the domain --- ignore any HGS domain dimensions that appear to be a plane
    dim_mask = list( (args.dims & v) and len(gla)>2 for v,gla in zip((1,2,4),_gl) )
    dim_mask = np.array(dim_mask, dtype=bool)

    eco = hgs_parse(pfx+'o.eco')
    unitM,unitL,unitt = eco.get_units()

    # some units for printout
    _ux = 'x y z'
    _udx = 'dx dy dz'
    _uyz = 'yz xz xy'
    _ux0x1 = " ".join(f"{a}{i}" for a,i in product("xyz",(0,1)))

    logger.info('Found dims '+' '.join(s for s,t in zip('xyz',dim_mask) if t))
    logger.info('Found zones:\n\t'+ '\n\t'.join(str(zn) for zn in args.zones))
    logger.info(f'Found units {unitM}-{unitL}-{unitt}')
    logger.info('\n'+80*'-'+'\n')

    shrink_val = 0.0
    if shrink_val != 0.0:
        logger.info(f'Modifying zones by shrinking by {shrink_val}')
        args.zones = [z.shrink(shrink_val) for z in args.zones]
        logger.info('Using zones:\n\t'+ '\n\t'.join(str(zn) for zn in args.zones))
        logger.info('\n'+80*'-'+'\n')

    distances = calc_distances(grid, args.zones)
    logger.info(f'Centroid-Centroid Distances ({_udx}):\n\t'+ '\n\t'.join(str(d) for d in distances))
    logger.info('\n'+80*'-'+'\n')

    areas, fluxes = calc_q(grid, args.zones, dim_mask)

    _opt = np.get_printoptions()
    A = np.vstack(areas)
    Q = np.stack((A,A),axis=-1).reshape(A.shape[0],-1) * np.vstack(fluxes)
    np.set_printoptions(precision=3,floatmode='fixed')
    logger.info(f'Plane areas ({_uyz}):\n'+_fmtlistarr(areas,'\t'))
    np.set_printoptions(sign=' ', precision=3, floatmode='fixed')
    logger.info(f'Average discharge (Q at planes {_ux0x1}):\n'+_fmtlistarr(Q,'\t'))
    logger.info(f'Average fluxes (q at planes {_ux0x1}):\n'+_fmtlistarr(fluxes,'\t'))
    logger.info('\n'+80*'-'+'\n')
    np.set_printoptions(**_opt)
    del _opt

    heads = calc_avg_heads(grid, args.zones, dim_mask)
    logger.info(f'Average heads (h at {_ux0x1}):\n'+_fmtlistarr(heads,'\t'))
    logger.info('\n'+80*'-'+'\n')

    gradients = calc_i(args.zones, distances, heads, dim_mask)
    logger.info(f'Average gradients (i_ {_ux}):\n'+_fmtlistarr(gradients,'\t'))
    logger.info('\n'+80*'-'+'\n')

    Kbulk = calc_Kbulk(args.zones, gradients, fluxes, dim_mask)
    logger.info(f'\nKbulk results ({unitL}/{unitt}):\n'+_fmtlistarr(Kbulk, '\t'))
    logger.info('\n'+80*'-'+'\n')

    # do final printouts
    pr = partial(print, file=sys.stdout)
    #np.set_printoptions(formatter={'float':'{:-.4f}'})
    np.set_printoptions(sign=' ',precision=4,floatmode='fixed')

    imprint = get_imprint(__file__, True)

    if args.json is None:
        pr('\n'+imprint)
        pr(f'\nKbulk results in {" ".join(v for v in compress("xyz",dim_mask))} ({unitL}/{unitt}):')
        for zn,Kb in zip(args.zones, Kbulk):
            pr(f'{str(Kb[dim_mask]):32} {zn}')

        if unitt != 's':
            pr(f'\nKbulk ({unitL}/s):')
            pr('\n'.join(f'{str(Kb[dim_mask]/convt2s[unitt]):32}'
                for zn,Kb in zip(args.zones, Kbulk)))

    else:
        _d = { 'meta':imprint, 'units':eco.get_units(), }
        _d.update(as_dict(args.zones, A, Q, fluxes, gradients, Kbulk, dim_mask))
        json_txt = json.dumps(_d, indent=2)

        # 1. Determine our context manager
        if isinstance(args.json, str) and args.json != '-':
            ctx = open(args.json, 'w')
        else:
            # nullcontext just "passes through" the object (stdout) 
            # and does NOT close it when the block ends.
            ctx = contextlib.nullcontext(sys.stdout)

        with ctx as fout:
            print(json_txt, file=fout)

 

if __name__ == '__main__':
    main()
    exit(0)
