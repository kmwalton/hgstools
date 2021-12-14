"""A collection of parsing functions for Hydrogeosphere binary-format files

The most useful function in here is `parse`, which attempts to determine the
correct parser to use based on the file name passed to it.

All parser functions return `dict` objects, where keys in the dictionary give
hints at the sections of types of data parsed out of each file. See
documentation for specialized `parse_*` functions, or set this module's logger
to INFO or INFO-1 level to identify and spew the documentation of which parsers
are used.

Example:
    Parse and store PM coordinates, and check which parser was used:

>>> from pyhgs.parser import parse
>>> parse_logger = logging.getLogger('pyhgs.parser')
>>> parse_logger.addHandler(logging.StreamHandler())
>>> parse_logger.setLevel(logging.INFO-2)
>>> pm_data = parse('my_simo.coordinates_pm')

Specification of binary file structure provided by Killian Miller
<kmiller@aquanty.com> to Ken Walton <kmwalton@g360group.org> by email on Feb.
24, 2021, and July 26, 2021.
"""
import os
import tempfile
import logging
from collections import OrderedDict

import numpy as np
from scipy.io import FortranFile

__docformat__ = 'numpy'

logger = logging.getLogger(__name__)

def parse_coordinates_pm(fn):
    """Parse o.coordinates_pm file and return a dict.

    Arguments:
        fn : str
            A filename, or the simulation prefix, or a directory and prefix.

    Returns:
        A dict with the following keywords

        1. nn: Total number of nodes (int4).
        2. ncoords: x(i), y(i), z(i), i=1,nn, Nodal xyz-coordinates for each
           node (real8).
        3. nx, ny, nz, nsptot: Number of grid lines in the x-, y-, and
           z-coordinates and the number species for transport (int4).
        4. ne2d: Number of triangles or rectangles in a 2-D node sheet (int4).
        5. tetramesh: Logical flag that is true if the mesh consists of
           tetrahedrons.
        6. gal_mode: Logical flag that is true if Galerkin mode is used for
           tetrahedral mesh.
        7. write_face_seg: Logical flag.
        8. nb2d: Maximum number of nodes connected to a node for a 2-D
           triangular mesh (int4).
    """

    if not os.path.exists(fn):
        fn = f'{fn}o.coordinates_pm'

    with FortranFile(fn,'r') as fin:
        nn = fin.read_ints()[0]
        coords = fin.read_reals()
        (nx,ny,nz,nsptot) = fin.read_ints()[0:4]
        ne2d = fin.read_ints()[0]
        tetramesh = fin.read_ints(np.bool_)[0]
        gal_mode = fin.read_ints(np.bool_)[0]
        write_face_seg = fin.read_ints(np.bool_)[0]
        nb2d = fin.read_ints()[0]

    coords.resize(nn,3)

    return OrderedDict( [
            ('nn',nn),
            ('ncoords',coords),
            ('nx',nx),
            ('ny',ny),
            ('nz',nz),
            ('nsptot',nsptot),
            ('ne2d',ne2d),
            ('tetramesh',tetramesh),
            ('gal_mode',gal_mode),
            ('write_face_seg',write_face_seg),
            ('nb2d',nb2d),
        ] )

def parse_elements_pm(fn):
    """Parse o.elements_pm file and return a dict.

    Arguments:
        fn : str
            A filename, or the simulation prefix, or a directory and prefix.

    Returns:
        A dict with the following keywords:

        1. nln: Number of nodes per element (6: triangular prism, 8: hexahedron)
           (int4).
        2. ne: Total number of elements (int4).
        3. inc(i, j) i=1,nln, j=1,ne: Node numbers of elemental incidences for
           each element (int4).
           *Node numbers are returned with 0-based indicies.*
        4. zone(i), i=1,ne: Element ID (zone number) for each element (int4).
    """

    if not os.path.exists(fn):
        fn = f'{fn}o.elements_pm'

    with FortranFile(fn,'r') as fin:
        nln = fin.read_ints()[0]
        ne = fin.read_ints()[0]
        inc = fin.read_ints()
        zone = fin.read_ints()

    inc.resize(ne,nln)
    np.add(inc,-1,out=inc)

    return OrderedDict( [
            ('nln',nln),
            ('ne',ne),
            ('inc',inc),
            ('zone',zone),
        ] )

def parse_coordinates_frac(fn):
    """Parse the o.coordinates_frac file and return a `dict`.

    Returns:
        A `dict` with the following keys/values.
        1. nnfrac: Total number of fracture nodes (int4)
        2. frac_scheme: 0 for common node or 1 for dual node (int4)
        3. link_frac2pm(i), i=1,nnfrac: Fracture node to PM node mapping (int4)
            *0-based indices returned.*
        4. link_pm2frac(i), i=1,nn: PM node to fracture node mapping (int4)
            *0-based indices returned.*
    """

    if not os.path.exists(fn):
        #assume a prefix or path/prefix given
        fn = f'{fn}o.coordinates_frac'

    with FortranFile(fn,'r') as fin:
        nnfrac = fin.read_ints()[0]
        frac_scheme = fin.read_ints()[0] 
        link_frac2pm = fin.read_ints()
        link_pm2frac = fin.read_ints()

    np.add(link_frac2pm,-1,out=link_frac2pm)
    np.add(link_pm2frac,-1,out=link_pm2frac)

    return OrderedDict( [
            ('nnfrac',nnfrac),
            ('frac_scheme',frac_scheme),
            ('link_frac2pm',link_frac2pm),
            ('link_pm2frac',link_pm2frac),
        ] )

def parse_elements_frac(fn):
    """Parse o.elements_frac file and return a `dict`.

    Returns:
        A `dict` with the following values
        1. nln: number of nodes per frac element (int4)
        2. nfe: total number of fracture elements (int4)
        3. inc(i,j), i=1,nln, j=1,nfe:
            node numbers of fracture element incidences (int4)
            *0-based indices returned.*
        4. zone(i), i=1,nfe: zone number for each fracture element (int4)
        5. face_map(i, j), i=1,nfe, j=1,2: fracture to element face mapping (int4)
        6. ap(i), i=1,nfe: aperture values (real8)
    """

    if not os.path.exists(fn):
        fn = f'{fn}o.elements_frac'

    with FortranFile(fn,'r') as fin:
        nln = fin.read_ints()[0]
        nfe = fin.read_ints()[0]
        inc = fin.read_ints()
        zone = fin.read_ints()
        face_map = fin.read_ints()
        ap = fin.read_reals()

    inc = np.reshape(inc,(nln,nfe),order='F').T
    #inc = np.reshape(inc, (nfe,nln)) # equivalent?
    np.add(inc,-1,out=inc)

    face_map.resize(nfe,2)

    return OrderedDict( [
            ('nln',nln),
            ('nfe',nfe),
            ('inc',inc),
            ('zone',zone),
            ('face_map',face_map),
            ('ap',ap),
        ] )


def _parse(fn, dtype, shape=None):
    """Return the timestamp and [possibly rehsaped] data in a dict.

        Arguments:
            fn : str
                filename to open

            dtype : numpy.datatype
                interperet the data as this type

            shape : tuple
                resize the data to a numpy.ndarray of this size, if provided
    """
    with FortranFile(fn,'r') as fin:
        ts = fin.read_ints(dtype=np.byte)
        d = fin.read_reals(dtype=dtype)

    if shape:
        d = d.reshape(shape)

    return OrderedDict( [
            ('ts',ts.tobytes().decode('UTF-8').strip()),
            ('data',d),
        ] )

def parse_1D_real8(fn):
    """(a) 1D real8 fields

File format: [int4][char80][int4][int4][real8]...[real8][int4]

o.ETEvap_olf.XXXX 
o.ETPmEvap3D_pm.XXXX
o.ETPmEvap_olf.XXXX
o.ETPmTranspire3D_pm.XXXX
o.ETPmTranspire_olf.XXXX
o.ETTotal_olf.XXXX
o.conc_dual.species.XXXX
o.conc_frac.species.XXXX
o.conc_olf.species.XXXX
o.conc_pm.species.XXXX
o.freeze_thaw_temp_pm.XXXX
o.head_chan.XXXX
o.head_dual.XXXX
o.head_frac.XXXX
o.head_olf.XXXX
o.head_pm.XXXX
o.head_well.XXXX
o.iconc_pm.species.XXXX
o.pet_olf.XXXX
o.rain_olf.XXXX
    """
    return _parse(fn, np.double)

def parse_1D_real4(fn):
    """(b) 1D real4 fields

File format: [int4][char80][int4][int4][real4]...[real4][int4]

o.ExchFlux_chan.XXXX 
o.ExchFlux_dual.XXXX
o.ExchFlux_olf.XXXX
o.ExchFlux_olf2_chan.XXXX
o.ExchFlux_pm2_chan.XXXX
o.ExchFlux_well.XXXX
o.ExchSolAdv_olf.species.XXXX
o.ExchSolDisp_olf.species.XXXX
o.exchsol_dual.species.XXXX
o.ice_sat_pm.XXXX
o.kxx.XXXX
o.kyy.XXXX
o.kzz.XXXX
o.por_pm.XXXX
o.sat_dual.XXXX
o.sat_frac.XXXX
o.sat_pm.XXXX
    """
    return _parse(fn, np.float)

def parse_2D_real8(fn):
    """(c) 2D real8 fields

File format: [int4][char80][int4][int4][real8][real8][int4]...[int4][real8][real8][int4]

o.friction_olf.XXXX 

    """
    raise NotImplementedError()

def parse_3D_real4(fn):
    """(d) 3D real4 fields

File format: [int4][char80][int4][int4][real4][real4][real4][int4]...[int4][real4][real4][real4][int4]

o.q_dual.XXXX 
o.q_pm.XXXX
o.tvk_pm.XXXX
o.v_chan.XXXX
o.v_dual.XXXX
o.v_frac.XXXX
o.v_olf.XXXX
o.v_pm.XXXX
o.v_well.XXXX
    """
    raise NotImplementedError()

def _is_fs_case_sensitive():
    """https://newbedev.com/check-if-file-system-is-case-insensitive-in-python
        Accessed 2021/07/26.
    """
    if not hasattr(_is_fs_case_sensitive, 'case_sensitive'):
        with tempfile.NamedTemporaryFile(prefix='TmP') as tmp_file:
            setattr(_is_fs_case_sensitive,
                    'case_sensitive',
                    not os.path.exists(tmp_file.name.lower()))
    return(_is_fs_case_sensitive.case_sensitive)

_file_to_parser_dict = {
    'o.coordinates_pm':parse_coordinates_pm,
    'o.elements_pm':parse_elements_pm,
    'o.coordinates_frac':parse_coordinates_frac,
    'o.elements_frac':parse_elements_frac,

    'o.ETEvap_olf':parse_1D_real8, 
    'o.ETPmEvap3D_pm':parse_1D_real8,
    'o.ETPmEvap_olf':parse_1D_real8,
    'o.ETPmTranspire3D_pm':parse_1D_real8,
    'o.ETPmTranspire_olf':parse_1D_real8,
    'o.ETTotal_olf':parse_1D_real8,
    'o.conc_dual':parse_1D_real8,
    'o.conc_frac':parse_1D_real8,
    'o.conc_olf':parse_1D_real8,
    'o.conc_pm':parse_1D_real8,
    'o.freeze_thaw_temp_pm':parse_1D_real8,
    'o.head_chan':parse_1D_real8,
    'o.head_dual':parse_1D_real8,
    'o.head_frac':parse_1D_real8,
    'o.head_olf':parse_1D_real8,
    'o.head_pm':parse_1D_real8,
    'o.head_well':parse_1D_real8,
    'o.iconc_pm':parse_1D_real8,
    'o.pet_olf':parse_1D_real8,
    'o.rain_olf':parse_1D_real8,

    'o.ExchFlux_chan':parse_1D_real4, 
    'o.ExchFlux_dual':parse_1D_real4,
    'o.ExchFlux_olf':parse_1D_real4,
    'o.ExchFlux_olf2_chan':parse_1D_real4,
    'o.ExchFlux_pm2_chan':parse_1D_real4,
    'o.ExchFlux_well':parse_1D_real4,
    'o.ExchSolAdv_olf':parse_1D_real4,
    'o.ExchSolDisp_olf':parse_1D_real4,
    'o.exchsol_dual':parse_1D_real4,
    'o.ice_sat_pm':parse_1D_real4,
    'o.kxx':parse_1D_real4,
    'o.kyy':parse_1D_real4,
    'o.kzz':parse_1D_real4,
    'o.por_pm':parse_1D_real4,
    'o.sat_dual':parse_1D_real4,
    'o.sat_frac':parse_1D_real4,
    'o.sat_pm':parse_1D_real4,

    'o.q_dual':parse_3D_real4, 
    'o.q_pm':parse_3D_real4,
    'o.tvk_pm':parse_3D_real4,
    'o.v_chan':parse_3D_real4,
    'o.v_dual':parse_3D_real4,
    'o.v_frac':parse_3D_real4,
    'o.v_olf':parse_3D_real4,
    'o.v_pm':parse_3D_real4,
    'o.v_well':parse_3D_real4,
}

def parse(fn):
    """Find a parser based on the given file name"""

    p = None

    # use lower case if filesystem is case insensitive
    _lowerfunc = lambda s: s
    if not _is_fs_case_sensitive:
        _lowerfunc = lambda s: s.lower()

    # linear search through _file_to_parser_dict
    # ick.
    fn = _lowerfunc(fn)
    for pat in _file_to_parser_dict:
        if _lowerfunc(pat) in fn:
            p = _file_to_parser_dict[pat]
            break

    if p:
        logger.info(f'Parsing {fn} using {p.__name__}')
        logger.log(logging.INFO-1,p.__doc__)
        return p(fn)
    else:
        raise RuntimeError(f'No parser for {fn}')

