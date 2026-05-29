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

Notes
-----
The `1D`, `2D`, and `3D` in the function names below refer to the number of
dimensions of each datum, **not** the spatial dimensions of the HGS domain.
E.g. concentration files must be parsed with the `1D` reader, and velocity
files must be parsed with the `3D` reader.

Binary file header format
-------------------------
All HGS binary output (.NNNN) files begin with a Fortran record containing an
80-character header string.  Two layouts exist:

**Legacy format** (prior to HGS r2853, July 2025)::

    Bytes  1–80  simulation timestamp, space-padded on the right

**New format** (HGS r2853+, July 2025)::

    Bytes  1–55  simulation timestamp, space ' ', metadata
    Bytes 56–60  file magic number ``#!HGS``
    Bytes 61–63  file format version number (e.g. ``001``)
    Bytes 64–65  endianness: ``LE``, ``BE``, or ``UK`` (unknown)
    Bytes 66–68  data type: ``I04``, ``R04``, ``R08``, or ``UNK``
    Bytes 69–70  number of components per value, zero-padded (e.g. ``01``)
    Bytes 71–80  number of values in the data record, zero-padded

Parsers that return a ``dict`` expose these fields when the new header is
detected: ``hgs_file_ver``, ``endian``, ``data_type``, ``n_components``,
``n_values``.  For legacy files all five fields are ``None``.

Source: Aquanty HGS r2853 release notes (aquanty.com/updates, July 2025) and
Aquanty webinar, January 2026.

"""
import os
import tempfile
import logging
import enum
from collections import OrderedDict, defaultdict

import numpy as np
from scipy.io import FortranFile, FortranEOFError

from .lst import LSTFile
from .eco import EcoFile
from .mprops import parse as MPropsFile
from .grok import parse as GrokFile

__docformat__ = 'numpy'
# see https://numpydoc.readthedocs.io/en/latest/format.html

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARN)
_console = logging.StreamHandler()
_console.setLevel(logging.WARN)
_console.setFormatter(
    logging.Formatter(f'pyhgs.parser %(levelname)s - %(message)s'))
logger.addHandler(_console)

_NEW_HEADER_MAGIC = '#!HGS'

def _decode_header_str(raw):
    """Decode a raw 80-byte Fortran header record into ``(ts, ver)``.

    Parameters
    ----------
    raw : bytes-like
        The 80 raw bytes of the header record — either a ``numpy`` array
        (e.g. from ``FortranFile.read_ints(dtype=np.byte)`` or
        ``np.fromfile(..., dtype=[..., ('time_str','a80'), ...])[...]``) or a
        plain ``bytes`` object.

    Returns
    -------
    ts : str
        The simulation timestamp token, stripped of surrounding whitespace.
    ver : str or None
        The version/magic token (e.g. ``'#!HGS001LER08010000001600'``) when
        present (HGS r2853+), otherwise ``None``.
    """
    if hasattr(raw, 'tobytes'):
        raw = raw.tobytes()
    text = raw.decode('UTF-8', errors='replace').strip()
    tokens = text.split()
    if len(tokens) == 2 and tokens[1].startswith(_NEW_HEADER_MAGIC):
        return tokens[0], tokens[1]
    return text, None

def _decode_new_header(ver_token):
    """Decode structured sub-fields from a new-format (r2853+) ``ver`` token.

    Parameters
    ----------
    ver_token : str or None
        The second whitespace-delimited token from the 80-char header record,
        as extracted by ``_parse`` / ``_parse_nd``.  Expected to start with
        ``#!HGS`` for new-format files.

    Returns
    -------
    dict or None
        ``None`` when *ver_token* is absent or not a new-format token.
        Otherwise a dict with keys:

        hgs_file_ver : str
            Three-character file format version (e.g. ``'001'``).
        endian : str
            Two-character endianness code: ``'LE'``, ``'BE'``, or ``'UK'``.
        data_type : str
            Three-character numeric type: ``'I04'``, ``'R04'``, ``'R08'``,
            or ``'UNK'``.
        n_components : int
            Number of scalar components per data value (e.g. 1 for scalars,
            3 for velocity vectors).
        n_values : int
            Total number of data values in the file's data record.
    """
    if not ver_token or not ver_token.startswith(_NEW_HEADER_MAGIC):
        return None
    # Token layout (0-based within the token string):
    #  [0:5]  '#!HGS'  magic
    #  [5:8]  file format version (3 chars)
    #  [8:10] endianness (2 chars)
    # [10:13] data type (3 chars)
    # [13:15] number of components, zero-padded (2 chars)
    # [15:25] number of values, zero-padded (10 chars)
    return {
        'hgs_file_ver': ver_token[5:8],
        'endian':       ver_token[8:10],
        'data_type':    ver_token[10:13],
        'n_components': int(ver_token[13:15]),
        'n_values':     int(ver_token[15:25]),
    }

def parse_coordinates_pm(fn):
    """Parse `o.coordinates_pm` file and return a dict.

Arguments
---------
fn : str
    A filename, or the simulation prefix, or a directory and prefix.

Returns
-------
A `dict` with the following keywords

nn: int
    Total number of nodes (int4).
ncoords: array
    Nodal xyz-coordinates for each node, i (real8).
    i=[0..nn) of [x_i, y_i, z_i]
nx, ny, nz, nsptot: int
    Number of grid lines in the x-, y-, and z-coordinates and the number species
    for transport (int4).
ne2d: int
    Number of triangles or rectangles in a 2-D node sheet (int4).
tetramesh: bool
    Logical flag that is true if the mesh consists of tetrahedrons.
gal_mode: bool
    Logical flag that is true if Galerkin mode is used for tetrahedral mesh.
write_face_seg: bool
    Logical flag.
nb2d: int
    Maximum number of nodes connected to a node for a 2-D triangular mesh
    (int4).
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

Arguments
---------
fn : str
    A filename, or the simulation prefix, or a directory and prefix.

Returns
-------
A `dict` with the following keywords:

nln: int
    Number of nodes per element (6: triangular prism, 8: hexahedron) (int4).
ne: int
    Total number of elements (int4).
inc(i, j): numpy.ndarray
    i=[0..ne), j=[0..nln)
    Node numbers of elemental incidences for each element (int4). *Node numbers
    are returned with 0-based indicies.*
zone(i): numpy.ndarray
    i=[0..ne)
    Element ID (zone number) for each element (int4). Zone numbers retain their
    1-based indexing.
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
    """Parse the *o.coordinates_frac* file and return a `dict`.

    **Note:**
    Fracture index values (in the `link_pm2frac` array) reflect all-domain
    node ordering. For meshes that involve only porous media and fracture domains,
    subtracting the number of porous media nodes from these fracture indices
    will yield the index relative to only the fracture domain (like the indices
    needed in the `link_frac2pm` array).

    Returns
    ---
    nnfrac : int
        Total number of fracture nodes (int4)
    frac_scheme : int
        0 for common node or 1 for dual node (int4)
    link_frac2pm(i), i=1,nnfrac : int
        Fracture node to PM node mapping (int4) *0-based indices returned.*
    link_pm2frac(i), i=1,nn : int
        PM node to fracture node mapping (int4) *0-based indices returned.*
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
    """Parse *o.elements_frac* file and return a `dict`.

    Returns
    -------
    nln: int
        number of nodes per frac element (int4)
    nfe: int
        total number of fracture elements (int4)
    inc(i,j), i=1,nfe, j=1,nln: numpy.ndarray
        node numbers of fracture element incidences (int4) *0-based indices
        returned.*
    zone(i), i=1,nfe: numpy.ndarray
        zone number for each fracture element (int4)
    face_map(i, j), i=1,nfe, j=1,2: numpy.ndarray
        fracture to element face mapping (int4)
    ap(i), i=1,nfe: numpy.ndarray
        aperture values (real8)
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


def _parse(fn, dtype, shape=None, with_timestamp=True):
    """Return the timestamp and [possibly rehsaped] data in a dict.

    `with_timestamp` causes this parser to read (or not read) a timestamp (byte
    array) at the beginning of the file. Errors will occur if a timestamp is
    requested but not present, and vice versa. A timestamp of `None` will be in
    the return `dict` if no timestamp is read.

    Parameters
    ---
    fn : str
        filename to open

    dtype : numpy.datatype
        interperet the data as this type

    shape : tuple
        resize the data to a numpy.ndarray of this size, if provided

    with_timestamp : bool
        Optionally, read a time stamp preceding the data. Default True. This
        **must** match the presence of a timestamp in the datafile.
    """
    ts, ver, d = None, None, None

    with FortranFile(fn,'r') as fin:
        if with_timestamp:
            ts, ver = _decode_header_str(fin.read_ints(dtype=np.byte))
        d = fin.read_reals(dtype=dtype)

    if shape:
        d = d.reshape(shape)

    header_info = _decode_new_header(ver)
    ret = OrderedDict([('ts', ts), ('ver', ver)])
    if header_info:
        ret.update(header_info)
    ret['data'] = d
    return ret

def _parse_nd(fn, dtype, shape=None):
    """Read a number of *n-tuple* data points from file

    E.g., This will read 3-tuple velocity vectors, (vx,vy,vz), for an
    unknown-length sequence of elemental values.

    This function is _much_ more efficient if the shape of the data is known a
    priori. Array copying for unknown-length final data is avoided in this case.

    Handles both the legacy AOS (Array of Structures) layout used by HGS prior
    to r2853 (July 2025) and the new SOA (Structure of Arrays) layout.  The
    format is detected from the ``#!HGS`` magic token in the file header.

    Parameters
    ---
    fn : str
        The data file name.

    dtype : `numpy.dtype`
        The numeric type of the data in the array.

    shape : tuple (n_fields, n_components), optional
        A hint at the shape of the output data, e.g. ``(n_elements, 3)`` for
        velocity.  Required for the fast path; if omitted the function reads
        until EOF (slower).

    Returns
    -------
        `numpy.ndarray` with shape (n_fields, n_components)
    """

    def _read_shape_unknown(fin):
        _d = fin.read_reals(dtype=dtype)
        while True:
            try:
                _rd = fin.read_reals(dtype=dtype)
            except FortranEOFError as e:
                break
            else:
                _d = np.vstack((_d,_rd,))

        return _d

    ts = None
    ver = None
    d = None

    if shape is None:
        with FortranFile(fn,'r') as fin:
            ts_raw = fin.read_ints(dtype=np.byte)
            d = _read_shape_unknown(fin)

        ts, ver = _decode_header_str(ts_raw)

        # SOA (modern): _read_shape_unknown built (n_components, n_fields);
        # AOS (legacy): it built (n_fields, n_components) — no transpose needed.
        is_soa = ver is not None
        ret_d = d.T if is_soa else d

    else:
        # Assumes HGS-"standard" offset of 80-character fortran array at the
        # beginning of the file.  Read the header first so the format (SOA vs
        # AOS) can be determined before reading the data records.
        with open(fn,'rb') as fin:
            ts_raw = np.fromfile(fin,
                    dtype=[('','i4'), ('time_str','a80',), ('','i4'), ],
                    count=1,)['time_str']

            ts, ver = _decode_header_str(ts_raw)
            is_soa = ver is not None
            if is_soa:
                # SOA: shape[1] records, each containing shape[0] values
                d = np.fromfile(fin,
                        dtype=[('','i4'), ('data', dtype, shape[0],), ('','i4',),],
                        count=shape[1])['data']
                ret_d = d.T
            else:
                # AOS: shape[0] records, each containing shape[1] values
                d = np.fromfile(fin,
                        dtype=[('','i4'), ('data', dtype, shape[1],), ('','i4',),],
                        count=shape[0])['data']
                ret_d = d
            logger.debug(f'Read {fin.tell()} bytes')

    header_info = _decode_new_header(ver)
    ret = OrderedDict([('ts', ts), ('ver', ver)])
    if header_info:
        ret.update(header_info)
    ret['data'] = ret_d
    return ret

def parse_1D_real8(fn, **kwargs):
    """(a) 1D real8 fields

    **Notes:**
    File expected file format

        [int4][char80][int4]
        [int4][real8]...[real8][int4]

   Used for:

    - *o.ETEvap_olf.NNNN*
    - *o.ETPmEvap3D_pm.NNNN*
    - *o.ETPmEvap_olf.NNNN*
    - *o.ETPmTranspire3D_pm.NNNN*
    - *o.ETPmTranspire_olf.NNNN*
    - *o.ETTotal_olf.NNNN*
    - *o.conc_dual.species.NNNN*
    - *o.conc_frac.species.NNNN*
    - *o.conc_olf.species.NNNN*
    - *o.conc_pm.species.NNNN*
    - *o.freeze_thaw_temp_pm.NNNN*
    - *o.head_chan.NNNN*
    - *o.head_dual.NNNN*
    - *o.head_frac.NNNN*
    - *o.head_olf.NNNN*
    - *o.head_pm.NNNN*
    - *o.head_well.NNNN*
    - *o.iconc_pm.species.NNNN*
    - *o.mass_sorbed_pm.NNNN*
    - *o.mass_sorbed_frac.NNNN*
    - *o.pet_olf.NNNN*
    - *o.rain_olf.NNNN*

    Parameters
    ---------
    fn : str
        File name

    kwargs : dict, optional
        Ignored.

    """
    return _parse(fn, np.float64)

def parse_1D_real8_nots(fn, **kwargs):
    """Like, `parse_1D_real8` but

    Used for:
    - *o.hen
    - *o.cen
    """
    return _parse(fn, np.float64, with_timestamp=False)

def parse_1D_real4(fn, **kwargs):
    """(b) 1D real4 fields

    Notes
    -----

    File format:

        [int4][char80][int4]
        [int4][real4]...[real4][int4]

    Used for:

    - *o.ExchFlux_chan.NNNN*
    - *o.ExchFlux_dual.NNNN*
    - *o.ExchFlux_olf.NNNN*
    - *o.ExchFlux_olf2_chan.NNNN*
    - *o.ExchFlux_pm2_chan.NNNN*
    - *o.ExchFlux_well.NNNN*
    - *o.ExchSolAdv_olf.species.NNNN*
    - *o.ExchSolDisp_olf.species.NNNN*
    - *o.exchsol_dual.species.NNNN*
    - *o.ice_sat_pm.NNNN*
    - *o.kxx.NNNN*
    - *o.kyy.NNNN*
    - *o.kzz.NNNN*
    - *o.por_pm.NNNN*
    - *o.sat_dual.NNNN*
    - *o.sat_frac.NNNN*
    - *o.sat_pm.NNNN*

    Parameters
    ---------
    fn : str
        File name

    kwargs : dict, optional
        Ignored.

    """
    return _parse(fn, np.float32)

def parse_2D_real8(fn, **kwargs):
    """(c) 2D real8 fields

Arguments
---------
    fn : str
        File name

    kwargs : dict, optional
        Ignored.

Notes
-----

File format:

    [int4][char80][int4]
    [int4][real8][real8]
    [int4]...
    [int4][real8][real8]
    [int4]

Used for:

    o.friction_olf.XXXX

    """
    raise NotImplementedError()

def parse_3D_real4(fn, **kwargs):
    """(d) 3D real4 fields

Arguments
---------
    fn : str
        File name

    count : int
        The number of 3-length vectors to be read

Notes
-----

File format:

    [int4][char80][int4]
    [int4][real4][real4][real4]
    [int4]...  # entries unknown
    [int4][real4][real4][real4]
    [int4]

Used for:

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

    if 'count' not in kwargs or kwargs['count'] is None:
        logger.warning(
                f'Parsing {fn} without "count". This might take a while.')
        return _parse_nd(fn,np.float32)

    _shp=(kwargs["count"],3,)
    logger.info(f'Parsing {fn} using shape={_shp}')
    return _parse_nd(fn,np.float32, shape=_shp)


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

class HGS_DATATYPE(enum.Enum):
    """Enumeration to describe types of data found in HGS output files"""
    UNSPEC = -1
    'Not yet specified'

    MULTI = 0
    'Combination of multiple sets of data'

    NODAL = 1
    'Data items for each node in the given domain'

    ELEMENTAL = 2
    'Data items for each element in the given domain'

_file_info = {
    'o.coordinates_pm'       :  (parse_coordinates_pm,    HGS_DATATYPE.MULTI,  ),
    'o.elements_pm'          :  (parse_elements_pm,       HGS_DATATYPE.MULTI,  ),
    'o.coordinates_frac'     :  (parse_coordinates_frac,  HGS_DATATYPE.MULTI,  ),
    'o.elements_frac'        :  (parse_elements_frac,     HGS_DATATYPE.MULTI,  ),
    'o.ETEvap_olf'           :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.ETPmEvap3D_pm'        :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.ETPmEvap_olf'         :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.ETPmTranspire3D_pm'   :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.ETPmTranspire_olf'    :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.ETTotal_olf'          :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.conc_dual'            :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.conc_frac'            :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.conc_olf'             :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.conc_pm'              :  (parse_1D_real8,          HGS_DATATYPE.NODAL,  ),
    'o.cen'                  :  (parse_1D_real8_nots,     HGS_DATATYPE.NODAL,  ),
    'o.freeze_thaw_temp_pm'  :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.head_chan'            :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.head_dual'            :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.head_frac'            :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.head_olf'             :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.head_pm'              :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.hen'                  :  (parse_1D_real8_nots,     HGS_DATATYPE.NODAL,  ),
    'o.head_well'            :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.iconc_pm'             :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.mass_sorbed_pm'       :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.mass_sorbed_frac'     :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.pet_olf'              :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.rain_olf'             :  (parse_1D_real8,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchFlux_chan'        :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchFlux_dual'        :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchFlux_olf'         :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchFlux_olf2_chan'   :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchFlux_pm2_chan'    :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchFlux_well'        :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchSolAdv_olf'       :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ExchSolDisp_olf'      :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.exchsol_dual'         :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.ice_sat_pm'           :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.kxx'                  :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.kyy'                  :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.kzz'                  :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.por_pm'               :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.sat_dual'             :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.sat_frac'             :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.sat_pm'               :  (parse_1D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.q_dual'               :  (parse_3D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.q_pm'                 :  (parse_3D_real4,          HGS_DATATYPE.ELEMENTAL,  ),
    'o.tvk_pm'               :  (parse_3D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.v_chan'               :  (parse_3D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.v_dual'               :  (parse_3D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.v_frac'               :  (parse_3D_real4,          HGS_DATATYPE.ELEMENTAL,  ),
    'o.v_olf'                :  (parse_3D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.v_pm'                 :  (parse_3D_real4,          HGS_DATATYPE.ELEMENTAL,  ),
    'o.v_well'               :  (parse_3D_real4,          HGS_DATATYPE.UNSPEC,  ),
    'o.lst'                  :  (LSTFile,                 HGS_DATATYPE.UNSPEC,  ),
    'o.eco'                  :  (EcoFile,                 HGS_DATATYPE.UNSPEC,  ),
    'mprops'                 :  (MPropsFile,              HGS_DATATYPE.UNSPEC,  ),
    'grok'                   :  (GrokFile,                HGS_DATATYPE.UNSPEC,  ),
}

def _find_file_info_entry(fn):
    # use lower case if filesystem is case insensitive
    _lowerfunc = lambda s: s
    if not _is_fs_case_sensitive:
        _lowerfunc = lambda s: s.lower()

    # linear search through FILE_TO_PARSER
    # ick.
    fn = _lowerfunc(fn)
    for pat,info in _file_info.items():
        if _lowerfunc(pat) in str(fn):
            return info

    raise RuntimeError(f'{fn} did not match any known HGS output file')

def get_datatype(fn):
    return _find_file_info_entry(fn)[1]

def parse(fn, **kwargs):
    """Find a parser based on the given file name

    Parameters
    ----------
    count : int, optional
        An integer denoting the quantity of scalar or vector data points in the
        file.  This will give the implementing parser function a hint at the
        size/shape of the data to expect in the file and may improve the
        performance of the data read by putting it in a container of an a
        priori-known size.

    """

    p = _find_file_info_entry(fn)[0]

    if p:
        logger.info(f'Parsing {fn} using {p.__name__}')
        logger.log(logging.DEBUG-1,p.__doc__)

        try:
            return p(fn, **kwargs)
        except BaseException as e:
            raise RuntimeError(f'When parsing {fn}') from e
    else:
        raise RuntimeError(f'No parser for {fn}')

def peek_NNNN_time(file_path):
    """Peek at the time string encoded in a prefixo.XXX.NNNN file

    Opens the file in binary mode, reads the first 88 bytes, and
    decodes a specific slice of those bytes as UTF-8 characters.

    Args:
        file_path (str): The path to the file.

    Returns:
        str: The decoded string, or None if an error occurs.
    """
    try:
        # Check if the file exists
        if not os.path.exists(file_path):
            print(f"Error: The file '{file_path}' does not exist.")
            return None

        # Open the file in binary read mode ('rb')
        with open(file_path, 'rb') as f:
            time_as_bytes = f.read(88)

            if len(time_as_bytes) < 88:
                print("Error: File is too small to read 88 bytes.")
                return None

            ts, _ver = _decode_header_str(time_as_bytes[4:84])
            return ts

    except Exception as e:
        print(f"An error occurred: {e}")
        return None
