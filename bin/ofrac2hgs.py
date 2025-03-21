#!/usr/bin/env python
"""Read DFN data and write HGS code for a discrete fracture grid.

Ken Walton
Feb 20, 2015
G360 Centre for Applied Groundwater Research

Upd: 14 Dec 2018, KW, output pickled data.
Upd: 18 Nov 2020, KW, new OFracGrid.merge interface.
Upd: 24 Nov 2020, KW, fixes to --force-extra-grid-lines
Upd: 02 Aug 2021, KW, fixes to --refine-near-fx-plane parsing
upd: 04 Nov 2021, KW, new OneLayerRFG
"""

import sys
import re
import datetime
import argparse
import pickle
import warnings
import decimal
from decimal import Decimal
from math import sqrt, floor, log10
from functools import reduce, partial
from collections import defaultdict

import numpy as np

# highly unconventional:
from hgstools.pyhgs._import_helpers import (
        _get_ofracs_module, _get_from_ofracs,)
ofracs_module = _get_ofracs_module()
locals().update(_get_from_ofracs(
   'D_AP', 'D_CO', 'D', 'OFrac', 'OFracGrid'))

from hgstools.pyhgs.gridgen_tools import printHGS

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.WARNING)
logger.verbose1 = partial(logger.log, logging.INFO)
logger.verbose2 = partial(logger.log, logging.INFO-1)




WISHLIST="""

Have RFG inherit from the OFracGrid object (or don't have an RFG class at all)

"""

# handy REs
_num_re_str = \
        r'[+-]?' \
        r'(?:(?:\d+\.\d*)|(?:\d*\.\d+)|(?:\d+))' \
        r'(?:[dDeE][+-]?\d+)?'

_num_re = re.compile(_num_re_str)

_gl_list_re_str = \
    r'\s*(?:(?P<axis>[xyzXYZ])\s*=\s*)?' \
    r'(?P<vals>(?:\s*'+_num_re_str+r'\s*,?){1,})\s*;?'

_gl_list_re = re.compile(_gl_list_re_str)

def apQuantize(v,n_sig):
    """Determine the aperture to n_sig significant figures."""

    #new_prec = Decimal('0.'+n_sig*'1')
    new_prec = Decimal(f'{10**(-n_sig):.{n_sig}f}')
    old_prec = D_AP(new_prec)

    if Decimal(v) < new_prec:
        return D_AP(new_prec/2)

    return Decimal(v).quantize(new_prec).quantize(old_prec)

class RFG:
    """Parse/store/manipulate/produce HGS-style grids and fracture definitions

    TODO: refactor and use more code from 'ofracs' module
    """

    def __init__(self,
            fnin,
            collapse_policy='warn-omit',
            domainSize=[],
            translate=None,
            forcedGridLines=[],
            nudgeTo=0.0,
            pmRefNearFx=[],
            regGlSpacing=[],
            maxGlSpacing=[]
        ):
        """fnin can be a string, of an iterable containing strings"""

        if type(fnin) == str:
            fnin = [ fnin, ]
        self.fnin = ', '.join(fnin)

        # merge input grids/dfns
        gridGen = iter(map(RFG._getGrid, fnin))

        self.fxnet = next(gridGen).merge(*gridGen)

        self.fxnet.collapse_policy = collapse_policy

        # success
        logger.verbose1('\nOriginal domain ' + str(self.fxnet))

        # start transformations
        logger.verbose2("\nTransforming...")

        if translate is not None:
            self.fxnet.translate(translate)
            logger.verbose2(f'\n...with translation {translate}')

        if domainSize:
            self.fxnet.setDomainSize( '(0,0,0)', domainSize )
            logger.verbose2(f'\n...with forced domain size: {self.fxnet!s}')

        self.nudgeTo=nudgeTo
        if nudgeTo > 0.0:
            self.fxnet.nudgeAll(nudgeTo)
            logger.verbose2(f'\n...with fracture nudging: {self.fxnet!s}')

        if forcedGridLines and sum(map(len,forcedGridLines))>0:
            # assume we have a triple of list-like things with numeric values
            for ia,gla in enumerate(forcedGridLines):
                for gl in gla:
                    self.fxnet.addGridline(ia,gl)

            logger.verbose2(f'\n...with forced grid lines: {self.fxnet!s}')

        if regGlSpacing and any(regGlSpacing):
            self.fxnet.addRegularGlSpacing(regGlSpacing)
            logger.verbose2(
                f'\n...with regular gridline spacing: {self.fxnet!s}')

        self.pmRefs = pmRefNearFx
        if pmRefNearFx:
            self.fxnet.refineNearFx(pmRefNearFx)
            logger.verbose2(f'\n...with PM refinement: {self.fxnet!s}')

        if maxGlSpacing and any(maxGlSpacing):
            self.fxnet.setMaxGlSpacing(maxGlSpacing)
            logger.verbose2(f'\n...with max gridline spacing: {self.fxnet!s}')

        logger.verbose1('\nTransformed domain ' + str(self.fxnet))

    @staticmethod
    def _getGrid(fn):
        return ofracs_module.parse(fn)

    def spewPreamble(self, moreMessages='', fout=sys.stdout):
        """Print a message about this grid and some of its parameters."""

        s = [ '!----- HGS code generated by {}'.format(__file__), ]
        s.append(f'! {datetime.datetime.now()}')
        s.append('! Dump of grid and explicit fractures defined in {}'.format(
                    self.fnin) )
        if self.nudgeTo:
            s.append(
                '! Fractures and gridlines "nudged" to {} increments.'.format(
                    self.nudgeTo) )
        if self.pmRefs:
            s.append(
                '! PM refined near fractures at increments: {}'.format(
                ", ".join( str(v) for v in self.pmRefs ) ) )

        if moreMessages:
            for l in moreMessages.splitlines():
                s.append( '! '+l )

        # assume that 'fout' is a file name string, a list of open files, or a
        # single open file
        if isinstance(fout,str):
            # a string, open the file and spew
            with open(fout,'w') as ftmp:
                print('\n'.join(l for l in s),end='\n\n', file=ftmp)
        elif hasattr(fout,'write'):
            # a file, write to it
            print('\n'.join(l for l in s),end='\n\n', file=fout)
        else:
            # assume a listlike, spew to each
            for flistitem in fout:
                print('\n'.join(l for l in s),end='\n\n', file=flistitem)

    def spewGrid(self, fout=sys.stdout):

        ps = partial(print, file=fout)

        if ( self.fxnet.isUniformGridSpacing(0) and
             self.fxnet.isUniformGridSpacing(1) and
             self.fxnet.isUniformGridSpacing(2) ):
            ps('generate uniform blocks')

            for a,(aStart,aEnd) in enumerate(self.fxnet.getBounds()):
                if aStart == 0.0:
                    ps('{:10.3f} {:5d}'.format(
                        aEnd, self.fxnet.getGridLineCounts()[a]-1))
                else:
                    ps('{:10.3f} {:5d} {:10.3f}'.format(
                        aEnd, self.fxnet.getGridLineCounts()[a]-1), aStart)

            ps('end grid generation')

        else:
            ps('generate variable blocks')
            for a in range(3):
                ps()
                printHGS(list(self.fxnet.iterGridLines(a)), fout)
            ps('end grid generation')


    def _pack_frac(self, f):
        # pack up the frac data, as it might be in one of two formats
        pack = None

        try:
            if f.__class__.__name__ == 'OFrac':
                pack = tuple( f.d ) + ( f.ap, )
            else:
                (pack, t) = f
        except Exception as e:
            raise ValueError(
                f'!!! Failed to pack fracture {f!s} of type {type(f)}') \
                from e

        return pack

    def spewFracs(self,
            fracsout=sys.stdout,
            fracpropasgt=sys.stdout,
            fpropsout=None,
            quantizeApertures=None,
        ):
        """Write out fracture definitions and "read properties".

        Args:
            fracsout : file-like
                The destination of the fracture definition text. Fracture
                definitions include the "choose faces" block and creating a new
                zone (or stating the aperture if no fpropsout file is specified).

            fracpropasgt : file-like
                The destination of "choose zone"/"read properties" text for
                aperture assignment. This can be the same as the preceding file
                object.

            fpropsout : file-like
                If not None...
                The destination of the fracture properties where apertures are
                assigned by fracture zone.

            quantizeApertures : int or None
                If not None...
                Reduce the precision in each fracture's aperture to be only
                'quantizeApertures' significant digits.
                This aims to result in fewer "read properties" calls in the
                fout because more fracture zones are grouped together with a
                single "read properties" statement.
        """

        # aliases for printing
        pfd = partial(print, file=fracsout) # fracture definition
        pfpa = partial(print, file=fracpropasgt) # fracture property asgt
        pprop = partial(print, file=fpropsout) # fracture properties

        prefixGuess = re.sub(r'o\.eco','',self.fnin)
        
        pfd( '\nuse domain type\n\tfracture\n\n'
            +'! begin explicit fractures\n')

        apQuant = defaultdict(list)
        apAssign = ''

        for i,f in enumerate(self.fxnet.iterFracs(), start=1):

            pack = self._pack_frac(f)

            s = '\n' + self._strFaceSelect(pack, withIndex=i)

            # never do this. It doesn't properly set the aperture
            #s += f'{self._strAperture(pack)}\n'

            # write a zone instead
            s += self._strNewZone(i)

            if quantizeApertures:
                # bin and write later
                apq = apQuantize(pack[-1],quantizeApertures)
                apQuant[apq].append(i)

            else:
                apAssign += '\nclear chosen zones\n'
                apAssign += self._strChooseZone(i)
                # write to buffer
                if fpropsout:
                    fracName=f'fracture_{i}'
                    apAssign += f'\nread properties\n\t{fracName}\n'
                    apStr = self._strAperture(pack[-1])
                    pprop(f'\n{fracName}\n{apStr}\nend')
                else:
                    apStr = self._strAperture(pack[-1])
                    s += f'\n{apStr}\n'
            pfd(s)

        pfd('! end explicit fractures\n\n')

        pfpa('! begin fracture property assignment\n')
        if fpropsout is not None:
            pfpa(f'properties file\n\t{fpropsout.name}\n\n')

        # print property assignments
        if quantizeApertures:
            for iq,(apq,fzonelist) in enumerate(apQuant.items()):
                s = '\nclear chosen zones\n'
                s += '\n'.join( self._strChooseZone(iz) for iz in fzonelist )

                apStr = self._strAperture(apq)

                if fpropsout:
                    apqName=f'fracture_ap_quantum_{iq}'
                    s += f'\nread properties\n\t{apqName}\n'
                    pprop(f'\n{apqName}\n{apStr}\nend')
                else:
                    s += f'\n{apStr}\n'

                pfpa(s)
        else:
            pfpa('\n\n'+apAssign)

        pfpa('!choose zones all\n!read properties\n!CommonFractureProperties')



    def _strFaceSelect(self, fracDataPack, withIndex=-1):
        """Make a string for face selection for a single fracture.

        Arguments:
            withIndex : int
                If this is >= 1, add to the string that adorns the "choose
                faces block" line.
        """

        (x1,x2,y1,y2,z1,z2,ap) = fracDataPack

        # adjust fracture coords to 3D block for face capture
        # Hope that the grid isn't spaced more finely than eps
        eps = Decimal('.001');
        x1 -= eps
        y1 -= eps
        z1 -= eps
        x2 += eps
        y2 += eps
        z2 += eps

        retStr = f'!---------- fracture{f" {withIndex}" if withIndex>0 else ""}\n'

        # build the text output
        retStr += 'clear chosen zones\n'\
            'clear chosen faces\n' \
            'choose faces block\n'

        fmt = '{: 7.3f}, {: 7.3f}\n'
        retStr += fmt.format(x1,x2)
        retStr += fmt.format(y1,y2)
        retStr += fmt.format(z1,z2)

        return retStr

    def _strNewZone(self, fracIndex):
        return f'new zone\n{fracIndex}'

    def _strChooseZone(self, fracZone):
        return f'choose zone number\n\t{fracZone}'

    def _strAperture(self, fracAperture):
        return f'aperture\n\t{D_AP(fracAperture)}'


    def pickleOFracGridToFile(self, fout):
        """dump to the open 'wb' file stream"""
        #import pdb ; pdb.set_trace()
        OFracGrid.pickleTo(self.fxnet, fout)

    def get_array_size_facts(self):
        """Return a string of array dimensioning facts"""

        ret = ''

        def product(l):
            return reduce((lambda x,y: x*y), l)

        ngl = self.fxnet.getGridLineCounts()
        nFxElements = self.fxnet.calcFxElementCount()
        nPMElements = product(map(lambda n: n-1,ngl))

        ret += '! minimum array_sizes.default values\n' \
            'fractures: 2d elements\n\t{}\nfractures: zones\n\t {}'.format(
                    nFxElements, self.fxnet.getFxCount())


        facts = [ ('mesh: x grid lines',ngl[0]),
                  ('mesh: y grid lines',ngl[1]),
                  ('mesh: z grid lines',ngl[2]),
                  ('pm: 3d brick elements',f'{nPMElements:,}'),
                  ('mesh: nnodes (excludes dual nodes)',f'{product(ngl):,}'),
                  ('mesh: total elements',f'{nPMElements+nFxElements:,}'),
        ]

        try:
            from scipy.stats import describe
            import numpy

            for a in 'xyz':
                spac = \
                    numpy.fromiter(iter(map(
                                    float, self.fxnet.iterGridLines(a))),float)

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    (nn,(mi,ma),amean,vari,skew,kurt) = describe(
                            spac[1:]-spac[:-1])

                facts.append(
                    (f'mesh spacing: {a} min, max, mean, stddev',
                     f'{mi:.3f}, {ma:.3f}, {amean:.3f}, {sqrt(vari):.3f}'))

        except ImportError:
            pass

        ret += '\n! other interesting facts ("!!" denotes non-HGS syntax)\n'
        ret += '\n'.join(list(f'! {l}\n!\t{v}' for l,v in facts[:3])) + '\n'
        ret += '\n'.join(list(f'!!{l}\n!!\t{v}' for l,v in facts[3:]))

        return ret

class OneLayerRFG(RFG):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._xz = self._make_xz_pairs()

        # hack the underlying OFracGrid to remove z-layers
        self.fxnet._gl[2] = [ 0., 1.,]

    def get_z(self,x):
        """Return the z of the given x"""
        ix = np.searchsorted(self._xz[:,0],x)
        return self._xz[ix,1]

    def _make_xz_pairs(self):
        """Finds xz pairs, assuming all fractures have zfrom=0.0"""
        x = np.fromiter(self.fxnet.iterGridLines(0),dtype=np.float32)
        #mnz = np.inf*np.ones(x.shape)
        mxz = np.zeros(x.shape)

        for f in self.fxnet.iterFracs():
            lbx = np.searchsorted(x,f.d[0])
            ubx = lbx + np.searchsorted(x[lbx:],f.d[1],side='right')

            zto = f.d[5]

            mxz[lbx:ubx] = np.fmax(zto,mxz[lbx:ubx])

        return np.hstack((x.reshape(-1,1),mxz.reshape(-1,1),))

    def _pack_frac(self, f):
        """Override parent method to use self.xz data"""

        pack = super()._pack_frac(f)

        lbx = np.searchsorted(self._xz[:,0], pack[0])
        ubx = lbx + np.searchsorted(self._xz[:,0], pack[1], side='right')

        pack = list(pack)
        pack[5] = D_CO(np.max(self._xz[lbx:ubx,1]))

        return tuple(pack)

    def spewGrid(self,fout=sys.stdout):
        """Print an xy- grid with 1-layer z- with undulating thickness"""

        # fall back to superclass grid output style if this grid does not
        # actually undulate
        if np.unique(self._xz[:,1]).size == 1:
            return super().spewGrid(fout)

        fw = fout.write # note, with write must add newlines at end
        ind = '  ' #indent
        il = 0 # indent level

        def _print_gls(ax):
            ngl = self.fxnet.getGridLineCounts()
            sep=' '
            fw(f'{ngl[ax]}\n')

            line = ''
            for i,g in enumerate(self.fxnet.iterGridLines(ax)):
                gstr = f'{sep}{g:8.3f}'
                if len(line)+len(gstr) >= 80:
                    fw(line+os.linesep)
                    line = gstr
                else:
                    line = line + gstr
            fw(line+os.linesep)

        fw('generate variable rectangles\n')
        _print_gls(0)
        _print_gls(1)

        il += 1
        fw(f'\n{il*ind}generate layers interactive\n')

        # base
        il += 1
        fw(f'{il*ind}base elevation\n')
        il += 1
        fw(f'{il*ind}elevation constant\n{(il+1)*ind}0.\n')
        il -= 1
        fw(f'{il*ind}end ! base elevation\n')

        # layer
        fw(f'\n{il*ind}new layer\n')
        il += 1

        fw(f'{il*ind}elevation from xz pairs\n')

        il += 1
        ilead = f'{il*ind}'
        for x,mxz in self._xz:
            fw(f'{ilead}{x:8.3f} {mxz:8.3f}\n')

        il -= 1
        fw(f'{il*ind}end ! elevation from xz pairs\n')

        fw(f'{il*ind}uniform sublayering\n{(il+1)*ind}1\n')
        il -= 1
        fw(f'{il*ind}end ! new layer\n')

        il -= 1
        fw(f'\n{il*ind}end ! generate layers interactive\n')
        fw('end ! generate variable rectangles\n')


    def get_array_size_facts(self):
        ret = super().get_array_size_facts()
        ret += f'\ngeneral: list\n\t{len(self._xz)}'
        return ret

def positive_decimal(v):
    """Return v as a `decimal.Decimal`"""
    try:
        d = Decimal(v)
    except decimal.InvalidOperation as e:
        raise ValueError('Could not convert "{v}" to Decimal')

    if d <= 0.0:
        raise ValueError(f'{v} is not > 0')
    return d

class _Triple_List(argparse.Action):
    """for parsing --reg-grid-space and --max-grid-space"""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def __call__(self, parser, namespace, values, option_string=''):
        ell = getattr(namespace,re.sub('-','_',option_string[2:]))
        for v in values:
            self._parse(v, ell, parser, option_string)
        return ell

    @staticmethod
    def _parse(v,ell, parser, opts):
        """Add v to the triple_list ell"""

        def _as_Decimal(val, ell=ell):
            d = Decimal(val)
            if any(ell):
                raise RuntimeError('Some values already assigned')
            for i in range(len(ell)):
                ell[i] = d

        def _as_ax_eq_val(s, ell=ell):
            """Find one or more axid=value string in `s`"""
            for m in re.findall(r'([xyz])\s*=\s*('+_num_re_str+')',s):
                ax, val = m[0], Decimal(m[1])
                iax = 'xyz'.index(ax.lower())
                if ell[iax] is not None:
                    parser.error(f'{ax} value assigned multiple times in {opts}')
                ell[iax] = val

        pl = [ _as_Decimal, _as_ax_eq_val, ]

        errret = []
        for ip,p in enumerate(pl):
            try:
                p(v)
            except Exception as e:
                errret.append(f'{p.__name__} gave {e!r}')
            else:
                break

        if len(errret) == len(pl):
            #raise RuntimeError(f'when parsing "{v}"\n'+'\n'.join(errret))
            parser.error(f'Invalid argument {f"{opts}=" if opts else ""}{v}')

        return ell

class _OneOrTwoFileNames(argparse.Action):
    '''For parsing --frac-out'''
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def __call__(self, parser, namespace, values, option_string=''):
        res = [sys.stdout, sys.stdout]
        if len(values) == 1:
            res[0] = open(values[0], 'w')
            res[1] = res[0]
        elif len(values) == 2:
            res[0] = open(values[0], 'w')
            res[1] = open(values[1], 'w')
        else:
            parser.error(f'{option_string} expecting one or two filenames.')

        setattr(namespace, self.dest, res)

def make_arg_parser():
    # set up command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument( '-v', '--verbosity',
            default=0,
            action='count',
            help='''Increase the verbosity with each "-v". Default is verbosity
            level 0, which only produces output related to node and gridline
            counts in the style of array_sizes.default; Higher levels give
            information on the operations of this script.''',
            )

    parser.add_argument(
            '-n', '--nudge-fx-coordinates-to',
            type=Decimal,
            default=Decimal('0'),
            help='Nudge each fracture to the nearest interval' )

    parser.add_argument( '-r', '--refine-near-fx-plane',
            nargs='+',
            type=positive_decimal,
            action='extend',
            help="""Add fracture-parallel gridline(s) in PM at the specified
            interval (or sequence of intervals) away from fracture planes.""")

    parser.add_argument(
        '--domain-size', default='', type=str,
        help='''Set the domain extent, as "X, Y, Z" or "(X,Y,Z)", and
        truncate the domain such that only the region between (0,0,0) and
        (X,Y,Z) remains.''',
    )

    parser.add_argument(
        '--translate', default=None, type=str,
        help='''Translate the domain by the given distances, as "Xt, Yt, Zt" or
        "(Xt,Yt,Zt)", by adding the given values to each coordinate in the
        domain/fractre network. Translation is performed before domain size
        truncation''',
        )


    parser.add_argument(
        '--force-extra-grid-lines',
        #nargs=1, # force list
        action='append',
        default=[],
        help='''Add in extra grid lines specified here, e.g., "x= .1, 250,
        399.9 ; z=0.1,199.9" If an axis is not specified, add this grid line to
        each dimension.''')

    parser.add_argument(
        '--reg-grid-space',
        nargs='+',
        metavar=('V | x=V1', 'y=V2'),
        action=_Triple_List,
        default=[None,None,None,],
        help='''Add regular spacing at the specified intervals from the origin.
        E.g. '5' for a regular grid spacing of 5 in all directions; 'x=5 z=1'
        to apply regular spacing to the x- and z- axes only'''
        )

    parser.add_argument(
        '--max-grid-space',
        nargs='+',
        metavar=('V | x=V1', 'y=V2'),
        action=_Triple_List,
        default=[None,None,None,],
        help='''The maximum space between gridlines. See format of
        --reg-grid-space'''
        )

    parser.add_argument(
            '--grid-out',
            default=sys.stdout,
            type=argparse.FileType('w'),
            help='Grid output file name')

    parser.add_argument( '--frac-out',
            nargs='+',
            #metavar='FRAC_DEFS_OUT [, FRAC_PROP_ASGT_OUT]',
            metavar=('FRAC_DEFS_OUT', 'FRAC_PROP_ASGT_OUT'),
            default=(sys.stdout, sys.stdout),
            action=_OneOrTwoFileNames,
            help='''Fracture definitions (and property assignment) output file
            name. Optionally, supply a second file name to store the property
            assignments. Default behaviour is to write to stdout, which is
            equivalent to '--frac-out - -'. This is useful if multiple .fprops
            files are needed, like an initial 'read properties generic.fprops'
            directive to apply generic properties to all fractures, then a
            second 'read properties FPROPS_OUT' file to assign/override
            apertures, etc.''',
            )

    parser.add_argument( '--fprops-out',
            default=None,
            type=argparse.FileType('w'),
            help='Fracture (property definitions) output file name')

    parser.add_argument( '--pickle-out',
            default=False,
            #type=argparse.FileType('wb'),
            type=str,
            help='Output file name for OFracGrid-type object data (pickled)')

    parser.add_argument( '--fx-collapse-policy',
            choices=getattr(ofracs_module,'__FX_COLLAPSE_POLICIES__'),
            default='warn-omit',
            help=f'''Define what to do if a fracture collapses when nudging,
            etc. Default "warn-omit".''',
            )

    parser.add_argument( '--quantize-apertures',
            default=None,
            metavar='S',
            type=int,
            help='''Round fracture apertures to at most S significant figures
            and assign all fracutres within that "aperture quantum" to the same
            .fprops material type (where material types only define apertures).
            This reduces the number of material zones and hence the number of
            .fprops file read-throughs that GROK will do. e.g., 20,000 fractures
            has been shown to take 1 hour to process!
            ''')

    return parser

if __name__ == "__main__":

    parser = make_arg_parser()

    parser.add_argument(
            'filename',
            metavar='FILE',
            help='filename or prefix of orthogonal fracture network source files',
            nargs='+')

    args = parser.parse_args()

    logger.setLevel(logging.INFO+1 - args.verbosity)

    setattr(ofracs_module,
            '__VERBOSITY__',
            max(0,args.verbosity-1) )
    setattr(ofracs_module,
            '__FX_COLLAPSE_POLICY__',
            args.fx_collapse_policy )

    # for parsing forced gridlines
    fglvals = {'x':set(), 'y':set(), 'z':set()}

    for fglarg in args.force_extra_grid_lines:
        for m in _gl_list_re.finditer(fglarg):

            vals = list(float(x) for x in _num_re.findall(m.group('vals')))

            if m.group('axis'):
                fglvals[ m.group('axis').lower() ].update(vals)
            else:
                for axSet in fglvals.values():
                    axSet.update(vals)

    def _to_decimal_list(s):
        if s in (None, '', False):
            return None
        else:
            return list(
                Decimal(v) for v in
                re.sub(r'[(),]',' ', s).split())

    # get parse grid and make modifications
    rfg = RFG(
        args.filename,
        args.fx_collapse_policy,
        domainSize=_to_decimal_list(args.domain_size),
        translate=_to_decimal_list(args.translate),
        forcedGridLines=list( fglvals.values() ),
        nudgeTo=args.nudge_fx_coordinates_to,
        maxGlSpacing=args.max_grid_space,
        regGlSpacing=args.reg_grid_space,
        pmRefNearFx=args.refine_near_fx_plane,
    )


    def quoteMultipartParam(s):
        if any(specialChar in s for specialChar in '$(), '):
            s = '"'+s+'"'
        return s

    otherMessages = """Whole Argument string {}\nat {}""".format(
                " ".join(map(quoteMultipartParam, sys.argv)),
                datetime.datetime.now())

    rfg.spewPreamble(
            moreMessages=otherMessages,
            fout=[args.grid_out, args.frac_out[0], args.fprops_out])
    if args.frac_out[0] != args.frac_out[1]:
        rfg.spewPreamble( moreMessages=otherMessages, fout=args.frac_out[1])

    rfg.spewGrid(fout=args.grid_out)
    rfg.spewFracs(fracsout=args.frac_out[0],
            fracpropasgt=args.frac_out[1],
            fpropsout=args.fprops_out,
            quantizeApertures=args.quantize_apertures,
        )


    if args.pickle_out:
        rfg.pickleOFracGridToFile(args.pickle_out)


    print(rfg.get_array_size_facts())

