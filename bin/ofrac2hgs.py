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
from math import sqrt
from functools import reduce
from collections import defaultdict

import numpy as np


from ofracs import (D_AP, OFrac)
from pyhgs.gridgen_tools import *
from pyhgs.parser.rfgen import *
from parser_fractran import *
from parser_hgs_rfgeneco import *

__VERBOSITY__ = 0


WISHLIST="""

Have RFG inherit from the OFracGrid object (or don't have an RFG class at all)

"""

# handy REs
_numberREStr = r'[+-]?(?:\d*\.)?\d+(?:[dDeE][+-]?\d+)?'
_numRE = re.compile(_numberREStr)
_gl_val_RE = re.compile(
    r';?\s*(?:(?P<axis>[xyzXYZ])\s*=\s*)?(?P<vals>(?:,?\s*'
            +_numberREStr+'){1,})')

def apQuantize(v,n_sig):
    """Determine the aperture to n_sig significant figures."""

    if Decimal(v) < Decimal('0.0000005'):
        warnings.warn(f'Encountered aperture < 1um: {v}')
        v = 1e-6
    #if Decimal(v) < Decimal('1e-6') or v > 1.:
    #    raise ValueError(f'v must be between zero and one, exclusive')

    least_sig = max(1e-6,10**(floor(log10(v))-n_sig+1))

    return Decimal(v).quantize(Decimal(f'{least_sig:.6f}'.strip('0')))

class RFG:
    """Parse/store/manipulate/produce HGS-style grids and fracture definitions

    TODO: refactor and use more code from 'ofracs' module
    """

    def __init__(self,
            fnin,
            domainSize=[],
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

        # success
        if __VERBOSITY__:
            print( '\nOriginal domain ' + str(self.fxnet), file=sys.stderr )
            if __VERBOSITY__ > 1:
                print("\nTransforming...")

        # start transformations
        if domainSize:
            self.fxnet.setDomainSize( '(0,0,0)', domainSize )

            if __VERBOSITY__ > 1:
                print(f'\n...with forced domain size: {self.fxnet!s}',
                        file=sys.stderr)

        self.nudgeTo=nudgeTo
        if nudgeTo > 0.0:
            self.fxnet.nudgeAll(nudgeTo)
            if __VERBOSITY__ > 1:
                print(f'\n...with fracture nudging: {self.fxnet!s}',
                        file=sys.stderr )

        if forcedGridLines and sum(map(len,forcedGridLines))>0:
            # assume we have a triple of list-like things with numeric values
            for ia,gla in enumerate(forcedGridLines):
                for gl in gla:
                    self.fxnet.addGridline(ia,gl)

            if __VERBOSITY__ > 1:
                print(f'\n...with forced grid lines: {self.fxnet!s}',
                        file=sys.stderr )

        if regGlSpacing and any(regGlSpacing):
            self.fxnet.addRegularGlSpacing(regGlSpacing)
            if __VERBOSITY__ > 1:
                print(f'\n...with regular gridline spacing: {self.fxnet!s}',
                        file=sys.stderr )

        self.pmRefs = pmRefNearFx
        if pmRefNearFx:
            self.fxnet.refineNearFx(pmRefNearFx)
            if __VERBOSITY__ > 1:
                print(f'\n...with PM refinement: {self.fxnet!s}',
                        file=sys.stderr )

        if maxGlSpacing and any(maxGlSpacing):
            self.fxnet.setMaxGlSpacing(maxGlSpacing)
            if __VERBOSITY__ > 1:
                print(f'\n...with max gridline spacing: {self.fxnet!s}',
                        file=sys.stderr )

        if __VERBOSITY__:
            print('\nTransformed domain ' + str(self.fxnet), file=sys.stderr)

    @staticmethod
    def _findParser(fn):
        errmsg = ''
        retParser = None

        # create a prioritized list of parser types
        parsers = [
             OFracGrid.PickleParser,
             RFGenOutFileParser,
             HGSEcoFileParser,
            ] \
            + list(iterFractranParsers())

        # try some different parsers
        for ParserClass in parsers:
            try:
                retParser = ParserClass(fn)
            except (
                  NotValidFractranInputError,
                  NotValidOFracGridError,
                  NotValidRFGenInputError,
                  NotValidHGSEcoInputError
                   ) as e:
                errmsg += '\n  '+ParserClass.__name__+\
                          ' did not work- {}'.format(str(e))
                retParser = None
            if retParser:
                return retParser

        if not retParser:
            print('Could not parse input file "{}":{}'.format(fn, str(errmsg)),
                    file=sys.stderr)
            sys.exit(1)

    @staticmethod
    def _getGrid(fn):
        return RFG._findParser( fn ).getOFracGrid()

    def spewPreamble(self, moreMessages='', fout=sys.stdout):
        """Print a message about this grid and some of its parameters."""

        s = [ '!----- HGS code generated by {}'.format(__file__), ]
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

        stdoutSave = sys.stdout
        sys.stdout = fout

        if ( self.fxnet.isUniformGridSpacing(0) and
             self.fxnet.isUniformGridSpacing(1) and
             self.fxnet.isUniformGridSpacing(2) ):
            print('generate uniform blocks')

            for a,(aStart,aEnd) in enumerate(self.fxnet.getBounds()):
                if aStart == 0.0:
                    print('{:10.3f} {:5d}'.format(
                        aEnd, self.fxnet.getGridLineCounts()[a]-1))
                else:
                    print('{:10.3f} {:5d} {:10.3f}'.format(
                        aEnd, self.fxnet.getGridLineCounts()[a]-1), aStart)

            print('end grid generation')

        else:
            print('generate variable blocks')
            for a in range(3):
                print()
                printHGS(list(self.fxnet.iterGridLines(a)))
            print('end grid generation')

        sys.stdout = stdoutSave

    def _pack_frac(self, f):
        # pack up the frac data, as it might be in one of two formats
        pack = None
        if type(f) == OFrac:
            pack = tuple( f.d ) + ( f.ap, )
        else:
            (pack, t) = f
        return pack

    def spewFracs(self,
            fout=sys.stdout,
            fpropsout=None,
            quantizeApertures=None,
        ):
        """Write out fracture definitions and properties.

        Args:
            fout : file-like
                The destination of the fracture definition text. Fracture
                definitions include the "choose faces" block and creating a new
                zone or stating the aperture.

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
        prefixGuess = re.sub('o\.eco','',self.fnin)
        more = '\n!use domain type\n!fracture\n\n!properties file\n'
        more +='!{}.fprops\n\n'.format( prefixGuess )
        more +='! begin explicit fractures\n'
        print( more, file=fout )

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
                    print(f'\n{fracName}\n{apStr}\nend',file=fpropsout)
                else:
                    s += f'{apStr}\n'
            print(s, file=fout)


        # print property assignments
        if quantizeApertures:
            for iq,(apq,fzonelist) in enumerate(apQuant.items()):
                s = '\nclear chosen zones\n'
                s += '\n'.join( self._strChooseZone(iz) for iz in fzonelist )

                apStr = self._strAperture(apq)

                if fpropsout:
                    apqName=f'fracture_ap_quantum_{iq}'
                    s += f'\nread properties\n\t{apqName}\n'
                    print(f'\n{apqName}\n{apStr}\nend',file=fpropsout)
                else:
                    s += f'{apStr}\n'

                print(s,file=fout)
        else:
            print('\n\n'+apAssign,file=fout)


        more ='! end explicit fractures\n\n'
        more += '!choose zones all\n!read properties\n!CommonFractureProperties'
        print( more, file=fout )



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
    """for parsing --reg-gl-space and --max-gl-space"""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def __call__(self, parser, namespace, values, option_string=''):
        ell = getattr(namespace,re.sub('-','_',option_string[2:]))
        for v in values:
            self._parse(v,ell)
        return ell

    @staticmethod
    def _parse(v,ell):
        """Add v to the triple_list ell"""

        def _as_Decimal(val, ell=ell):
            d = Decimal(val)
            if any(ell):
                raise RuntimeError('Some values already assigned')
            ell = 3*[d,]

        def _as_ax_eq_val(s, ell=ell):
            m=re.match(r'([xyz])\s*=\s*('+_numberREStr+')',s)
            ax, val = (m.group(1),Decimal(m.group(2)),)
            iax = 'xyz'.index(ax.lower())
            if ell[iax] is not None:
                raise RuntimeError(f'{ax} value already assigned')
            ell[iax] = val

        pl = [ _as_Decimal, _as_ax_eq_val, ]

        errret = []
        for ip,p in enumerate(pl):
            try:
                p(v)
            except Exception as e:
                errret.append(f'{p}: {e!s}')
            else:
                break

        if len(errret) == len(pl):
            raise RuntimeError('\n'.join(errret))

        return ell

def make_arg_parser():
    # set up command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument( '-v', '--verbosity',
            default=0,
            action='count',
            help='increase the verbosity with each "-v"' )

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
        help='domain dimensions, as "dx, dy, dz" or "(dx,dy,dz)' )

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
        action=_Triple_List,
        default=[None,None,None,],
        help='''Add regular spacing at the specified intervals from the origin.
        E.g. '5' for a regular grid spacing of 5 in all directions; 'x=5 z=1' to
        apply regular spacing to the x- and z- axes only'''
        )

    parser.add_argument(
        '--max-grid-space',
        nargs='+',
        action=_Triple_List,
        default=[None,None,None,],
        help='The maximum space between gridlines. E.g. "5", "x=5 z=1", ...')

    parser.add_argument(
            '--grid-out',
            default=sys.stdout,
            type=argparse.FileType('w'),
            help='Grid output file name')

    parser.add_argument( '--frac-out',
            default=sys.stdout,
            type=argparse.FileType('w'),
            help='Fracture (spatial definitions) output file name')

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
            choices=getattr(sys.modules['ofracs'],'__FX_COLLAPSE_POLICIES__'),
            default=getattr(sys.modules['ofracs'],'__FX_COLLAPSE_POLICY__'),
            help=f'''Define what to do if a fracture collapses when nudging, etc
            Default {getattr(sys.modules['ofracs'],'__FX_COLLAPSE_POLICY__')}'''
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

    #import pdb ; pdb.set_trace()

    __VERBOSITY__ = args.verbosity

    setattr(
            sys.modules['ofracs'],
            '__VERBOSITY__',
            max(0,__VERBOSITY__-1) )
    setattr(
            sys.modules['ofracs'],
            '__FX_COLLAPSE_POLICY__',
            args.fx_collapse_policy )

    # for parsing forced gridlines
    fglvals = {'x':set(), 'y':set(), 'z':set()}

    for fglarg in args.force_extra_grid_lines:
        for m in _gl_val_RE.finditer(fglarg):

            vals = list(float(x) for x in _numRE.findall(m.group('vals')))

            if m.group('axis'):
                fglvals[ m.group('axis').lower() ].update(vals)
            else:
                for axSet in fglvals.values():
                    axSet.update(vals)

    # get parse grid and make modifications
    rfg = RFG(
            args.filename,
            domainSize=list(
                Decimal(v) for v in
                re.sub(r'[(),]',' ',args.domain_size).split() ),
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
            fout=[args.grid_out, args.frac_out, args.fprops_out])

    rfg.spewGrid(fout=args.grid_out)
    rfg.spewFracs(fout=args.frac_out,
            fpropsout=args.fprops_out,
            quantizeApertures=args.quantize_apertures,
        )


    if args.pickle_out:
        rfg.pickleOFracGridToFile(args.pickle_out)


    print(rfg.get_array_size_facts())
