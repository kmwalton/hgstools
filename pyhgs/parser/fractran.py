#!/usr/bin/env python3
"""A parser for Fractran fracture networks and grids

Parses 
- *o.lst
- *o.xyc
- *o.fracture_lengths
- *o.hfz and *o.vfz
- *o.inc
- *o.CNN
- *o.HNN


Manipulates parsed data to form an OFracGrid object.

Access different parser types by the iterFractranParsers() method. Just pass a
filename or prefix to the returned objects -- those objects will auto-detect the
required files.
"""

VERBOSE = False

import sys,os,io,itertools,re
from scipy.io import FortranFile
import numpy as np
from math import log10

from ofracs import *

class NotValidFractranInputError(Exception):
   """Exception raised if input file is not valid"""
   pass


class Parser_FractranLengthsFile():
   """Parse a Fractran o.fracture_lengths file"""

   def __init__(self, fnin):
      """Read-in the datafile"""
      # check if just prefix
      if not fnin.endswith('o.fracture_lengths'):
          fnin+='o.fracture_lengths'

      #data file given on command line
      self._fnin = fnin

      # get lists vertical and horizontal fractures
      ( verticalFracs, horizontalFracs ) = (None,None)

      try:
          with open( fnin, 'r' ) as fin:
              ( verticalFracs, horizontalFracs ) = self._readFracLists(fin)
      except FileNotFoundError:
          raise NotValidFractranInputError(f'Could not open {fnin}.')


      # merge the fracture lists
      fracsHere = []

      for ( zfrom, zto, length, xpos, aperture ) in verticalFracs:
         fracsHere.append( ( xpos, xpos, 0.,1., zfrom, zto, aperture ) )

      for ( xfrom, xto, length, zpos, aperture ) in horizontalFracs:
         fracsHere.append( ( xfrom, xto, 0.,1., zpos, zpos, aperture ) )

      self.nFracsInEachOrientation = (len(verticalFracs),len(horizontalFracs))
      self.fxgrid = OFracGrid(fx=fracsHere)

   def getOFracGrid(self):
      """return a deep copy of the fracture network"""
      return copy.deepcopy( self.fxgrid )

   def getFracCounts(self):
      return self.nFracsInEachOrientation

   def _readFracLists(self, fin):
      """Read successive lists of fractures, put fractures in 'bins' based on header line text"""

      # REs to detect the following header lines
      # i zbot    ztop   length   xpos  aperture
      # i xleft   xright length   zpos  aperture
      vfHeader = re.compile(
         '^{}$'.format('\s+'.join(
             ['i','zbot','ztop','length','xpos','aperture'])))

      hfHeader = re.compile(
         '^{}$'.format('\s+'.join(
             ['i','xleft','xright','length','zpos','aperture'])))

      atEOF = False
      foundAnyHeader = False
      vfList = list()
      hfList = list()

      # look for header+lists until the end of the file
      while not atEOF:
          foundVFHdr = False
          foundHFHdr = False

          # find a header line, or fail
          line = fin.readline()
          while True:
             if line == '': atEOF = True; break
             if vfHeader.match(line.strip()): foundVFHdr = True; break
             if hfHeader.match(line.strip()): foundHFHdr = True; break
             line = fin.readline()

          if foundVFHdr:
             foundAnyHeader = True
             # make list objects to exhaust the generators made in _readFracList
             vfList.extend( list(self._readFracList(fin)) )
          elif foundHFHdr:
             foundAnyHeader = True
             hfList.extend( list(self._readFracList(fin)) )

      # fail if we found no headers
      if not foundAnyHeader:
         raise NotValidFractranInputError("Could not find any required fracture list header line in "+self._fnin)
      # fail if we found no fractures
      if not vfList and not hfList:
         raise NotValidFractranInputError("Could not find any fractures"+self._fnin)

      return (vfList, hfList)



   def _readFracList(self, fin):
      """Read and yield a list of fractures (as tuples of float values)"""

      reInt='[+-]?[0-9]+' # RE for a signed integer
      reFix='[+-]?[0-9]*\.?[0-9]+' # RE for a fixed point number

      fracEntry= \
         re.compile('^{0}\s+({1})\s+({1})\s+({1})\s+({1})\s+({1})\s*$'.format(reInt,reFix))

      lastPos = fin.tell()

      while True:
         l = fin.readline()

         # end of file
         if l == '':
            raise StopIteration


         m = fracEntry.match(l.strip())
         if m:
            yield tuple( float(v) for v in m.groups() )
         else:
            fin.seek(lastPos)
            raise StopIteration

         lastPos = fin.tell()


class DotLST_File(object):
    """Parse the .lst file and retrieve some useful inputs/outputs
        
    e.g.
    
         foo = DotLST_File("somePrefix")
         phi = foo['porosity']
    or
         foo = DotLST_File("somePrefix")
         (nGridLine_x,nGridLine_z) = foo.getShape()
         KxxZone2 = foo.Kxx[2]

    Currently implemented keys in this dictionary can be listed by
    DotLST_File.printMoreHelp()
"""

    
    PARAM_SEARCHES = [
        ('g',  'int',  'nx',          "nodes in x-direction"),
        ('g',  'int',  'nz',          "nodes in z-direction"),
        ('zl', 'float','porosity',    "effective porosity"),
        ('zl', 'float','Kxx',         "hydraulic conductivity in x-direction"),
        ('zl', 'float','Kzz',         "hydraulic conductivity in z-direction"),
        ('zl', 'float','bulk_density',"soil bulk density"),
        ('zl', 'float','retardation', r"^\s*retardation factor(?! of)"),
        ('g',  'float','g',           'gravity constant'),
        ('g',  'float','grav_accel',  'gravity constant'),
        ('g',  'float','fluid_viscosity','fluid viscosity'),
        ('g',  'float','fluid_density','fluid density'),
        ('g',  'float','porosity_fx','fracture porosity'),
        ('g',  'float','type1_inflow',re.compile(r'inflow at 1st-type nodes:(.+)',flags=re.I)),
        ('g',  'float','type1_outflow',re.compile(r'outflow at 1st-type nodes:(.+)',flags=re.I)),
        ('g',  'float','type2_inflow',re.compile('inflow at 2nd-type nodes:(.+)',flags=re.I)),
        ('g',  'float','type2_outflow',re.compile('outflow at 2nd-type nodes:(.+)',flags=re.I)),
        ]
    """Tuples of quantity (zl=zone list; g=global)
    the data type (int, float)
    the parameter name for use in retrieval later
    the regular expression string used in the search."""
 
    @staticmethod
    def printMoreHelp():
        """Print the list of keywords available"""
        print('\n'.join(pName for (q,t,pName,stuff) in DotLST_File.PARAM_SEARCHES))


    def __init__(self,prefix=None):
        """Read the file and put useful data into the dictionary"""

        lines=[]

        if not prefix:
            # look in current directory for a single .pre file and infer the prefix
            # Use the first one (hope there's only one)
            # Grab all but the extension.
            prefix=glob.glob('*.pre')[0][0:-4]

        # Read some of the output file into the lines array.
        with open( prefix+'o.lst', 'r') as fin:
            lines = fin.readlines() # get the whole thing

        # store the prefix
        self.prefix = prefix
        self.directoryName = os.path.split(
                os.path.split(
                    os.path.abspath(prefix+'o.lst'))[0])[1]

        lines = list(map(lambda l:l.lower().strip(), lines))

        self.d = {}

        #import pdb ; pdb.set_trace()

        for (qty,typ,pName,srchLine) in DotLST_File.PARAM_SEARCHES:

            valueIter = None
            if type(srchLine) == str:
                # assume space delimits search line and the required value. get last word
                pat = re.compile(srchLine)
                valueIter = map(
                    lambda l: l.split()[-1],
                    filter(pat.search, lines))

            elif type(srchLine) == type(re.compile("dummy")):
                #import pdb ; pdb.set_trace()
                valueIter = map(
                    lambda m: m.group(1),
                    filter(bool, map(srchLine.search, lines)))
            else:
                raise RuntimeError("Invalid PARAM_SEARCHES item found.")

            if qty == 'zl':
                # 'list' stores the sequence of values from the matching lines
                # 'eval(typ)' does the appropriate type cast

                # split grabs the right-most, space separated "word" from teh
                # matching line
                self.d[pName] = list(map(eval(typ), valueIter))

            elif qty == 'g':
                # 'eval(typ)' does the appropriate type cast
                # 'next' looks for only the first line, l, matching the search
                # split grabs the right-most, space separated "word" from teh
                # matching line
                self.d[pName] = eval(typ)(next(valueIter))

            else:
                raise Exception('Unexpected "Quantity" in PARAM_SEARCHES')


        self.d['units']='?'

        units = {
            'L':[ (.001,'mm'),
                  (.01 ,'cm'),
                  (1.0 ,'m'),
                  (1000.0 ,'km'), ],
            'T':[ (1.0, 's'),
                  (1.0/86400, 'd'),
                  (1.0/86400/365, 'yr'), ],
            'M':[ (1.0 ,'kg'),
                  (0.001,'g'),
                  (0.000001,'mg'), ],
        }

        # try to guess the input units
        bench_g = 9.80665 # m/s^2
        bench_mu = 0.001308 # kg/m/s @ 10C 
        gUnit = []
        muUnit = []

        for il,(fl,sl) in enumerate(units['L']):
            for it,(ft,st) in enumerate(units['T']):
                test_g = self.d['grav_accel']/fl*pow(ft,2)
                gUnit.append(
                    (abs(log10(test_g/bench_g)), il, it) )


        (logdiff, bestL, bestT) = sorted(gUnit)[0]
                
        if logdiff > 0.05: # assume g varies by less than 5%
            print("Warning: could not guess the units for L and T",file=sys.stderr)
        
        for im,(fm,sm) in enumerate(units['M']):
            test_mu = self.d['fluid_viscosity'] / fm * units['L'][bestL][0] * units['T'][bestT][0]
            muUnit.append( (abs(log10(test_mu/bench_mu)), im) )

        (logdiff, bestM) = sorted(muUnit)[0]
        if logdiff > 0.5: # viscosity can vary by .5 OoM in 10-40 deg C
            print("Warning: could not guess the units for M",file=sys.stderr)

        self.d['units'] = '{} {} {}'.format(
                    units["M"][bestM][1],
                    units["L"][bestL][1],
                    units["T"][bestT][1] )
        self.d['unit'] = { 'M':units["M"][bestM][1],
                           'L':units["L"][bestL][1],
                           'T':units["T"][bestT][1] }


        #import pdb ; pdb.set_trace()
        if VERBOSE:
            print(f'Found units: {self.d["units"]}')
                


    def __getattr__(self,key):
        """For access to allParamsObj.d['param'] as allParmsObj.param"""
        if key in self.d:
            return self.d[key]
        return super().__getattribute__(key)
        
    def __getitem__(self,key):
        """For access to allParamsObj.d['param'] as allParmsObj['param']"""
        return self.d[key]

    def __str__(self):
        """Print out the parameter names and values extracted from the
            underlying data file"""

        s=[]

        KW=max(map(lambda t: len(t[2]),DotLST_File.PARAM_SEARCHES))
        SEP=' > '

        for (q,t,key,junk) in DotLST_File.PARAM_SEARCHES:
            if q=='zl':
                vals=", ".join(map(str, self.d[key]))
                s.append(f'{key:{KW}}{SEP}{vals}')
            elif q=='g':
                s.append(f'{key:{KW}}{SEP}{self.d[key]}')


        return '\n'.join(s)

    def getShape(self):
        """Return the number of gridlines in the domain as a tuple"""
        return (self.nx,self.nz)

class DotFZ_File:
    """Parse and store fracture index->(hzone, zone, node0, node1, aperture)"""

    def __init__(self, fn):
        self.orientation = fn[-3]
        self.fn = fn
        try:
            with FortranFile( fn, 'r' ) as fin:
                self.n = fin.read_ints()[0]
                iHZone = fin.read_ints()
                iZone = fin.read_ints()
                iNF = fin.read_ints().reshape(self.n,2)
                rAps = fin.read_reals()
        except:
            raise RuntimeError("Couldn't open or parse the FortranFile "+fn)

        # allocate
        self.d = len(rAps) * [ (-1,-1,-1,-1,0.0,), ]

        # subtract 1 from the index specified in this file to go from Fortran's
        # 1-base index to python 0-based
        for i,v in enumerate(zip( iHZone, iZone, iNF, rAps )):
            self.d[i] = (v[0],v[1],v[2][0]-1,v[2][1]-1,v[3],)

    def __getitem__(self, i):
        """return a tuple (H-Zone, Zone, n0, n1, aperture) for fracture element i
            
        Node indicies have been corrected from 1-based to 0-based
        """
        return self.d[i]

    def __len__(self):
        return len(self.d)

    def __str__(self):

        s = '{} {} elements in {}'.format(
            str(self.n), self.orientation, self.fn)

        #for arr in [ self.iHZone, self.iZone, self.iNF, self.rAps, ] :
        #    s += '\n{}, n={}'.format(str(arr),len(arr))

        return s

class DotCNN_File:
    """Parse and store data from 'o.cNN' files"""

    def __init__(self, fn, shape=None):
        """Parse the o.cNN or o.h01 file. Set the shape of the data to (nx, nz) if
            specified.
        """
        if not re.search(r'o\.[ch]\d\d$', fn):
            fn += 'o.c01'

        self.fn = fn
        fin = FortranFile( fn, 'r' )

        # title strings should look like either:
        # Heads at steady-state
        # Concentration at 10000.0
        self.title = fin.read_record([('bText','50a')])['bText'][0].decode('ascii').strip()
        self.time = self.title.split()[-1]

        (self.nx, self.nz) = (None,None)
        if shape:
            (self.nx, self.nz) = shape
            self.cNodal = fin.read_reals().reshape(self.nx, self.nz)
        else:
            self.cNodal = fin.read_reals()

        fin.close()

    def get(self, ix, iz):
        """Assume that shape is set. If not, this will fail"""
        return self.cNodal.item((iz,ix))

    def __str__(self):
        return "{}, # nodal data points={}".format(self.title, self.cNodal.size)

    def prettyPrint(self):
        return str(self) + '\n' + str(self.cNodal)

class DotHNN_File(DotCNN_File):
    """Simple alias to DotCNN_File, as these have the same format"""

    def __init__(self, fn, shape=None):
        """Parse the o.h01 file. Set the shape of the data to (nx, nz) if
            specified.
        """

        if not re.match(r'\S+o\.h\d\d', fn):
            fn += 'o.h01'

        super().__init__(fn, shape)


#class DotAPS_File():
#    """This is an ascii file."""
#    def __init__(self, fn):
#        self.fn=fn
#
#        fin = FortranFile( fn, 'r' )
#
#        self.v = fin.read_ints()[0]
#        self.ra = fin.read_reals()[0]
#
#    def __str__(self):
#        return self.fn

class DotXYC_File:
    """Parse&store node index->(x,z) + useful data extractions"""
    def __init__(self, fn, shape=None):
        self.fn = fn
        with FortranFile(fn, 'r') as fin:
            self.n = fin.read_ints()[0]
            self.d = fin.read_reals().reshape(self.n,2)

        self.s = shape
        self.gl = None


    def __iter__(self):
        for i in range(self.n):
            yield (self.d[i][0], self.d[i,1], )

    def __getitem__(self, i):
        """return a coordinate (x, z) of node i"""
        return self.d[i]

    def __str__(self):
        return '{} has {} node (x,z) pairs'.format(self.fn, self.n)


    def setShape(self, shape):
        self.s = shape

    def getGridLines(self):
        """Get the ([x grid lines],[z grid lines]) in this object"""
        # store
        if not self.gl:
            self.gl = (self.d[::self.s[1],0],self.d[0:self.s[1],1])
        return self.gl

    def getElemCentroidGridLines(self):
        """Get the centroid ([x lines],[z lines]) in this object"""
        (xgl,zgl) = self.getGridLines()
        return ((xgl[1:]+xgl[:-1])/2.0, (zgl[1:]+zgl[:-1])/2.0)

    def whereIs(self, fractranNId):
        """Returns the location of a given Node Id (1-based,Fortran)"""
        nid = fractranNId-1
        if self.s:
            ix = int(nid / self.s[1])
            iz = nid % self.s[1]
            return (nid, list(self.d[nid]), [ix,iz,])
        else:
            return (nid, list(self.d[nid]))


class DotINC_File:
    """Parse and store PM element->node index"""
    def __init__(self, fn):
        self.fn = fn
        with FortranFile(fn, 'r') as fin:
            self.n = fin.read_ints()[0]
            d1_based = fin.read_ints().reshape(self.n,4)
            
            self.d = self.n*[ (0,0,0,0,), ]

            for i,tup in enumerate(d1_based):
                self.d[i] = tuple( map( lambda v: int(v-1), tup ) )

    def __iter__(self):
        for i in range(self.n):
            yield self.d[i]

    def __getitem__(self, i):
        """return a node indices (n0, n1, n2, n3) bounding element i

        Node indices have been corrected from 1-based to 0-based
        """
        return self.d[i]

class _FxElemContainer():
    """A struct for collecting Fx elements into "whole" fractures"""

    grid = None

    def __init__(self, vhContainer, firstElemInd):
        """Make a Fx container with one element"""
        # store container DotFZ_File
        # (will be the horizontal set, or vertical set)
        self.myCont = vhContainer

        # store list of element indices
        self.myE = [firstElemInd,]

        # copy-in Element Data
        ( self.hzone, self.zone, self.startN, self.endN, self.ap ) = \
            self.myCont[firstElemInd]

    def mergeWith(self, other):
        """Attempt to merge two containers
        if merge successful, modified self and return True
        if merge unsuccessful, return false"""

        # compare attribute values
        if (   other.hzone != self.hzone
                or other.zone != self.zone
                or other.ap != self.ap
                or other.myCont != self.myCont
        ):
            return False

        # compare+add to end
        if self.endN == other.startN :
            self.endN = other.endN
            self.myE = list( self.myE + other.myE )
        # compare+add to start
        elif other.endN == self.startN :
            self.startN = other.startN
            self.myE = list( other.myE + self.myE )
        else:
            return False

        return True

    def toTuple(self, y0=0., y1=1.):
        s = _FxElemContainer.grid.n[self.startN]
        e = _FxElemContainer.grid.n[self.endN]
        return ( s[0], e[0], y0, y1, s[1], e[1], self.ap, )

    def __repr__(self):
        tup = self.toTuple()
        return "_FxElemContainer( #el={:3d}, l=({}), ap={:.6f} )".format(
                len(self.myE),
                ','.join( '{:8.3f}'.format(v) for v in tup[0:6] ),
                tup[6]
            )

    def strHeader(self):
        hdr = ''
        s = _FxElemContainer.grid.n[self.startN]
        e = _FxElemContainer.grid.n[self.endN]
        if s[0] == e[0]:
            hdr = '       zbot       ztop        length      xpos        aperture'
        else:
            hdr = '      xleft      xright      length       zpos        aperture'

        return hdr

    def __str__(self):
        """o.fracture_length file headings:
  Vertical fractures
      i       zbot       ztop        length      xpos        aperture
      1     45.2571     47.2186      1.9615    208.7779      0.0001
      2      0.0000      7.0535      7.0535    349.3857      0.0001
...
  Horizontal fractures
      i      xleft      xright      length       zpos        aperture'
      1    286.5989    337.2598     50.6609     38.1940      0.0001
      2    234.0111    318.6383     84.6272     27.7040      0.0001
...
"""
        s = _FxElemContainer.grid.n[self.startN]
        e = _FxElemContainer.grid.n[self.endN]

        if s[0] == e[0]:
            # same x, vertical fracture
            return f'{s[1]:12.4f}{e[1]:12.4f}{e[1]-s[1]:12.4f}{s[0]:12.4f}{self.ap:12f}'
        elif s[1] == e[1]:
            # same z, horizontal fracture
            return f'{s[0]:12.4f}{e[0]:12.4f}{e[0]-s[0]:12.4f}{s[1]:12.4f}{self.ap:12f}'
        else:
            raise RuntimeError(f'Degenerate fracture: {self!r}')



class Parser_FractranBinaryGrid(OFracGrid):

    __DEBUG = False
    #__DEBUG = True

    def __init__(self, prefix):
        """Parse and store the relevant Fractran output files to make a grid"""

        # set element grid
        _FxElemContainer.grid = self

        # save file prefix
        self.prefix = prefix

        try:
            # nodes
            self.n = DotXYC_File( prefix+'o.xyc' )
            # PM elements
            # self.ePM = DotINC_File( prefix+'o.inc' )
            # zones and nodes for Fx elements
            self.eHFx = DotFZ_File( prefix+'o.hfz' )
            # zones and nodes for Fx elements
            self.eVFx = DotFZ_File( prefix+'o.vfz' )
        except:
            raise NotValidFractranInputError(
                "Couldn't open or parse one/all of {0}o.xyc, {0}o.hfz, {0}o.vfz".format(self.prefix) )

        # assume that Fractran always lists coordinates as
        # for x in x-gridlines:
        #     for z in z-gridlines:
        #         (x,z)
        nzgl = 0
        for i in range(self.n.n-1):
            if self.n[i][1] > self.n[i+1][1]:
                nzgl = i+1
                break

        glx = list( v[0] for v in self.n[::nzgl] )
        gly = [ 0., 1., ]
        glz = list( v[1] for v in self.n[:nzgl] )


        #import pdb; pdb.set_trace()

        # make horizontal fractures
        hfx = []
        if self.eHFx.n > 0:
            hfx = [ _FxElemContainer(self.eHFx,0), ]
            for i in range(1,len(self.eHFx)):
                tmp = _FxElemContainer( self.eHFx, i )
                if not hfx[-1].mergeWith( tmp ):
                    hfx.append( tmp )

        # make vertical fractures
        vfx = []
        if self.eVFx.n > 0:
            vfx = [ _FxElemContainer(self.eVFx,0), ]
            for i in range(1,len(self.eVFx)):
                tmp = _FxElemContainer( self.eVFx, i )
                if not vfx[-1].mergeWith( tmp ):
                    vfx.append( tmp )

        #import pdb; pdb.set_trace()

        if Parser_FractranBinaryGrid.__DEBUG:
            # mimic printing the lengths file
            print('  Vertical fractures\n      i'+vfx[0].strHeader())
            for i,f in enumerate(vfx,start=1):
                print(f'{i:6}{f!s}')
            print('\n  Horizontal fractures\n      i'+hfx[0].strHeader())
            for i,f in enumerate(hfx[0:],start=1):
                print(f'{i:6}{f!s}')

            #self._printInputs(hfx,vfx,5)
            #self._printInputs(hfx,vfx)


        self.fxgrid = OFracGrid(
            domainOrigin=(glx[0], gly[0], glz[0]),
            domainSize=(glx[-1]-glx[0], gly[-1]-gly[0], glz[-1]-glz[0]),
            fx=map( lambda fx: fx.toTuple(gly[0],gly[1]), itertools.chain(hfx,vfx) ),
            gl=[ glx, gly, glz, ]
        )


    def _printInputs(self, hfx, vfx, abbr=0):

        width=int(1+log10(max(len(hfx),len(vfx))))

        for k,fxSet in (('h',hfx), ('v',vfx),):
            nSet = len(fxSet)
            b = nSet

            if abbr > 0:
                b = abbr

            for i in range(b):
                print( '{0}{1: <{w}} {2}'.format(k,i,fxSet[i],w=width) )

            if abbr > 0:
                print('...')
                for i in range(nSet-b,nSet):
                    print( '{0}{1: <{w}} {2}'.format(k,i,fxSet[i],w=width) )

    def getOFracGrid(self):
        return copy.deepcopy( self.fxgrid )


def iterFractranParsers():
    """Generate parser Class objects available in this module"""
    yield Parser_FractranBinaryGrid
    yield Parser_FractranLengthsFile

