"""A parser for RFGen grids from RFGen text output"""

import sys
import re
import copy
from ofracs import OFracGrid

class NotValidRFGenInputError(Exception):
    """Exception raised if input file is not valid"""
    pass

class RFGenFracListParser:

    def __init__(self, fnin):
        
        with open(fnin,'r') as fin:
            fracsHere = RFGenFracListParser.read_fracs(fin)

        self.fxnet = OFracGrid(
            fx=list( f[0][:7] for f in fracsHere ) )

    def getOFracGrid(self):
        return copy.deepcopy( self.fxnet )


    @staticmethod
    def read_fracs(fin):
        """Read a list of fractures with an 'xfrom xto ...'-style header and list


        Arguments:
            fin : file-like object
                An open file pointing near where the header should be
        """

        line = fin.readline()
        header_match = None

        # read through to header definition
        # the (first) index column, 'i', is optional
        # the (last) type column 'ifractyp' or 'type' is optional
        fracListHdr = re.compile(
         r'(i)?\s*(xfrom\s+xto\s+yfrom\s+yto\s+zfrom\s+zto)\s*(type|ifractyp)?',
         flags=re.I)

        while True:
           header_match = fracListHdr.search(line)
           if header_match:
              break
           line = fin.readline()
           if line == '':
              break

        if not header_match:
           raise NotValidRFGenInputError(
              "Could not find the required fracture list header line in "
              +fin.name)

        # see if there's a column for fracture index number
        _readI = header_match.groups()[0] == 'i'
        _readType = header_match.groups()[-1] == 'type' \
                    or header_match.groups()[-1] == 'ifractyp'

        #
        #  Read the axis-aligned, 2d fractures in this data file
        #  Store them in a list of tuples:
        #
        # fracsHere =
        #    [ ...
        #      ( [ Xfrom,Xto, Yfrom,Yto, Zfrom,Zto, Aperture ], type ),
        #      ... ]
        #
        fracsHere = []

        # TODO:
        #  Use REs for reading fracture lines.
        #  break when !fin or no match ... more robust for different invocations
        #  of RFGen (RFGen vs HGS embedded)

        for line in fin.readlines():
           line = line.strip()

           # strip out comments
           COMMENT = ['//', '!'] # c-style, and HGS-style
           for COM in COMMENT:
               if COM in line:
                  line = line.partition(COMMENT)[0].strip()

           # skip blank lines
           if line == "" : continue

           # a hack to exit when reading a HGS eco file
           if line.startswith('GRID'): break   

           # split up the coordinate strings into an array
           t = line.split()
           
           if _readI: t = t[1:] # junk the fracture index

           o = -1 # set to invalid value
           if _readType:
               o = int(t[-1]) 
               t = t[:-1]

           fracsHere.append( ( tuple( map(float, t) ), o, ) )

        return fracsHere



class RFGenOutFileParser:

   def __init__(self, fnin):
      #data file given on command line
      fin = None
      try:
          fin = open( fnin, 'r' )
      except FileNotFoundError:
          extramsg = ''
          if fnin.lower() != "report-rfgen.txt":
              extramsg = " Did you mean to use 'report-rfgen.txt'?"
          raise NotValidRFGenInputError(f'Could not open {fnin}.{extramsg}')

      self._fnin = fnin

      m = None    # re match object
      line = ""   # line of input file

      gridlines = [ [], [], [] ]

      gridEntry= re.compile(r'^[0-9]+\s+([+-.0-9]+)\s+([+-.0-9]+)\s+([+-.0-9]+)$')
      #                            i         xi            yi           zi

      #import pdb ; pdb.set_trace()

      # find the first grid line entry
      line = fin.readline()
      while fin and line:
         m = gridEntry.match(line.strip())
         if m:
            break
         line = fin.readline()
         if line == '':
            break

      if not m:
         raise NotValidRFGenInputError("Could not find the grid line list in "+fnin)

      # read grid lines until the end of the list (signified by any line that
      # doesn't match the re pattern
      spot = 0
      while fin and m:
         for xyz in range(3):
            gridlines[xyz].append(float(m.group(xyz+1)))
         spot = fin.tell()
         line = fin.readline().strip()
         m = gridEntry.match(line)

      fin.seek(spot) # rewind to the line that did not match the pattern

      # fix the lists of gridlines (remove the zero values that pad the end of
      # two of the three lists)
      for xyz in range(3):
         listEndIndex=len(gridlines[xyz])
         for i in range(listEndIndex-1):
            if gridlines[xyz][i+1] < gridlines[xyz][i]:
               listEndIndex = i+1
               break
         gridlines[xyz] = gridlines[xyz][0:listEndIndex]

      fracsHere = RFGenFracListParser.read_fracs(fin)

      fin.close()

      # make new grid lines so that all fractures fall on a gridline
      for a in range(3):
         gla = set( gridlines[a] )
         sizeBefore = len(gla)
         for f in fracsHere:
            gla.update( f[0][2*a:2*a+2] )
         gridlines[a] = list( gla )
         gridlines[a].sort()
         if sizeBefore != len(gridlines[a]):
            print(
               "Fracture coordinates made additional gridlines in {}.".format(
                  'XYZ'[a],
               file=sys.stderr ) )


      self.fxnet = OFracGrid(
         gl=gridlines,
         fx=list( f[0][:7] for f in fracsHere ) )

      self.fxnet.metadata['originaldatafile'] = self._fnin

   def getOFracGrid(self):
      return copy.deepcopy( self.fxnet )

