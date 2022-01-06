"""A parser for RFGen grids from a grok .eco file"""

import re,sys,copy
from ofracs import *

class NotValidHGSEcoInputError(Exception):
    """Exception raised if input file is not valid"""
    pass


class HGSEcoFileParser:

    
   def __init__(self, fnin):
      #data file given on command line
      if not fnin.endswith('o.eco'):
            fnin += 'o.eco'
      fin = None
      try:
          fin = open( fnin, 'r' )
      except FileNotFoundError:
          raise NotValidHGSEcoInputError(f'Could not open {fnin}')
      self._fnin = fnin

      m = None    # re match object
      line = fin.readline()

      gridlines = [ [], [], [] ]
      gridEntry= re.compile('^[0-9]+\s+'+r'\s+'.join(3*['([+.0-9-]+)',])+'$')

      #
      # READ THE GRID
      #
      # find the grid line dataheader
      # "x,y,z  grid line locations"
      while fin and line:
        if line.strip() == 'x,y,z  grid line locations':
            break
        line = fin.readline()

      # find the first grid line entry
      line = fin.readline() # first grid line
      m = gridEntry.match(line.strip())

      if not m:
        fin.close()
        raise NotValidHGSEcoInputError("Could not find the grid line list in "+fnin)
 
      # read grid lines until the end of the list (signified by any line that
      # doesn't match the re pattern
      while fin and m:
         for xyz in range(3):
            gridlines[xyz].append(float(m.group(xyz+1)))
         line = fin.readline().strip()
         m = gridEntry.match(line)

      # fix the lists of gridlines (remove the zero values that pad the end of
      # two of the three lists)
      for xyz in range(3):
         listEndIndex=len(gridlines[xyz])
         for i in range(listEndIndex-1):
            if gridlines[xyz][i+1] < gridlines[xyz][i]:
               listEndIndex = i+1
               break
         gridlines[xyz] = gridlines[xyz][0:listEndIndex]

      #
      # READ THE Fractures
      #
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

      # read through to header definition
      # the (first) index column, 'i', is optional
      # the (last) type column 'ifractyp' or 'type' is optional
      fracListHdr = re.compile(
         '(i)?\s*(xfrom\s+xto\s+yfrom\s+yto\s+zfrom\s+zto)\s+aperture\s*(type|ifractyp)?')
      while fin:
         m = fracListHdr.search( line.lower() )
         if m:
            self.originalHeader = line
            break
         line = fin.readline()

      if not m:
         raise NotValidRFGenInputError(
            "Could not find the required fracture list header line in "+fnin)

      # see if there's a column for fracture index number
      _readI = m.groups()[0] == 'i'
      _readType = m.groups()[-1] == 'type' or m.groups()[-1] == 'ifractyp'
 

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
         
         o = -1 # set to invalid value
         if _readType:
             o = int(t[-1]) 
             t = t[:-1]

         # these will be made into OFrac objects, which want 
         # xfrom, xto, yfrom, yto, zfrom, zto, ap
         fracsHere.append( tuple(map(float, t)) )
 

      fin.close()

      # make the OFracGrid object
      self.fxnet = OFracGrid( gl=gridlines, fx=fracsHere )

   def getOFracGrid(self):
      return copy.deepcopy( self.fxnet )
