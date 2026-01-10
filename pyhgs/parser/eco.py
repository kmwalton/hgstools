"""A parser for RFGen grids from a grok .eco file





"""

import os
import re
import copy
from decimal import Decimal as D


from .._import_helpers import _get_from_ofracs

class NotValidHGSEcoInputError(Exception):
    """Exception raised if input file is not valid"""
    pass

_tags = {
    'zones_pm_prop':(
        lambda txt:
           list( map( 
             lambda v: (int(v[0]), v[1].strip(),),
             re.findall(
              r'Porous media zone\s+(\d+)\s+' +
              r'has properties of(?ms:\s)+(\S.*)',
              txt, flags=re.I)
           ))
        ),
    'zones_pm':(lambda txt:
        max([1,] + list(map(int, re.findall(
            r'Total number of porous media zones now\s+(\d+)',
            txt, flags=re.I))))
        ),
    'units':(lambda txt:
        re.search(r'GROK:\s+units:\s+(\w+?)-(\w+?)-(\w+?)',
            txt, flags=re.I).group(1,2,3)
        ),
}

_unit_abbr={
    'kilogram':'kg',
    'kg':'kg',

    'meter':'m',
    'metre':'m',
    'm':'m',

    'day':'d',
    'd':'d',
    'hour':'hr',
    'hr':'hr',
}

class EcoFile:
    
   def __init__(self, fnin):
      #data file given on command line
      if not fnin.endswith('o.eco'):
            fnin += 'o.eco'

      if not os.path.isfile(fnin):
         raise NotValidHGSEcoInputError(f'Could not open {fnin}')

      self._fnin = fnin

      # load whole file into memory
      with open(self._fnin) as fin:
         self.txt = fin.read()

   def get_n_zones(self):
      '''Return the number of 'new zones created.'''
      return _tags['zones_pm'](self.txt)

   def get_pm_zone_properties(self):
      '''Return a list of (zone id#, material name)'''
      return _tags['zones_pm_prop'](self.txt)

   def get_output_times(self):
      '''Return a list of output times'''

      # older hydrogeosphere versions

      reOTimes = re.compile(r'OUTPUT TIME:\s+([0-9.]+)')
      # TODO add lookahead assertion that 'e' is not part of 'end'
      #reOTimes = re.compile(r'^OUTPUT TIME: *([0-9.eEdD+-]+)')

      times_raw = reOTimes.findall(self.txt)

      if not times_raw:
        times_raw = self._find_output_times_v2829(self.txt)

      self.outputTimes = list(map(D, times_raw))

      return self.outputTimes

   def get_units(self):
       """Return the (mass, length, time) units of the problem as a tuple"""

       u = _tags['units'](self.txt)

       return [ _unit_abbr[_] for _ in u ]

   def _find_output_times_v2829(self, txt):
      """Finds the 'OUTPUT TIMES' block and extracts the table data.

      This function looks for the '# OUTPUT TIMES' header, then finds the
      'Step # Timestep Time' line. It captures all text after that line
      until it encounters a line of dashes ('-------------------').

      The captured data is stored in a named group 'table_values'.

      Args:
         txt: The string content to search within.

      Returns:
         A list of the times it found, or an empty list if none found
      """
      # This pattern uses re.DOTALL (the 's' flag), which makes '.' match
      # newline characters. This is crucial for matching across multiple lines.
      #
      # Breakdown:
      # # OUTPUT TIMES       - Matches the literal starting header.
      # .*?                  - Lazily matches any character (including newlines)
      #                        until the next part of the pattern is found.
      # Step #\s+Timestep\s+Time\s*\n - Matches the table header line and the
      #                                  newline after it.
      # (?P<table_values>.*?) - Starts a named capture group 'table_values'.
      #                        It lazily captures all characters until...
      # \n-------------------  - It finds a newline followed by the dash
      #                        delimiter you specified.

      pattern = r"# OUTPUT TIMES.*?Step #\s+Timestep\s+Time\s*\n(?P<table_values>.*?)\n-------------------"

      # re.DOTALL makes the '.' special character match any character,
      # including a newline.
      match = re.search(pattern, txt, re.DOTALL)

      times = []

      if match:
          times = list(
             float(l.strip().split()[-1])
               for l in match['table_values'].strip().split('\n'))

      return times

   def getOFracGrid(self):

      # delay the fetching/import of this the ofrac module in case the user does
      # not have it.
      locals().update(_get_from_ofracs('OFracGrid'))

      fin = open(self._fnin, 'r')

      m = None    # re match object
      line = fin.readline()

      gridlines = [ [], [], [] ]
      gridEntry= re.compile(r'^[0-9]+\s+'+r'\s+'.join(3*['([+.0-9-]+)',])+'$')

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
        raise NotValidHGSEcoInputError(
            "Could not find the grid line list in "+self._fnin)
 
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
         r'(i)?\s*(xfrom\s+xto\s+yfrom\s+yto\s+zfrom\s+zto)\s+aperture\s*(type|ifractyp)?')
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
      return copy.deepcopy( self.fxnet )

