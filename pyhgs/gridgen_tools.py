#!/usr/bin/env python3
"""Generates grid in one dimension for HGS or CF format-input

Originally written for ssfl, 1D, unsaturated zone column simulations
Thu, Feb  4, 2016  4:56:22 PM

"""
import sys,os,datetime,argparse
from math import log10,floor
from itertools import zip_longest
from functools import partial

def regularLines(spacing=0.1, start=0.0, end=30.0, irregGLs=[]):
   """Grid lines from 0 m to 30m in irregular intervals
   """
   gl = irregGLs + [start]

   nextgl = start+spacing
   while nextgl < end-1.0e-4:
      gl.append(nextgl)
      nextgl += spacing

   gl.append(end)

   gl.sort()
   return gl


def printHGS(gl, fout=sys.stdout):
   """Print the given array, plus metadata/comments, in HGS format

   Parameters
   ----------
   gl : array-like
      An array of ascending numeric values representing the location of grid
      lines.
   fout : file-like
      An open file/stream where output text will be written.
   """

   pfout = partial(print, file=fout)


   pfout("! Grid lines printed by {}, {}".format(
               os.path.basename(__file__), datetime.datetime.now()) )
   pfout("! {} grid lines; {} total distance ".format(len(gl), gl[-1]-gl[0]))

   gl.reverse()

   WIDTH=9 # sign + 4-digits + decimal + 3-digits
   PRECISION=3
   FMT="{{:{}.{}f}}".format(WIDTH,PRECISION)

   #HGS line buffer 300 chars, less \r\n, less one number
   # I'm no sure where '300' came from. From experience,
   # the maximum line width is much lower. Hence, 80.
   MAX_LINE=80 - 2 - WIDTH

   pfout( len(gl) )

   while gl:
      pfout( FMT.format(gl.pop()), end='' )
      cc=WIDTH # character count
      while gl and cc<MAX_LINE and not float(gl[-1]) % 5.0 < 1e-6 :
         pfout( " "+FMT.format(gl.pop()), end='' )
         cc += 1+WIDTH

      pfout()

def printRFGen(gl):
   """print grid lines (in format of RFGen); a 'i x y z'-shaped list, in columns"""


   # print i x y z list of gridlines

   maxGLdim=0
   for k in gl:
      print("! {} {}-grid lines; {} total distance ".format(
               len(gl[k]), k, gl[k][-1]-gl[k][0]))
      maxGLdim = max(maxGLdim, len(gl[k]))


   WIDTH=9
   PRECISION=3
   FFMT="{{:{}.{}f}}".format(WIDTH,PRECISION)
   FMT = "{{:-{}d}} {} {} {}".format(floor(log10(maxGLdim))+1, FFMT, FFMT, FFMT)

   MAX_LINE=80 - 2 - WIDTH # keep less than 80 characters

   for i,(x,y,z) in \
      enumerate(zip_longest(gl['x'], gl['y'], gl['z'], fillvalue=0.0),start=1 ):
      print( FMT.format( i, x, y, z ) )


def printCF(gl):
   """print grid lines (suitable for CF)"""

   INDW=int(floor(log10(len(gl))) + 1)
   GLPREC=3 #TODO: scan for the actual input precision (assume 3 decimal digits enough)
   GLDW=int(floor(log10(gl[-1]))+1)+GLPREC+1 # assume largest magnitude at end
   FMT=" {{:{0}d}} {{:{0}d}} {{:{1}.{2}f}}".format(INDW,GLDW,GLPREC)

   print("// Grid lines printed by {}, {}".format(
               os.path.basename(__file__), datetime.datetime.now()) )
   print("// {} CVs; {} total distance ".format(len(gl)-1, gl[-1]-gl[0]))

   def genCvSize( g ):
      for i in range(len(g)-1):
         dg = g[i+1]-gl[i]
         if dg == 0.0:
            raise ValueError("Found 0.0 distance between adjacent grid lines {} and {}+1".format(i) )
         yield dg

   dgen=genCvSize(gl)

   istart=0
   dlast=next(dgen)
   dsum=dlast
   for i,d in enumerate(dgen,start=1):
      dsum+=d
      if abs(d - dlast) < 1e-7:
         pass
      else:
         print( FMT.format(istart+1,i,dlast) )
         istart=i
         dlast=d
   print( FMT.format(istart+1,len(gl)-1,dlast) )
   print("end\n// total output distance {}".format(dsum))


if __name__ == "__main__":
   #gridlines = problem_1DfracturedColumn_z()
   #printHGS(gridlines)
   #printCF(gridlines)

   parser = argparse.ArgumentParser()
   parser.add_argument( '--nudge-fx-coordinates-to', type=float,
         help='Nudge each fracture to the nearest interval' )
   args = parser.parse_args()

   defspacing=0.25
   zstart = -122.0
   zend = -2.0

   gl = { 'x':[], 'y':[], 'z':[] }
   gl['x']=regularLines(
         spacing=defspacing, start=0.0,end=0.3, irregGLs=[] )
   gl['y']=regularLines(
         spacing=defspacing, start=0.0,end=0.3, irregGLs=[] )
   gl['z']=regularLines(
         spacing=defspacing,
         start=zstart,end=zend,
         irregGLs=[zstart+0.01, zend-0.02, zend-0.01] )

   #printRFGen(gl)
   #printCF(gl)
   for c in 'xyz':
      printHGS(gl[c])


   sys.exit(0)
