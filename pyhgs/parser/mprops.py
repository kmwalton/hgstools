#!/usr/bin/env python
"""Parse Hydrogeosphere .mprops files

Parses material types and their properties into dict objects.

Limitations: only supports simple property=float -type properties. Table-type
properties will cause an error.

A main method will spew the properties to the console or write the dictionary to
a python pickle file.

"""

import sys
import os
import re
import io
import pprint
import argparse
import pickle
from more_itertools import grouper

def _parse_mat(mat_str):
   """Parse simple name\\nvalue pairs"""
   d = dict()

   for k,v in grouper(mat_str.strip().split('\n'),2):
      d[k.lower()] = float(v)

   return d

def parse(file, do_skips=False):
   """Return a dict of material types with dicts of properties


   Material types dict is keyed by material name
   Properties dicts are keyed by property name, where property names are
   converted to lowercase and spaces are turned into underscores.

   Arguments:
      file : str or file-like

      do_skips : bool
         Ignore lines in skip on/skip off blocks, as HGS would

   """

   file_text = ''
   if type(file) is str:
      with open(file,'r') as fin:
         file_text = fin.read().strip()
   elif isinstance(file, io.TextIOBase):
      file_text = file.read().strip()
   else:
      raise ValueError('Expecting string or text stream in "file",'\
              f' not {type(file)}')

   # get rid of comments
   file_text = re.sub('\!.*?\n','\n',file_text)

   # get rid of trailing whitespace, leading space, multiple newlines
   file_text = re.sub(r'\s*\n\s*','\n',file_text)

   # parse-out skip on/skip off lines so non-standard definitions prevail
   if do_skips:
       file_text = re.sub(
               r'\s*skip\s+on.*?skip\s+off\s*', '\n',
               file_text,flags=re.DOTALL)
   else:
       file_text = re.sub(r'\s*skip\s+(on|off)\s*','\n',file_text)

   mats = dict( map(
         lambda t: (t[0],_parse_mat(t[1]),),
         re.findall(
            r'^(?P<name>\S+)\n(?P<body>.*?)end material$',
            file_text,
            re.MULTILINE|re.DOTALL) ) )

   return mats


if __name__=='__main__':
    
    argp = argparse.ArgumentParser()
    argp.add_argument('FILE',type=argparse.FileType('r'),)
    argp.add_argument('FILE_OUT', nargs='?',
            default=None,
            type=argparse.FileType('wb'),)
    args = argp.parse_args()

    mats = []

    while True:

        try:
            m = Material(args.FILE)
        except:
            break
        else:
            mats.append(m)

    mats = parse(args.FILE)

    fn = os.path.basename(args.FILE.name)
    print(f'Found {len(mats)} materials in {fn}:')
    pprint.pp(mats)

    if args.FILE_OUT:
        pickle.dump(d)

    sys.exit(0)


