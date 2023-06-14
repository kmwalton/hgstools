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
import copy
from more_itertools import grouper

# idea from
# https://stackoverflow.com/questions/3296499/case-insensitive-dictionary-search
class CaseInsensitiveDict(dict):
    def __init__(self, d=dict()):
        self._d = d
        self._s = dict((k.lower(), k) for k in d.keys())
        if len(self._s) != len(self._d):
            raise ValueError('>=2 keys in d map to the same value')
    def __contains__(self, k):
        return k.lower() in self._s
    def __len__(self):
        return len(self._s)
    def __iter__(self):
        return iter(self._s)
    def __getitem__(self, k):
        return self._d[self._s[k.lower()]]
    def actual_key_case(self, k):
        return self._s.get(k.lower())
    def __setitem__(self, k, v):
        self._d[k] = v
        self._s[k.lower()] = k
    def __str__(self):
        return str(self._d)

def _parse_mat(mat_str):
   """Parse simple name-value pairs"""
   d = dict()

   for k,v in grouper(mat_str.strip().split('\n'),2):
      d[k] = float(v)

   return CaseInsensitiveDict(d)
   #return d

def parse(file, do_skips=False):
   """Return a dict of material types with dicts of properties

   Material types dict is keyed by material name, (and delivered in a case
           insensitive, dictionary-like datatructure)
   Properties dicts are keyed by property name, where property names are
   converted to lowercase and spaces are turned into underscores.

   Arguments:
      file : str or file-like

      do_skips : bool
         Ignore lines in skip on/skip off blocks, as HGS would.

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

   mats = dict(map(
         lambda t: (t[0].strip(),_parse_mat(t[1]),),
         re.findall(
            r'^(?P<name>\S.*)\n(?ms:(?P<body>.*?))^\s*(?i:end material)$',
            file_text,
            re.MULTILINE) ) )

   mats = CaseInsensitiveDict(mats)

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


