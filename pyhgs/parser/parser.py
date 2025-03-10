#!/usr/bin/env python
"""
main method herein will dump a binary file to console output.

.TODO..

This old functionality should be merged in with `pyhgs.parser.eco` and `pyhgs.parser.grok`
"""

import sys
import os
import re
import argparse
import pprint
import logging
from itertools import count
from bisect import bisect_left,bisect
from decimal import Decimal as D

import numpy as np
from scipy.io import FortranFile

import tabulate

from . import parse

__docformat__ = 'numpy'

class ParserHGSEcoFile:

    def __init__(self,p):
        self.prefix = p

        lines = []
        with open(f'{p}o.eco','r') as fin:
            lines = fin.readlines()

        
        reOTimes = re.compile(r'^OUTPUT TIME: *([0-9.]+)')
        # TODO add lookahead assertion that 'e' is not part of 'end'
        #reOTimes = re.compile(r'^OUTPUT TIME: *([0-9.eEdD+-]+)')

        self.outputTimes = list(
            D(reOTimes.match(l).group(1))
            for l in lines if reOTimes.match(l)
            )


class GrokParser():

    def __init__(self,prefix=''):
        """
            Arguments:
                prefix : str
                    A prefix, or path/prefix, combined.
                    Or the empty string will look for the prefix in batch.pfx
        """
        if not prefix:
            with open('batch.pfx','r') as batchfin:
                self.dir = '.'
                self.pfx = batchfin.readline().strip()
        else:
            self.dir = os.path.dirname(prefix)
            self.pfx = os.path.basename(prefix)
            if self.pfx.lower().endswith('.grok'):
                self.pfx = self.pfx[:-5]

        self.grokfn = os.path.join(self.dir,f'{self.pfx}.grok')

        # TODO
        #   make a line number+file index for keeping track of includes
        with open(self.grokfn,'r') as fin:
            self._groktxt = fin.read()

        inc_fn = list( re.findall(r'include\s+(\S+)',self._groktxt,
                    flags=re.I|re.MULTILINE) )

        self._includes = dict((m,'',) for m in inc_fn)

        #for fn in self._includes:
        #    with open(fn,'r') as fin:
        #        self._includes[fn] = fin.read()
        
    def __str__(self):
        s = f'{self.grokfn}'
        if self._includes:
            s += ' with includes: '+', '.join(self._includes.keys())
        return s

    def get(self, what):
        return None

def _parse_print_binary(fn):

    fdata = parse(fn)

    for k in fdata:
        # header line
        ss=f''
        if hasattr(fdata[k],'shape') and len(str(fdata[k].shape))>2:
            ss = f' shape={fdata[k].shape}'

        # options
        if args.unabridged:
            np.set_printoptions(threshold=np.inf)

        # data
        print(f'{k}{ss}:\n{fdata[k]}')
  
def _parse_print_grok(fn):
    gp = GrokParser(fn)
    print(gp)

if __name__ == '__main__':

    argp = argparse.ArgumentParser()
    argp.add_argument('-v','--verbose',
        action='count',
        default=0,
        help='Increase verbosity.',
        )
    argp.add_argument('-u','--unabridged',
        default=False,
        action='store_true',
        help='Print full arrays of values. Default: print truncated arrays.'
        )
    # TODO add output filename
    # TODO add pickle vs text spew option
    argp.add_argument('FILE_NAME',
        help='The name of the file to parse')
    args = argp.parse_args()

    plogger = logging.getLogger('pyhgs.parser')
    plogger.addHandler(logging.StreamHandler())
    if args.verbose > 2:
        plogger.setLevel(logging.DEBUG)
    elif args.verbose > 1:
        plogger.setLevel(logging.INFO-2)
    elif args.verbose > 0:
        plogger.setLevel(logging.INFO)
        

    errmsg = ''
    try:
        _parse_print_binary(args.FILE_NAME)
    except RuntimeError as e:
        errmsg = str(e)

    if errmsg:
        try:
            _parse_print_grok(args.FILE_NAME)
            errmsg = '' # success flag
        except Exception as e:
            errmsg += '\n\n' + str(e)

    if errmsg:
        print(errmsg, file=sys.stderr)
        sys.exit(1)


    sys.exit(0)

