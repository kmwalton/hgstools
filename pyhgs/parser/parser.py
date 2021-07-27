#!/usr/bin/env python
"""
main method herein will dump a binary file to console output.
"""

import sys
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

from pyhgs.parser import parse

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
        

    fdata = None
    try:
        fdata = parse(args.FILE_NAME)
    except Exception as e:
        print(e,file=sys.stderr)
        sys.exit(1)

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

    sys.exit(0)

