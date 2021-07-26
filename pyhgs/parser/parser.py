#!/usr/bin/env python
"""
main method herein will dump a binary file to console output.
"""

import sys
import re
import argparse
import pprint
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
    argp.add_argument('FILE_NAME',
        help='The name of the file to parse')
    # TODO add output filename
    # TODO add pickle vs text spew option
    args = argp.parse_args()

    fdata = None
    try:
        fdata = parse(args.FILE_NAME)
    else:
        sys.exit(1)

    # TODO print unabridged arrays
    pprint.pprint(fdata)
    sys.exit(0)

