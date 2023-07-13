#!/usr/bin/env python

import os
import argparse

from pyhgs.mesh import HGSGrid
from pyhgs.calcs avconc
import pyhgs.cli

if __name__=='__main__':
    print('hello world')

    argp = argparse.ArgumentParser()

    argp.add_argument('PATH_TO_PREFIX',
        metavar='path/to/prefix',
        action=pyhgs.cli.PathToPrefix,
        help='path/to/prefix of the HGS simulation',
    )

    args = argp.parse_args()

    # load grid
    path_to, prefix, leftovers = PathToPrefix.split(args.PATH_TO_PREFIX)
    g = HGSGrid(path_to+os.sep+prefix)

    # parse spatial location

    # parse weighting scheme
    # compute weights


    # get static data

    # parse times for computation

    # iterate through time series


    # produce output
