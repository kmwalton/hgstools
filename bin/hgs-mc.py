#!/usr/bin/env python
"""Run Monte Carlo-style realization of a Hydrogeosphere problem

Examples:
    Run 5 variants on ../base_sim/. Creates ./mc1, ..., ./mc5. 
    $ python hgs-mc.py 5 ../base_sim

Dependencies:
    - hgs-copy-inputs.py (for copying base_sim)
    - pre- and post-processing scripts in base_sim to set up a new realization

    - pycf.ancillary.TextManip

"""
import sys
import os
import argparse
import shutil
import shlex

from multiprocessing import Pool

import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

# Useful:
#https://wiseodd.github.io/techblog/2016/06/13/parallel-monte-carlo/

class HGS_MCRunner():
    
    def genInputs(self, n):
        for i in range(n):
            s=f'mc{i:0>3d}'
            logger.info(f'making {s}')
            yield s

    def runSingle(self, d):
        logger.debug(f'running {d}')
        return 0


if __name__ == '__main__':

    ap = argparse.ArgumentParser()

    ap.add_argument( 'nRuns', type=int )

    ap.add_argument( 'base_sim', type=str,
        help="Directory containing the base simulation inputs.")

    ap.add_argument( '-v', '--verbose', default=1, action='count',
        help='Set verbosity. Default 1. Set to 0 for quiet, or higher for '\
            'increasing levels of chatter.'
        )

    ap.add_argument( '--num-processes', default=1, type=int,
        help="The maximum number of concurrent processes to use.")

    ap.add_argument( '--copy-command',
        default=shutil.which('hgs-copy-inputs.py'),
        type=str,
        help='The name of the script or a shell command that clones the ' \
            'base simulation inputs in a MC directory.')

    ap.add_argument( '--run-command',
        default=shutil.which('hgs-runall.ps1'),
        type=str,
        help='The name of the script or a shell command that that invokes ' \
            'preprocessing, HGS, and postprocessing/cleanup',
        )

    args = ap.parse_args()

    if args.verbose < 1:
        logger.setLevel(logging.WARN)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)

    # check for base dir
    if not os.path.isdir(args.base_sim):
        print(f'Could not find base simulation directory "{args.base_sim}".',
            file=sys.stderr)
        sys.exit(-1)
    # autodetect the batch.pfx file
    if not 'batch.pfx' in os.listdir(args.base_sim):
        print(f'Could not find "{args.base_sim}{os.path.sep}batch.pfx".',
            file=sys.stderr)
        sys.exit(-1)

    # Check if other executable files exist
    # May need posix=False here, if retaining quotes in individual
    # parameters is important
    for f in [ args.copy_command, args.run_command ]:
        f = shlex.split(f.strip())[0]
        if not shutil.which(f):
            print('Could not find {f}.', file=sys.stderr)
            sys.exit(-1)

    mc = HGS_MCRunner()

    # run
    if args.num_processes > 1:
        with Pool(args.num_processes) as p:
            p.map(mc.runSingle, mc.genInputs(args.nRuns))
    else:
        for r in mc.genInputs(args.nRuns):
            mc.runSingle(r)
