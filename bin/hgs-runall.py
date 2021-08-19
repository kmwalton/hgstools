#!/usr/bin/env python3
"""Launches the extended Hydrogeosphere toolchain.

Basic Toolchain:
1) grok
2) phgs
3) hsplot

Extensions include:
0) preprocess* | sort
1..3)
4) postprocess* | sort

Note: this program captures the standard output and errors streams. This may
interfere with a debugger.
"""

import sys
import os
import glob
import argparse
import subprocess

from pyhgs.runner import HGSToolChainRun

if __name__ == '__main__':

    argp = argparse.ArgumentParser()

    argp.add_argument('-d','--sim-dir',
        default='.',
        type=str,
        help='Directory containing simulation inputs',
        )

    argp.add_argument('--dry-run',
        action='store_true',
        help='List the steps in the toolchain. Verify executability of each '\
         'step, check file read/write permissions, etc.'

    args = argp.parse_args()

    tc = HGSToolChainRun(args.sim_dir)

    if args.dry_run:
        print(tc)
        (retval,msgs) = tc.check_tools()
        if retval != 0:
            print( '\n'.join(msgs), file=sys.stderr )
        sys.exit(retval)
    
