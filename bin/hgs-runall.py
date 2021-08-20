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
        help='List the steps in the toolchain and verify executability of each '\
         'tool; check file read and write permissions in the simulation ' \
         'directory'
        )

    argp.add_argument('-sw', '--start-with',
        type=str,
        default=None,
        help="""Start the toolchain with tool N. See --dry-run for integer
        index values for each tool, or specify a word matching a tool name, like
        'grok' or 'postprocess'. The chain will begin with the matching tool with
        the lowest number.
        """
        )

    argp.add_argument('-ew', '--end-with',
        type=str,
        default=None,
        help="""End the toolchain after tool N finishes. See --dry-run for
        index values for each tool, or specify a word matching a tool name, like
        'hsplot' or 'preprocess'. The chain will end after the matching tool with
        the highest number.
        """
        )

    args = argp.parse_args()

    tc = HGSToolChainRun(args.sim_dir)

    tc.set_start(args.start_with)
    tc.set_end(args.end_with)

    if args.dry_run:
        print(tc)
        (retval,msgs) = tc.check_tools()
        if retval != 0:
            print( 'Error with one or more executables in the toolchain:\n' \
                + '\n'.join(msgs), file=sys.stderr )
        sys.exit(retval)

    
