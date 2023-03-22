#!/usr/bin/env pythin
"""Erase all files that *are not* in the keep file's list"""

import argparse
import os
from pyhgs.runner import PyPowerShellRunner


if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description='''Erase all files in the current directory that *are not*
        in the keep file's list''',
    )
    ap.add_argument( 'keep_file',
        type=str,
        help='The name of a file (relative or absolute path) that contains a '\
            'list of file names and/or glob-style patterns of files to '\
            'retain in the present directory.',
        )
    args = ap.parse_args()

    # check for keep file
    if args.keep_file and not os.path.isfile(args.keep_file):
        ap.error(f'Could not find keep file {args.keep_file}')

    keep_file_list = PyPowerShellRunner.read_keep_file(args.keep_file)

    runner = PyPowerShellRunner(None)
    runner.preexistingfiles = []
    runner.eraseRunOutputs(keep=keep_file_list)
