#!/usr/bin/env python
"""Erase all files that *are not* in the keep file's list"""

import argparse
import os
from pathlib import Path
from operator import itemgetter
from itertools import (count, filterfalse)
from fnmatch import fnmatch

from hgstools.pyhgs.runner import BaseRunner

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

def ravel_merge_nested(nested, flat_dest, sep='.', _depth=0):
    """Ravel the nested dict of files to a new flat dict"""

    ret = {}

    for cat, v in nested.items():
        if isinstance(v, set):
            for fn in v:
                if not fn in flat_dest: # don't overwrite!
                    ret[fn] = cat
        elif isinstance(v, dict):
            rret = ravel_merge_nested(v, flat_dest, sep, _depth+1)
            for fn, rcat in rret.items():
                rret[fn] = cat + sep + rcat
            ret.update(rret)
        else:
            raise NotImplementedError()

    if _depth == 0:
        flat_dest.update(ret)
        return flat_dest

    # else
    return ret

def _key_by_lst_mtime(runner):

    lst = Path(runner.simdir)/(runner.prefix+'o.lst')
    if lst.is_file():
        return lst.stat().st_mtime

    return 0



if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description='''Erase all files in the current directory that *are not*
        in the keep file's list''',
    )

    ap.add_argument( '--keep-file',
        type=str,
        default=None,
        metavar='path/to/KEEP_FILE',
        dest='keepfile_name',
        help='''\
            The name of a file (relative or absolute path) to a "keep file".
            A "keep file" is a text file that specifies simulation files to
            be retained. Lines of this file may be: a relative or abolute file
            name, a glob-style pattern, or a category of simulation file
            (prefixed by "cat="). Categories of files are identified in a
            --dry-run; they can be specified in plain text, e.g., "cat=input", or as
            a glob-style pattern, e.g., "cat=output.tecplot*".
        '''
        )

    ap.add_argument('--dry-run', action='store_true', default=False,
        help='''Show the categorization of all files in the current directory and
        whether whether they will be deleted, given the KEEP_FILE and other
        options' values.''',
        )

    keep_inputs_grp = ap.add_mutually_exclusive_group()
    keep_inputs_grp.add_argument('--keep-inputs', action='store_true',
        dest='keep_inputs',
        default=True,
        help='''Keep files in the 'input' category (default behaviour).
        Equivalent to a line with 'cat=input' in the KEEP_FILE.''',
        )
    keep_inputs_grp.add_argument('--no-keep-inputs', action='store_false',
        dest='keep_inputs',
        help='''Delete files in the 'input' category.''',
        )


    ap.add_argument( '--make-snippits', action='store_true', default=False,
        help='''Reduce prefixo.eco and prefixo.lst to a few lines.'''
        )

    ap.add_argument('-v', '--verbose',
        default=0,
        action='count',
        help='''Increase the verbosity. Default 0, no output.'''
        )

    ap.add_argument('pfx',
        metavar='PREFIX',
        nargs='*',
        default=['.',],
        help='''Simulation prefixes within the current directory. Default
        behaviour (without any prefix specified) is to get the prefix from
        ./batch.pfx. Otherwise, the file list is generated as the merged list from
        multiple simulation prefixes in this directory.'''
        )
    args = ap.parse_args()

    logger.setLevel(logging.INFO + (args.verbose-1)) # Default is INFO+(1-0)

    # check for keep file
    if args.keepfile_name and not os.path.isfile(args.keepfile_name):
        ap.error(f'Could not find keep file {args.keepfile_name}')

    keep_cats, keep_files = BaseRunner.read_keep_file(args.keepfile_name)
    if args.keep_inputs:
        keep_cats.append('input')
    else:
        keep_cats = [ c for c in keep_cats if c != 'input' ]

    runners = []
    for p in args.pfx:
        _runner = BaseRunner(None, simdir=p)
        runners.append(_runner)
        if _runner.simdir != '.':
            ap.error(f'Prefix {p} points outside the current directory.')
    runners.sort(key=_key_by_lst_mtime, reverse=True) # newest first

    logger.debug('Processing order:\n'+'\n'.join(f'{i:2} {r.prefix}' for i,r in zip(count(1),runners)))

    # use a single list of files where the each prefix below can
    # 'claim'/'categorize' its files and pass back those that have not been
    # categorized
    files = set(Path('.').iterdir())
    catfiles = dict()

    # use the runners in reverse-order of invocation categorize the files.
    for runner in runners:
        runner.preexistingfiles = []
        _catfiles, files = runner.categorize_files(files)
        catfiles = ravel_merge_nested(_catfiles, catfiles)

    # Determine whether each file will be deleted. Create a dict of bool,
    # whether the file (key) will be kept (value)
    keepmaskfiles = dict()
    for fn,cat in catfiles.items():
        _k = False # default, delete
        _k |= any(fnmatch(cat, c) for c in keep_cats)
        _k |= any(fnmatch(fn, g) for g in keep_files)
        keepmaskfiles[fn] = _k

    # print the listing
    if args.dry_run:
        s = ''
        cw = max(map(len, catfiles.values()))
        lastcat = ''
        for fn,cat in sorted(catfiles.items(), key=itemgetter(1)):
            keep=keepmaskfiles[fn]
            if cat == lastcat:
                cat = '' # avoid redundant output
            else:
                lastcat = cat
            s += f'{"keep" if keep else "delete":6} {cat:{cw}} {fn!s}\n'

        logger.info(f'Action {"Category":{cw}} File\n'+s)
        exit(0)

    # do the deletion
    for fn,keep in filterfalse(itemgetter(1), keepmaskfiles.items()):
        logger.info(f'Deleting {fn}')
        os.remove(fn)

