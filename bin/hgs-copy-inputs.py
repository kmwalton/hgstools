#!/usr/bin/env python
"""Copy an HGS simulation's input files to the specified destination directory

Copy all the normally-named HGS input files, plus aptly named pre- and
post-processing files, plus any other explictly-named, glob-, or regex-matched
files from the current directory to another directory.

Includes:
*.grok
*.[mfdw]props
parallelindex.dat
array_sizes.default
*.control
batch.pfx

Other includes:
preprocess*
postprocess*

Excludes:
File names beginning with a '.' (Unix hidden files)
File names ending with a '~' (Some file editor temporary files)
"""

import sys
import os
import re
import argparse
import shutil
from glob import glob
from itertools import chain
from math import log10

#helpers
def matchesAny( s, patterns, exacts ):
    """Returns True if s matches anything in patterns or exacts"""
    if s in exacts:
        return True
    for p in patterns:
        if p.match(s):
            return True
    return False

def dummyCopy(f,dest):
    """Return True if f exists and dest is writeable."""
    return os.path.isfile(f) and os.access(dest,os.W_OK)

# set up argument parser
p = argparse.ArgumentParser( description=__doc__.split('\n')[0],
       epilog='\n'.join(__doc__.split('\n')[1:]),
       formatter_class=argparse.RawDescriptionHelpFormatter)

vgrp = p.add_mutually_exclusive_group()
vgrp.add_argument('-v','--verbose',
        default=1,
        action='count',
        help='Increase verbosity (up to 2). Default: verbose=1.',
        )
vgrp.add_argument('-q','--quiet',
        action='store_const',
        const=0,
        dest='verbose',
        help='Silence all output. (Sets verbose=0).',
        )

p.add_argument( '--dry-run',
        action='store_true',
        default=False,
        help='Do a dry-run; do not copy anything. Default: Non-dry-run; copy.',
        )

p.add_argument('-e','--regexp',
        metavar='REGEXP',
        dest='include_regex',
        action='append',
        default=[],
        help='Additional regular expressions for matching file names to copy.',
        )

#positional parameters
p.add_argument('globs',
        nargs='*',
        metavar='FILE | GLOB',
        default=[],
        help='File names or glob-style strings specifying extra files to copy',
        )
p.add_argument('dest', type=str, help='The destination directory')

args = p.parse_args()

# set up exclusion patterns
exclude_regex = [r'^\..*', f'.*\~$', ]
excludePatterns = list( map( lambda s: re.compile(s,re.IGNORECASE), exclude_regex ) )
del exclude_regex

# set up inclusion patterns / files / globs
include_regex = [
    r'.*\.grok',
    r'.*\.[mfdw]props',
    r'parallelindx\.dat', 
    r'array_sizes\.default',
    r'.*\.control',
    r'batch\.pfx',
    r'preprocess.*',
    r'postprocess.*'
]

includePatterns = list(
    map(
        lambda s: re.compile(s,re.IGNORECASE),
        include_regex
        + args.include_regex
    ) )
del include_regex

includeFiles = list( chain(*map(glob, args.globs)))

# make list of files to be copied
fileList = list(
    filter(
        lambda f:
            matchesAny(f,includePatterns,includeFiles)
            and not matchesAny(f, excludePatterns, []),
        os.listdir('.') ) )


#
# do some work
#

# set up copy function and ensure destination exists
copy_func = None
if args.dry_run:
    # set dummy copy function
    copy_func = dummyCopy
else:
    # set real copy function
    copy_func = shutil.copy2
    # make destination directory
    os.makedirs( args.dest, exist_ok=True )


# do copying
for f in fileList:
    if copy_func( f, args.dest ):
        if args.verbose > 1:
            print( "Copied {}".format(f) )
if args.verbose:
    print(f'Copied {len(fileList)} files to {args.dest}.')

# exit with failure only if something above results in an Exception
sys.exit(0)
