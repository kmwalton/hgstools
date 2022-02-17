#!/usr/bin/env python
"""Parses information out of a grok file and makes a dict

Limitations:
VERY limited --- will only work on a few keywords

"""

import sys
import os
import re
import io
import pprint
import argparse
import pickle
from more_itertools import grouper

_NUM_RE = r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?'

_tags = {
    'files_mprops':(lambda txt:
        list( re.findall(
            r'^use domain type\nporous media.*?properties file\n(.*?)\n',
            txt, flags=re.S|re.M|re.I)
        )),
    'files_fprops':(lambda txt:
        list( re.findall(
            r'^use domain type\nfracture.*?properties file\n(.*?)\n',
            txt, flags=re.S|re.M|re.I)
        )),
    'files_include':(lambda txt:
        set( re.findall(
            r'^include\s+(.*?)\n',
            txt, flags=re.M|re.I)
        )),
    'zoned distribution coefficient':(lambda txt:
        list(map(float,
            re.search(f'zoned distribution coefficient((?:\\s{_NUM_RE})+)',
                txt, flags=re.M|re.I).group(1).strip().split())
        )),

}

def _get_as_stream(file, mode='r'):
    if type(file) is str:
        return open(file,mode)
    elif isinstance(file,io.IOBase):
        return file
    else:
        raise ValueError(f'Unhandled type {type(file)}')

def _clean(text_in):
    """Clean out hgs comments and extra whitespace"""

    # get rid of comments
    text = re.sub('\!.*?\n','\n',text_in)

    # get rid of trailing whitespace, leading space, multiple newlines
    text = re.sub(r'\s*\n\s*','\n',text)

    return text


def parse(file, do_includes=False, do_skips=False):
    """Return a dictionary of grok file data
    """

    d = {}

    file = _get_as_stream(file)
    file_name = file.name # hope the underlying stream has this attribute
    file_text = file.read()
    file_text = _clean(file_text)

    original_includes = _tags['files_include'](file_text)
    if do_includes and original_includes:
        # use the 'include files' by appending them to the text
        for f in original_includes:
            # resolve file name
            if not os.path.isabs(f):
                f = os.path.join(os.path.dirname(file_name), f)

            with open(f,'r') as fin:
                file_text += fin.read()

        # clean again
        file_text = _clean(file_text)

    # parse-out skip on/skip off lines so non-standard definitions prevail
    if do_skips:
         file_text = re.sub(
                    r'\s*skip\s+on.*?skip\s+off\s*', '\n',
                    file_text,flags=re.DOTALL)
    else:
         file_text = re.sub(r'\s*skip\s+(on|off)\s*','\n',file_text)

    # make final dictionary
    for item,itemfinder in _tags.items():
        try:
            d[item] = itemfinder(file_text)
        except (AttributeError, ValueError) as e:
            pass

    d['files_include'] = original_includes

    return d

if __name__=='__main__':
    
    argp = argparse.ArgumentParser()
    argp.add_argument('-es', action='store_true',default=False,
        help='Enforce skip on/off in the file. Default: keep skipped text.')
    argp.add_argument('-ii', action='store_true',default=False,
        help='Include "include" files. Default: do not include these.')
    argp.add_argument('FILE',type=argparse.FileType('r'),)
    argp.add_argument('PICKLE_FILE_OUT', nargs='?',
        default=None,
        type=argparse.FileType('wb'),
        help='Optional output file name for pickled results.')
    args = argp.parse_args()

    # work
    grokdata = parse(args.FILE, do_skips=args.es, do_includes=args.ii)

    # output
    fn = os.path.basename(args.FILE.name)
    print(f'Data found in {fn}')
    pprint.pp(grokdata)

    sys.exit(0)
