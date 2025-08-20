"""Miscellaneous tools to useful in command line interfaces (CLI)."""
import argparse
import os
import re
import warnings

__docformat__ = 'numpy'

def parse_path_to_prefix(s='.'):
    """Determine the `(path_to, prefix, leftovers)` for any given string.

        path_to, prefix, leftovers =
            pyhgs.cli.parse_path_to_prefix(args.PATH_TO_PREFIX)

    `path_to`, `prefix` and `leftovers` are determined by the following
    ways:

    1. `s` is examined for drive and directory components. Any components that
       are found will be returned in `path_to`, or `path_to` will be '.'.
    2. If `path_to/batch.pfx` exists, then the contents of *batch.pfx* will
       be returned in `prefix` and the balance of `s` in `leftovers`
    3. if `o.` or `.` in the non-`path_to` part of `s`, then `prefix` and and
       `leftovers` take the strings split before this marker.
    4. if the non-directory part of `s` is a single word, that word is
       returned in `prefix` and `leftovers` is the empty string.

    Parameters
    ----------
    s : str or Path-like
        The string or Path to examine. Default '.', the the current working
        directory.
    """

    (pth, pfx, more) = ('', '', '')

    if os.path.isdir(s):
        # look for '.' or '..', which will get path.split into the tail
        pth = s
    else:
        (pth, more) = os.path.split(s)

    if not pth:
        pth = '.'

    batchpfx = ''
    batchfn = pth+os.path.sep+'batch.pfx'
    try:
        with open(batchfn,'r') as fin:
            batchpfx = fin.read().strip()
    except FileNotFoundError:
        # should this halt the program?!
        #parser.error( f'Could not find {batchfn}; '\
        #        'cannot auto-detect problem prefix')
        pass

    if batchpfx and more and more.startswith(batchpfx):
        pfx = batchpfx
        more = more[len(batchpfx):]
    elif batchpfx and not more:
        pfx = batchpfx
        more = ''
    #elif not batchpfx:
    else:
        pats = [
            # Probable hgs output file: break at 'o.'
            r'(.*?)(o\..*)',
            # filename with a dot extension: find nothing as the prefix
            # and find everything as 'more'
            r'()(.*\..*)',
            # No dot extension: find everything as the prefix
            r'([^.]*)()',
        ]

        for p in pats:
            m = re.match(p, more)
            if m:
                pfx, more = m.groups()
                break

    return (pth, pfx, more)



class PathToPrefix(argparse.Action):
    """Custom `argparse.Action` to find `PATH_TO_PREFIX` for a HGS problem

    The `__call__` method of this Action will take any string its `value` input,
    including the empty string, a directory name (relative or absolute), or the
    combination of a directory name and an HGS problem prefix.

    It defines, in 'dest' in the ArgumentParser object's namespace, a relative
    (or absolute) path and the problem prefix. This can be split using
    `os.path.split` as in the example below.

    If it does not appear as part of the argument's `value`, `PREFIX` is
    determined by finding and reading a *batch.pfx* HGS file in the current
    working directory.

    This action is best used with a positional variable, like:

        argp.add_argument('PATH_TO_PREFIX',
            metavar='path/to/prefix',
            action=pyhgs.cli.PathToPrefix,
            help='path/to/prefix of the HGS simulation',
            )

        args = argp.parse_args()

        path_to, filename = os.path.split(args.PATH_TO_PREFIX)
        path_to, prefix, leftovers = PathToPrefix.split(args.PATH_TO_PREFIX)
        
    """

    def __init__(self, *args, **kwargs):
        """Pass all arguments to `Action.__init__`"""
        super().__init__(*args, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        (pth, pfx, more) = parse_path_to_prefix(values)

        # TODO What if nargs != 1? This should be able to return either a string
        # value or a list, per the convention of nargs

        setattr(namespace, self.dest, os.path.join(pth,pfx+more))


    @staticmethod
    def parse_path_to_prefix(s):
        '''Alias for `pyhgs.cli.parse_path_to_prefix`'''
        warnings.warn('Use pyhgs.cli.parse_path_to_prefix', DeprecationWarning)
        return parse_path_to_prefix(s)


    @staticmethod
    def split(s):
        '''Alias for `pyhgs.cli.parse_path_to_prefix`'''
        return parse_path_to_prefix(s)
