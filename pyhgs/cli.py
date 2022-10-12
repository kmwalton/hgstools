"""Miscellaneous tools to useful in command line interfaces (CLI)."""
import argparse
import os

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

        path_to, prefix = os.path.split(args.PATH_TO_PREFIX)
        
    Alternatively, the static method `parse_path_to_prefix` may be used to
    extract the same as above from any string (and outside the task of parsing
    command line arguments).
    
        path_to, prefix = PathToPrefix.parse_path_to_prefix('path/to/pfx')
    """

    def __init__(self, *args, **kwargs):
        """Pass all arguments to `Action.__init__`"""
        super().__init__(*args, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        (pth, pfx) = PathToPrefix.parse_path_to_prefix(values)

        # TODO What if nargs != 1? This should be able to return either a string
        # value or a list, per the convention of nargs

        setattr(namespace, self.dest, os.path.join(pth,pfx))


    @staticmethod
    def parse_path_to_prefix(s):
        """Determine the path_to and prefix for any given string.

        Assume the current directory if no path is apparent in `s`.
        Look for *batch.pfx* if no problem prefix is apparent in `s`.
        """

        (pth, pfx) = (None,None)

        if os.path.isdir(s):
            # look for '.' or '..', which will get path.split into the tail
            (pth, pfx) = (s, None)
        else:
            (pth, pfx) = os.path.split(s)

        if not pth:
            pth = '.'

        if not pfx:
            batchfn = pth+os.path.sep+'batch.pfx'
            try:
                with open(batchfn,'r') as fin:
                    pfx = fin.read().strip()
            except FileNotFoundError:
                parser.error( f'Could not find {batchfn}; '\
                    'cannot auto-detect problem prefix')

        return (pth, pfx)
