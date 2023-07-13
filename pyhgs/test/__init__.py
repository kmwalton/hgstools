import os
"""Module to provide access to `test` directory adjacent to `pyhgs`"""

_TESTPATH = os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..', '..', 'test'))
'Path to the test directory beside pyhgs'

def get_sims_dir():
    """Return an absolute path to the test_sims directory"""
    return _TESTPATH+os.sep+'test_sims'

def sims_join(*args):
    """Equivalent to `os.path.join(get_sims_dir(), *args)`

    A convenience method to get /path/to/test_sims/sim_dir/prefix
    """
    return get_sims_dir() + os.sep + os.path.join(*args)
