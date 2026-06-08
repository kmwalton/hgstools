"""Supports testing of pyhgs by unittest and/or with HGS simulation data"""

__docformat__ = 'numpy'

import os

def skip_if_no_sim_output(prefix, suffixes=['o.eco','o.lst',]):
    """Return True if any files with given prefix+suffixes are missing"""
    for s in suffixes:
        if not os.path.isfile(f'{prefix}{s}'):
            return True
    return False


