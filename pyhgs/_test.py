"""Supports testing of pyhgs by unittest and/or with HGS simulation data"""

import os

def skip_if_no_sim_output(prefix, suffixes=['o.eco','o.lst','o.hsplot.eco',]):
    """Return True if any files with given prefix+suffixes are missing"""
    for s in suffixes:
        if not os.path.isfile(f'{prefix}{s}'):
            return True
    return False


