"""Utilities for testing pyhgs with HGS simulation data/outputs"""

import os

def skip_if_no_sim_output(prefix, suffixes=['o.eco','o.lst','o.hsplot.eco',]):
    """Return false if any files with the given suffixes are missing"""
    for s in suffixes:
        if not os.path.isfile(f'{prefix}{s}'):
            return True
    return False

