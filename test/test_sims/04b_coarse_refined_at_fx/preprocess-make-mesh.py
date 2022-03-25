#!/usr/bin/env python

import os
import shutil
import subprocess

if __name__ == '__main__':

    child_env = os.environ.copy()
   
    cp = subprocess.run(
        [ 'python.exe', shutil.which('ofrac2hgs.py'),
            '--grid-out', 'gridinclude.grok',
            '--frac-out', 'fracinclude.grok',
            '--refine-near-fx-plane', '0.05', '0.10',
            '--max-grid-space', '2.5',
            '--', # end the list; begin positional arguments
            'fracs.rfd',
        ],
        env=child_env,
        )
    cp.check_returncode()
