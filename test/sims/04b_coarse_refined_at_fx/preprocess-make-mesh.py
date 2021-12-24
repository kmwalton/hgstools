#!/usr/bin/env python

import subprocess

if __name__ == '__main__':
   
    cp = subprocess.run(
        [ 'ofrac2hgs.py',
            '--grid-out', 'gridinclude.grok',
            '--frac-out', 'fracinclude.grok',
            '--refine-near-fx-plane', '0.05', '0.10',
            '--max-grid-space', '2.5',
            '--', # end the list; begin positional arguments
            'fracs.rfd',
        ],
        )
    cp.check_returncode()
