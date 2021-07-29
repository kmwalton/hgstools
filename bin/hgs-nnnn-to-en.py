#!/usr/bin/env python
"""Convert a .nnnn-suffix file to a .Xen file so it may be used as a restart.

e.g. Convert an arbitrary concentration file to a concentration restart file:
    hgs-nnnn-to-en.py prefixo.conc_pm.species.0042 prefixo.restart.cen
"""

import sys
import argparse

argp = argparse.ArgumentParser()
argp.add_argument('FILE_IN',
    help='Input file',
    )
argp.add_argument('FILE_OUT',
    help='Output file',
    )
args = argp.parse_args()

with open(args.FILE_OUT, 'wb') as fout:
    with open(args.FILE_IN, 'rb') as fin:
        fin.seek(4)
        time = bytes(fin.read(80)).decode('UTF-8').strip()
        fin.seek(88)
        fout.write(fin.read())

print(f'{args.FILE_IN} at t={time} converted successfully.')
sys.exit(0)
