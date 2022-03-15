#!/usr/bin/env python

import os
import re

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

ERRORS_FILTER = [
]

class LSTFileParser:
    """Parse and return data from phgs output, the .lst file"""

    def __init__(self, fnin):
        """Load the .lst file text, given file name or HGS prefix"""


        # Resolve the file name
        # in case a prefix was given
        if not fnin.endswith('o.lst'):
            logger.debug(f'Adding suffix o.lst to {fnin}')
            fnin += 'o.lst'
        self._fnin = fnin

        # checks
        if not os.access(fnin, os.R_OK):
            raise RuntimeError(
                f'Could not find/read {fnin} in {os.getcwd()}.')
        if os.path.getsize(fnin) > 100 * (1 << 20):
            raise RuntimeError('File is > 100 MiB (arbitrary). Not parsing.')

        logger.debug(f'Accepting read-ability and size of {fnin}')

        # injest text
        with open( fnin, 'r' ) as fin:
            self._txt = fin.read()


        # index the timestep locations
        self._tsloc = []
        _ts = {}

        for m in re.finditer(
                r'^ +SOLUTION FOR TIMESTEP +(\d+)',
                self._txt,
                flags=re.M, ):
            itime = int(m.group(1))

            if itime in _ts:
                if isinstance(_ts[itime],int):
                    _ts[itime] = [ _ts[itime], m.span()[0], ]
                else:
                    _tsloc.append(m.span()[0])
            else:
                _ts[itime] = m.span()[0]

        if any( isinstance(v,list) for v in _ts.values() ):
            raise NotImplementedError('Found timestep reported twice')

        self._tsloc = list( v for k,v in sorted(_ts.items()) )
        logger.debug(f'Found {len(self._tsloc)} timesteps.')

    def ss_flow(self):
        """Return True if the flow system is steady state"""
        m = re.search('Steady-state simulation(No accumulation in domain)',
            self._txt[self._tsloc[0]:self._tsloc[1]],
            flags = re.M|re.S)

        return bool(m)

    def get_n_ts(self):
        """Return the number of timesteps"""
        return len(self._tsloc)

    def get_ec(self):
        """Return 0 for NORMAL EXIT, 1 otherwise."""

        exit_phrase = re.search(r'SIMULATION TIME REPORT.*---- (.*?) ----',
            self._txt[self._tsloc[-1]:],
            flags = re.M | re.S )

        if exit_phrase.group(1) == 'NORMAL EXIT':
            return 0

        return 1

    #def iter_errors(self):

    #    for tss,tse in zip(


