#!/usr/bin/env python

import os
import re

from itertools import count

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

ERRORS_FILTER = [
    re.compile(r'^(?P<name>W0SOLV): (?P<message>.*)$', flags=re.M),
]
"""A list of regular expressions for error messages found in the .lst output
"""

_NUM_RE = r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?'
"""RegEx for a floating point number"""

FILTER = {
    'timestepsize':re.compile(
        r'(?ms:^Global target time.*Tnext)' \
        r'(^(?:\s+'+_NUM_RE+'){3}(?:\s+.*)$)+'),
    'accepted_solution_time':re.compile(
        r'^ Accepted solution at time\s+('+_NUM_RE+')', flags=re.M),
    'initial_time':re.compile(r'Initial time =\s+('+_NUM_RE+')', flags=re.M),

}
"""Dictionary of compiled re objects for various chunks of a .lst file"""


class LSTFileParser:
    """Parse and return (select) data from phgs output, the .lst file"""

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
        # HGS labels its timesteps starting from 1.
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

        # bookend the list of timesteps with 0 for the beginning of the file and
        # the max length
        self._tsloc = [0,] \
                  + list( v for k,v in sorted(_ts.items()) ) \
                  + [len(self._txt),]
        logger.debug(f'Found {len(self._tsloc)} timesteps.')

    def ss_flow(self):
        """Return True if the flow system is steady state"""

        # 'steady-state' statements should be found in the preamble (failure
        # case) or before the beginning of the second timestep
        _end = self._tsloc[-1]
        if self.get_n_ts() > 1:
            _end = self._tsloc[2]

        m = re.search('Steady-state (?:flow )?simulation',
            self._txt[:_end],
            flags = re.M|re.S)

        return bool(m)

    def get_n_ts(self):
        """Return the number of timesteps

        Valid timestep indices are 1 to this number, inclusive.
        """
        # subtract 0 from the "0" timestep, which is all the text before the
        # first timestep
        return len(self._tsloc)-2

    def get_ts_time(self, itime='all'):
        """Get the simulation time, as of the end of the itime'th timestep.

        Arguments:
            itime : int or 'all'
                Return the simulation time (in simulation time units) of the
                specified timestep (int). Or, If itime is the keyword 'all',
                return a generator of all simulation time for all timesteps.

        Returns
            A singleton float value if the time of one timestep is requested
            Or, a generator of all times. This includes the initial time of the
            simulation in timestep zero.

        This probably won't work for restarted simulations, 
        """

        if not itime == 'all':
            # allow ValueError to be raised if necessary
            its = int(itime)
            # allow bad index errors to be raised
            tss = self._tsloc[its]
            tse = self._tsloc[its+1]

            if its == 0:
                match = FILTER['initial_time'].search(
                        self._txt[tss:tse])
            else:
                match = FILTER['accepted_solution_time'].search(
                        self._txt[tss:tse])

            return float(match.group(1))

        else:
            def _gen_times():

                match = FILTER['initial_time'].search(
                        self._txt[:self._tsloc[1]])
                yield float(match.group(1))

                itss = self._tsloc[1:-1]
                itse = self._tsloc[2:]
                for i, tss, tse in zip(count(1), itss, itse):
                    match = FILTER['accepted_solution_time'].search(
                            self._txt[tss:tse])
                    if not match:
                        raise RuntimeError(
                            f'No accepted solution found in timestep {i}')
                    yield float(match.group(1))

            return _gen_times()


    def get_ec(self):
        """Return 0 for NORMAL EXIT, 1 otherwise."""

        exit_phrase = re.search(r'SIMULATION TIME REPORT.*---- (.*?) ----',
            self._txt[self._tsloc[-2]:],
            flags = re.M | re.S )

        if exit_phrase and exit_phrase.group(1) == 'NORMAL EXIT':
            return 0

        return 1

    def iter_errors(self,itimestep=None):
        """Yield errors (timestep index, error name, error message)

        Arguments
        ---------

        itimestep : int, optional
            Yield errors for the given timestep, only. Default: yield errors for
            all timesteps.
        """

        itss = self._tsloc[:-1]
        itse = self._tsloc[1:]

        if itimestep is not None:
            itss = [self._tsloc[itimestep],]
            itse = [self._tsloc[itimestep+1],]


        for itime,tss,tse in zip(count(), itss, itse):
            for err_re in ERRORS_FILTER:
                m = err_re.search(self._txt,tss,tse)
                if m:
                    yield (itime, m.group('name'), m.group('message'),)


    def get_fluid_balance(self, itimestep=1):
        """Return a dictionary of { BC_name:(Q_in, Q_out, Q_net, Domain), .. }
        """

        fbre = re.compile(
            r'RATE OF FLUID EXCHANGE\s+IN\s+OUT\s+TOTAL.*\n' \
            r'(?s:(?P<bc_data>.*))\n'\
            r'TOTAL\S+(?P<total>(?:\s+'+_NUM_RE+'){3}).*\n',
            flags=re.M)

        m = fbre.search(self._txt,
                self._tsloc[itimestep],
                self._tsloc[itimestep+1])

        ret = {}

        # get BCs
        for l in m.group(1).split('\n'):
            l = l.strip().split()
            ret[l[0]] = (float(l[1]), float(l[2]), float(l[3]), l[4],)

        # get total
        l = m.group(2).strip().split()
        ret['TOTAL'] = (float(l[0]), float(l[1]), float(l[2]), None,)

        return ret


    def get_flow_solver_iterations(self, itimestep=1):

        # for flow solver iterations, search the preamble as well as timestep 1
        # if timestep 1 is requested.

        _s = self._tsloc[0]
        if itimestep > 1:
            _s = self._tsloc[itimestep]

        _e = self._tsloc[-1]
        if itimestep <= self.get_n_ts():
            _e = self._tsloc[itimestep+1]

        m = re.search(
            r'^\s*Number of flow matrix solver iterations = +(\d+)',
            self._txt[_s:_e],
            flags=re.M)

        if m:
            return int(m.group(1))
        return 0

    def get_transport_solver_iterations(self, itimestep=1):

        _s = self._tsloc[0]
        if itimestep > 1:
            _s = self._tsloc[itimestep]

        _e = self._tsloc[-1]
        if itimestep <= self.get_n_ts():
            _e = self._tsloc[itimestep+1]

        m = re.search(
            r'^\s*Number of transport matrix solver iterations: +(\d+)',
            self._txt[_s:_e],
            flags=re.M)

        if m:
            return int(m.group(1))
        return 0
