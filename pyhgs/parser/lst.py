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
    'timestep_start':re.compile(
        r'^-{82}\n\s+SOLUTION FOR TIMESTEP\s+(\d+)',
        flags=re.M),
    'timestepsize':re.compile(
        r'(?ms:^Global target time.*Tnext)' \
        r'(^(?:\s+'+_NUM_RE+'){3}(?:\s+.*)$)+'),
    'accepted_solution_time':re.compile(
        r'^ Accepted solution at time:?\s+('+_NUM_RE+')',
        flags=re.M),
    'initial_time':re.compile(
        r'Initial time =\s+('+_NUM_RE+')',
        flags=re.M),
    'met_global_target':re.compile(
        r'^ Met global target time\ +('+_NUM_RE+')',
        flags=re.M),
    'simulation_time_report':re.compile(
        r'^-{9} SIMULATION TIME REPORT\s*\n(?P<report>(?ms:.*?))\n-{58}',
        flags=re.M),

    'mass_balance': re.compile(
        r'^RATE OF MASS EXCHANGE\s+IN\s+OUT\s+TOTAL.*\n' \
        r'(?s:(?P<bc_data>.*?))\n' \
        r'   (?P<tot>NET1 EXCHANGE RATE.*?)\s+(?P<totval>(?:'+_NUM_RE+'))\s*$',
        flags=re.M),

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
        self._n_ts = -1
        self._has_sim_report = False

        _ts = {}
        for m in FILTER['timestep_start'].finditer(self._txt):
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

        self._n_ts = len(_ts)

        # bookend the list of timesteps with 0 for the beginning of the file,
        # and the simulation report + the max length
        self._tsloc = [0,] + list( v for k,v in sorted(_ts.items()) )

        # find the simulation report
        sim_rpt = FILTER['simulation_time_report'].search(
                self._txt[self._tsloc[-1]:])
        if sim_rpt:
            self._tsloc.append(self._tsloc[-1]+sim_rpt.span()[0])
            self._has_sim_report = True

        # end of input
        self._tsloc.append(len(self._txt))

        logger.debug(f'Found {len(self._tsloc)} timesteps.')

    def _iter_timestep_text_bounds(self, itimestep=None):
        """Generate the start and end character numbers of timestep text blocks

        Arguments
            itimestep : int or None
                Get an iterator over only the requested timestep. Default None,
                implying all timesteps

        Yields tuples of (timestep index, text start position, text end
        position), where positions are absolute character indices in self._txt.
        """

        itss = self._tsloc[1:-1-self._has_sim_report]
        itse = self._tsloc[2:]

        if itimestep is not None:
            itss = [self._tsloc[itimestep],]
            itse = [self._tsloc[itimestep+1],]

        for i, tss, tse in zip(count(1), itss, itse):
            yield (i,tss,tse,)

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

    def get_global_target_times(self):
        """Return a list of ( target time, itime when target achieved)-tuples.
        """
        ret = []

        for i, tss, tse in self._iter_timestep_text_bounds():
            # look in last 256 chars of a timestep for this regex
            m = FILTER['met_global_target'].search(self._txt[tse-256:tse])

            if m:
                ret.append((float(m.group(1)), i,))

        return ret

    def get_n_ts(self):
        """Return the number of timesteps

        For parameters to other functions in this module, valid timestep indices
        are 1 to this number, inclusive. (Hydrogeosphere reports 1-based time
        indices.)
        """
        return self._n_ts

    def get_ts_time(self, itime='all'):
        """Get the simulation time at the end of the itime'th timestep.

        Arguments:
            itime : int or 'all'
                A specific timestep index or the keyword 'all' to signify
                all simulation solution times are desired.

        Returns
            A singleton float value.
            Or, a generator of simulation times over all timesteps. (This
            includes the initial time of the simulation in timestep 0, so the
            length of this generated sequence is get_n_ts()+1.)
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

                for i, tss, tse in self._iter_timestep_text_bounds():
                    match = FILTER['accepted_solution_time'].search(
                            self._txt[tss:tse])
                    if not match:
                        raise RuntimeError(
                            f'No accepted solution found in timestep {i}')
                    yield float(match.group(1))

            return _gen_times()

    def get_ts_dtime(self, itime='all'):
        """Return the time delta of this step

        Arguments:
            itime : int or 'all'
                A specific timestep index or the keyword 'all' to signify
                all simulation solution times are desired.

        Returns
            A singleton float value.
            Or, a generator of time deltas over all timesteps. (The
            initial time delta in this sequence is the time between the 'Initial
            time' of the simulation and the result of `get_ts_time(1)`.
        """

        if itime <= 0:
            raise ValueError('itime must be >= 1')

        elif not itime == 'all':
            return float(self.get_ts_time(itime)-self.get_ts_time(itime-1))

        else:

            def _gen_dt():

                itime = self.get_ts_time('all')
                itm1 = next(itime)

                for it in itime:
                    yield(it-itm1)
                    itm1 = it

            return _gen_dt()


    def get_ec(self):
        """Return 0 for NORMAL EXIT, 1 otherwise."""

        if self._has_sim_report and \
                re.search(r'- +NORMAL EXIT +-', self._txt[self._tsloc[-2]:],
                flags=re.M ):
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

        for itime,tss,tse in self._iter_timestep_text_bounds(itimestep):
            for err_re in ERRORS_FILTER:
                m = err_re.search(self._txt,tss,tse)
                if m:
                    yield (itime, m.group('name'), m.group('message'),)


    _fluid_balance_regex = [
        # revision ~1827 steady-state simulation
        re.compile(
            r'RATE OF FLUID EXCHANGE\s+IN\s+OUT\s+TOTAL.*\n' \
            r'(?s:(?P<bc_data>.*))\n'\
            r'TOTAL\S+(?P<total>(?:\s+'+_NUM_RE+'){3}).*\n',
            flags=re.M),

        # revision ~2469 transient flow simulation
        re.compile(
            r'RATE OF FLUID EXCHANGE.*\n' \
            r'Boundary condition name\s+IN\s+OUT\s+NET.*' \
            #r'(?P<bc_data>(?:\n\S.+\s+(?:+'+_NUM_RE+'\s+){3}\S.*){1,})'\
            r'(?s:(?P<bc_data>.*?))\n'\
            r'TOTAL:\[[^\]]+\]\s+(?P<total>(?:\s+'+_NUM_RE+'){3})\s*\n',
            flags=re.M)
    ]

    def get_fluid_balance(self, itimestep=1):
        """Return a dictionary of { BC_name:(Q_in, Q_out, Q_net, Domain), .. }
        """


        # some steady state simulations have no timesteps
        if self.get_n_ts() == 0:
            itimestep = 0

        m = None
        for fbregex in LSTFileParser._fluid_balance_regex:
            m = list( fbregex.finditer(self._txt,
                    self._tsloc[itimestep],
                    self._tsloc[itimestep+1]) )

            if m:
                # Return the final fluid balance block within the timestep.
                # There may be more than one such block if the flow or
                # transport solution breaks the "multiplier" criterion and the
                # step is repeated with a reduced step size
                m = m[-1]
                break

        if not m:
            raise RuntimeError(f'No fluid balance match found at ts={itimestep}')

        ret = {}

        # get BCs
        for l in m.group(1).strip().split('\n'):
            l = l.strip().split(maxsplit=4)
            ret[l[0]] = (float(l[1]), float(l[2]), float(l[3]), l[4],)

        # get total
        l = m.group(2).strip().split()
        ret['TOTAL'] = (float(l[0]), float(l[1]), float(l[2]), None,)

        return ret

    def get_mass_balance(self, itimestep=1):
        """Return a dictionary of { BC_name:(Q_in, Q_out, Q_net, Domain), .. }
        """

        m = list( FILTER['mass_balance'].finditer(self._txt,
                self._tsloc[itimestep],
                self._tsloc[itimestep+1]) )

        if not m:
            raise RuntimeError(f'No fluid balance match found at ts={itimestep}')

        # Return the final fluid balance block within the timestep.
        # There may be more than one such block if the flow or
        # transport solution breaks the "multiplier" criterion and the
        # step is repeated with a reduced step size
        m = m[-1]

        ret = {}

        # get BCs
        for l in m['bc_data'].strip().split('\n'):
            mm = re.match(r'\s*(.*?)'+2*('\s+('+_NUM_RE+')') +'\s*', l)
            ret[mm[1]] = (float(mm[2]), float(mm[3]),)

        # get total
        ret[m['tot']] = float(m['totval'])

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
