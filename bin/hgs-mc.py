#!/usr/bin/env python
"""Run Monte Carlo-style realization of a Hydrogeosphere problem

Examples:
    Run 5 variants on ../base_sim/. Creates ./mc1, ..., ./mc5.
    $ python hgs-mc.py 5 ../base_sim

Dependencies:
    - hgs-copy-inputs.py (for copying base_sim)
    - hgs-runall.ps1 (for running an instance)
    - pre- and post-processing scripts in base_sim to set up a new realization

"""
import sys
import os
import argparse
import shutil
import shlex
import time
import subprocess
import datetime
import tempfile
import enum
import re
from multiprocessing import Pool
from multiprocessing.shared_memory import ShareableList
from math import ceil,log10

from pyhgs.runner import PyPowerShellRunner

import logging


def _logger_setup(l):
    l.verbose1 = lambda msg : l.log(logging.INFO-1,msg)
    l.verbose2 = lambda msg : l.log(logging.INFO-2,msg)
    l.verbose3 = lambda msg : l.log(logging.INFO-3,msg)

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
_logger_setup(logger)

# Useful:
#https://wiseodd.github.io/techblog/2016/06/13/parallel-monte-carlo/

class InstanceStatusEnum(enum.IntEnum):
    NOT_STARTED = enum.auto()
    SET_UP_RUNNING = enum.auto()
    SET_UP_FAILED = enum.auto()
    SET_UP_SUCCESS = enum.auto()
    RUNNING = enum.auto()
    FINISHED_FAILED = enum.auto()
    FINISHED_SUCCESS = enum.auto()

class HGS_MCRunner():
    """Manage and run a batch of HGS instances, Monte Carlo-style."""

    def __init__(self, copy_command, tc_command, base_dir, keep_file, stop,
            start=0,):
        """Initialize the runner and make a temporary copy of base_dir.

        This run operates out of a copy of base_dir to avoid errors caused by
        the user "changing something" during the MC run.

        Arguments:
            copy_command : str
                A string representing the shell command and its arguments that
                will be used to copy base_sim to its MC instance directory

            tc_command : str
                A string representing the shell command and its arguments that
                will be used to run the Hydrogeosphere tool chain

            base_dir : str
                A a string representing the path to the base simulation instance

            keep_file : str or None
                A path to a file with a list of filenames or globs representing
                files to retain in MC instance directories after the instance is
                complete

            stop : int
                An integer representing the stopping index value. (Usually this
                is the number of instances; the final index will be (stop-1).)

            start : int
                Optional. Default zero. The starting instance index.
        """

        self.RUN_START = datetime.datetime.now()
        """Time of object creation. Assuming that one HGS_MCRunner object is
        made per monte carlo run, this timestamp can be a unique identifier for
        this run."""

        self.DRY_RUN = False
        """Skip running Hydrogeosphere, hgs-runall. This may (will) have have
        consequences in subsequent analysis steps, as there may not (will not)
        be any output to analyze."""

        self._start = start
        self._stop = stop

        self._inst_status = ShareableList(
                self._stop*[InstanceStatusEnum.NOT_STARTED.value,])
        """Status of each instance"""

        self._inst_ec = ShareableList(self._stop*[-1,])
        """Exit code of each instance"""

        self.copy_command = shlex.split(copy_command)
        """Program and arguments needed to copy the simulation inputs"""
        if self.copy_command[0].lower().endswith('.py'):
            self.copy_command.insert(0,'python')
            self.copy_command[1] = shutil.which(self.copy_command[1])

        self.tc_command = tc_command
        """Program and arguments needed to run the tool chain."""

        self.base_sim_dirn = base_dir
        """Directory containing base simulation that will be copied to MC
        directories"""

        self.keep_file_list = []
        """A list of files/globs, relative to simulation instance directory,
        to retain in each instance directory."""

        if keep_file:
            if not os.path.isfile(keep_file):
                raise RuntimeError(f'Could not find keep file {keep_file}')
            self.keep_file_list = PyPowerShellRunner.read_keep_file(keep_file)
        else:
            logger.warn('No keep file specified. Most simulation outputs '\
                    'will be discarded')

        _bdn = os.path.split(os.path.abspath(base_dir))[1]
        self._bd_tmp = tempfile.TemporaryDirectory(
            prefix=f'{self._make_mc_run_tag()}_{_bdn}_',
            dir=os.getcwd(),)
        """tempfile.TemporaryDirectory, a copy of base_dir in the CWD"""

        self.base_temp_dirn = self._bd_tmp.name
        """alias to temporary dir name"""

        shutil.copytree(self.base_sim_dirn,self.base_temp_dirn,
                dirs_exist_ok=True)

    def __del__(self):
        self._inst_status.shm.close()
        self._inst_status.shm.unlink()
        self._inst_ec.shm.close()
        self._inst_ec.shm.unlink()

        if self._bd_tmp:
            del self._bd_tmp


    def gen_mc_instances(self):
        """
        """
        m,n = (self._start,self._stop)
        width = int(ceil(log10(n)))
        for i in range(m,n):
            s=self._make_mc_run_tag('{ii:0>{w}d}'.format(ii=i,w=width))
            yield s

    @staticmethod
    def _setup_instance(d, copy_command, base_temp_dirn, isnam, iecnam):
        """
            Arguments:
                d : str
                    The value generated by `gen_mc_instance`

                copy_command : str
                    The copy command and parameters

                base_temp_dirn : str
                    The temporary directory to be copied

                isnam : str
                    name of the `multiprocessing.shared_memory.ShareableList`
                    block

                iecnam : str
                    name of the `multiprocessing.shared_memory.ShareableList`
                    block
        """

        ii = int(re.search(r'_(\d+)$',d).group(1))

        inst_status = ShareableList(name=isnam)
        inst_ec = ShareableList(name=iecnam)

        if inst_status[ii] != InstanceStatusEnum.NOT_STARTED.value:
            logger.info(f'Instance {d} already set up')
            return
        else:
            inst_status[ii]=InstanceStatusEnum.SET_UP_RUNNING.value

        logger.info(f'Setting up MC run instance {d}')

        # assume last two arguments in the copy command will be the source and
        # target directories
        #import pdb ; pdb.set_trace()
        cp = subprocess.run(copy_command+[base_temp_dirn,d,])

        if cp.returncode == 0:
            inst_status[ii] = InstanceStatusEnum.SET_UP_SUCCESS.value
        else:
            inst_status[ii] = InstanceStatusEnum.SET_UP_FAILED.value

        inst_ec[ii] = cp.returncode

        inst_status.shm.close()
        inst_ec.shm.close()

        return cp.returncode

    @staticmethod
    def _run_instance(d, copy_command, base_temp_dirn, tc_command,
            keep_file_list, isnam, iecnam, log_level,):
        """Launch an HGS Run instance.

        Returns the exit code of the copy command if copying was unsuccessful.
        Othewise, it returns the exit code of the tool chain.

        Arguments:
            d : str
                The value generated by `gen_mc_instance`

            copy_command : str
                The copy command and parameters

            base_temp_dirn : str
                The temporary directory to be copied

            tc_command : str
                The toolchain running command and arguments

            keep_file_list : list-like
                list of strings or globs to retain in the output directory

            isnam : str
                name of the `multiprocessing.shared_memory.ShareableList`
                block
                Status of this instance is stored at index `d`

            iecnam : str
                name of the `multiprocessing.shared_memory.ShareableList`
                block
                Exit code of this instance is stored at index `d`

            log_level : int
                Logging level for this instance's subprocess
        """

        # Local variables for shared memory
        inst_status = ShareableList(name=isnam)
        inst_ec = ShareableList(name=iecnam)

        # instance index
        ii = int(re.search(r'_(\d+)$',d).group(1))

        # remake logger object in this thread
        mylogger = logging.getLogger(__name__)
        _logger_setup(mylogger)
        mylogger.setLevel(log_level)

        # send feedback to user
        mylogger.info(f'Launching instance {ii+1}')
        mylogger.verbose1(f'Setting up {d}')

        HGS_MCRunner._setup_instance(d, copy_command, base_temp_dirn, isnam,
                iecnam)

        mylogger.verbose2(
            f'Set up status for {d}: ' \
            f'{InstanceStatusEnum(inst_status[ii]).name}')

        if inst_status[ii] != InstanceStatusEnum.SET_UP_SUCCESS.value:
            return inst_ec[ii]
        
        mylogger.verbose1(f'Running {d}')
        inst_status[ii] = InstanceStatusEnum.RUNNING.value

        runner = PyPowerShellRunner(
                tc_command,
                simdir=d,
                timeout=600000.)

        returncode,stdout,stderr = runner.runSim()

        if returncode == 0:
            # try to keep all input files generated by makeMesh
            runner.eraseRunOutputs(keep=keep_file_list)
            inst_status[ii] = InstanceStatusEnum.FINISHED_SUCCESS.value
        else:
            inst_status[ii] = InstanceStatusEnum.FINISHED_FAILED.value

        inst_ec[ii] = returncode

        mylogger.verbose2(
            f'Run exit status for {d}: '\
            f'{InstanceStatusEnum(inst_status[ii]).name}')

        inst_status.shm.close()
        inst_ec.shm.close()

        return returncode

    def run(self, num_processes):
        """Run series instances."""
        _N = self._stop - self._start
        # run
        results = []
        if num_processes > 1:
            with Pool(num_processes) as pool:
                results_objs = []
                for ii,inst in enumerate(mc.gen_mc_instances()):
                    logger.info(f'Enqueuing instance {ii+1} of {_N}')
                    results_objs.append(pool.apply_async(
                            HGS_MCRunner._run_instance, (inst,
                                self.copy_command, self.base_temp_dirn,
                                self.tc_command,
                                self.keep_file_list,
                                self._inst_status.shm.name,
                                self._inst_ec.shm.name,
                                logger.getEffectiveLevel(),
                            )))

                    # There is some unknown issue that occurs (in powershell 'batch
                    # file cannot be found') when many (>=4) processes are launced
                    # at the same time. Insert a delay in the first batch of
                    # processes to try to avoid the issue.
                    if ii < args.num_processes:
                        time.sleep(2)

                results = [r.get() for r in results_objs]

        else:
            for inst in mc.gen_mc_instances():
                results.append(
                        HGS_MCRunner._run_instance(inst,
                            self.copy_command, self.base_temp_dirn,
                            self.tc_command,
                            self.keep_file_list,
                            self._inst_status.shm.name,
                            self._inst_ec.shm.name,
                            logger.getEffectiveLevel(),
                        ))

        return all( r==0 for r in results )

    def _make_mc_run_tag(self,suff=''):
        """Return a name for this run (no suff) or MC instance (with suff)"""

        name = f'mc{"_dry" if self.DRY_RUN else ""}'\
               f'_{self.RUN_START.strftime("%Y%m%d_%H%M")}'

        if suff:
            name += f'_{suff}'

        return name

    def __str__(self):
        return self._make_mc_run_tag()


    def getRunReport(self):
        """Return a string with the count of instances in each status category
        """

        s = f'N={self._stop-self._start}; ' \
            f'start..stop={self._start}..{self._stop-1}\n'

        count = dict( (k.value,0) for k in InstanceStatusEnum )

        for ist in range(self._start,self._stop):
            st = self._inst_status[ist]
            count[st] += 1

        for k in InstanceStatusEnum:
            if count[k.value]:
                s += f'N_{k}={count[k.value]}\n'

        return s[:-1]


if __name__ == '__main__':

    ap = argparse.ArgumentParser()

    ap.add_argument( 'n_inst',
        metavar='N_INSTANCES',
        type=int,
        help='The number of MC instances to produce',
        )

    ap.add_argument( 'base_sim', type=str,
        help='Directory containing the base simulation inputs.'
        )

    ap.add_argument( '-v', '--verbose', default=1, action='count',
        help='Set verbosity. Default 1. Set to 0 for quiet, or higher for '\
            'increasing levels of chatter.'
        )

    ap.add_argument( '--num-processes', default=1, type=int,
        metavar='N',
        help='The maximum number of concurrent processes to use.'
        )

    ap.add_argument( '--keep-file',
        default=None,
        type=str,
        help='The name of a file (relative to base_sim) that contains a '\
            'list of file names and/or glob-style patterns of files to '\
            'retain in each MC instance directory.',
        )

    ap.add_argument( '--copy-command',
        default='hgs-copy-inputs.py',
        type=str,
        help='The name of the script or a shell command that clones the ' \
            'base simulation inputs in a MC directory. ' \
            'Note: in Windows environments, the PATHEXT environment ' \
            'variable may need to be modified for "uncommon" executable ' \
            'file extensions',
        )

    ap.add_argument( '--run-command',
        default='hgs-runall.ps1',
        type=str,
        help='The name of the script or a shell command that that invokes ' \
            'preprocessing, grok, pghs, hsplot, and postprocessing/cleanup' \
            'Note: in Windows environments, the PATHEXT environment ' \
            'variable may need to be modified for "uncommon" executable ' \
            'file extensions',
        )

    ap.add_argument( '-st', '--start-at',
        dest='start_at',
        metavar='I_INSTANCE',
        default=0,
        type=int,
        help='Start at the given instance index number. Default 0.')

    args = ap.parse_args()

    if args.verbose < 1:
        logger.setLevel(logging.WARN)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)

    # check for base dir
    if not os.path.isdir(args.base_sim):
        print(f'Could not find base simulation directory "{args.base_sim}".',
            file=sys.stderr)
        sys.exit(-1)
    # autodetect the batch.pfx file
    if not 'batch.pfx' in os.listdir(args.base_sim):
        print(f'Could not find "{args.base_sim}{os.path.sep}batch.pfx".',
            file=sys.stderr)
        sys.exit(-1)

    # check for keep file
    if args.keep_file and not os.path.isfile(args.keep_file):
        ap.error(f'Could not find keep file {args.keep_file}')

    # Check if other executable files exist
    # May need posix=False here, if retaining quotes in individual
    # parameters is important
    #import pdb ; pdb.set_trace()
    for f in [ args.copy_command, args.run_command ]:
        exe = shlex.split(f)[0]
        if os.access(exe,os.X_OK):
            pass
        elif shutil.which(exe):
            pass
        else:
            print(f'Could not find {exe}.', file=sys.stderr)
            sys.exit(-1)

    mc = HGS_MCRunner(args.copy_command, args.run_command, args.base_sim,
            args.keep_file, args.n_inst, args.start_at)

    mc.run(args.num_processes)

    print(f'Run Report:\n{mc.getRunReport()}')

    del mc
    sys.exit(0)
