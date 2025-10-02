"""Tools for running and managing the Hydrogeosphere (HGS) toolchain.

This module provides the core `HGSToolChainRun` class responsible for:
1. Discovering and sequencing external HGS executables (grok, phgs, hgs2vtu, preprocess/postprocess scripts).
2. Asynchronously executing the toolchain using Python's `asyncio`.
3. Providing real-time, filtered console output for monitoring simulation progress.
4. Buffering all raw output (using a bounded `deque`) for immediate dumping upon error or user interruption.

The module uses modular filter classes (e.g., `PhgsFilter`) to suppress verbose
solver output and display crucial status updates and performance metrics (like
simulation speed and timestep failures).
"""
import sys
import os
import glob
import shutil
import asyncio
import datetime
import re
import logging
import shlex
from collections import deque
from itertools import chain
from pathlib import Path

from hgstools.pyhgs.cli import parse_path_to_prefix

logger = logging.getLogger(__name__)

# --- Filter Classes for Console Output ---
# These filter specific tool output to suppress noise and highlight status.

class BaseOutputFilter:
    """Base class for output filtering.

    The base class implements common rules for all tools, ensuring important
    messages (like errors and tool banners) are always printed to the console.

    The `should_print` method returns a string if a custom summary is generated
    (e.g., for PHGS progress), False if the line should be suppressed, or
    True if the raw line should be printed.
    """

    def __init__(self, filtering_on=True):
        """
        Parameters
        ----------
        filtering_on : bool, optional
            If False, all output is passed through without filtering,
            by default True.
        """
        self._filtering_on = filtering_on
        self._EXIT_FILTER = re.compile(r'-+ \w* exit -+', re.I)

    def set_filtering_status(self, status):
        """Sets the filter status.

        Parameters
        ----------
        status : bool
            True to enable filtering, False to disable (pass all output).
        """
        self._filtering_on = status

    def _tool_specific_rules(self, line_str):
        """Implement tool-specific filtering logic here.

        Parameters
        ----------
        line_str : str
            The raw line of text from the subprocess stdout.

        Returns
        -------
        bool or str
            - bool: True to print line, False to suppress.
            - str: A formatted summary line to print instead of the raw line.
        """
        # Suppress everything by default unless explicitly allowed by rules
        return False

    def should_print(self, line_str):
        """Determines if a line should be printed to the console.

        Parameters
        ----------
        line_str : str
            The raw line of text from the subprocess stdout.

        Returns
        -------
        bool or str
            - If filtering is OFF, returns True (to print the raw line).
            - If filtering is ON:
                - bool: True to print line, False to suppress.
                - str: A formatted summary line to print instead of the raw line.
        """
        # If filtering is off, immediately return True to pass all output through
        if not self._filtering_on:
            return True

        # Strip whitespace for consistent comparison
        line = line_str.strip()

        # 1. Keep separator banners as a break between tools
        # Updated to check for '@@' to catch all parts of the HGS tool banners.
        if line.startswith('@@'):
            return True

        # 2. Keep final exit/summary lines
        # uhoh - this also catches 'Step' lines. Omitting for now.
        #if line.startswith('---') or line.startswith('****'):
        if self._EXIT_FILTER.match(line):
            return True

        # 3. Keep critical status messages
        if any(keyword in line.upper() for keyword in ['WARNING', 'ERROR', 'FAILED', 'FATAL', 'RAISED']):
            return True

        # 4. Apply tool-specific rules
        return self._tool_specific_rules(line_str)


class GrokFilter(BaseOutputFilter):
    """Filtering rules specific to the GROK pre-processor.

    This filter suppresses all general output from GROK, allowing only banners,
    warnings, errors, and exit status to pass through.
    """
    def _tool_specific_rules(self, line_str):
        line = line_str.strip()
        # Allow lines containing specific keywords (in addition to base filter)
        if any(keyword in line.upper() for keyword in ['WARNING', 'ERROR', 'RAISED']):
            return True
        # Suppress everything else by default
        return False


class PhgsFilter(BaseOutputFilter):
    """Filtering rules specific to the PHGS solver.

    This filter tracks transient simulation progress (time, dt, failures) and outputs a
    single summary line per accepted timestep for performance monitoring, while
    suppressing verbose iteration details.
    """
    # Note: Regex patterns are constructed dynamically in __init__
    _PHGS_ACCEPTED_TIME_RE = re.compile(r'Accepted\s+solution\s+at\s+time:\s+([\d.E+-]+)')

    def __init__(self, prefix, filtering_on=True):
        """
        Parameters
        ----------
        prefix : str
            The simulation's file prefix (e.g., 'module1b'), used to construct
            the dynamic regular expression for step capture.
        filtering_on : bool, optional
            If False, all output is passed through without filtering,
            by default True.
        """
        super().__init__(filtering_on=filtering_on)
        self._prefix = prefix
        self._step_data = {}  # Tracks data for the current step
        self._last_wall_time = None
        self._in_time_report = False # State flag for Time Report processing

        # Dynamic Regex construction using the simulation prefix
        # Matches lines like: "----------------------module1b Step: 4566------------------------"
        self._PHGS_STEP_RE = re.compile(
                rf'^-+\s*{re.escape(self._prefix)}\s+Step:\s*(\d+)\s*-+\s*$', flags=re.I)

        # Matches lines like: " 1.9000598038E-05 1.0000843927E+01 Accept timestep"
        self._PHGS_DT_RE = re.compile(r'([\d.E+-]+)\s+[\d.E+-]+\s+Accept\s+timestep')

    def _format_summary(self):
        """Generates the single summary string for an accepted timestep."""
        step = self._step_data.get('step', 'N/A')
        sim_time = self._step_data.get('time', 0.0)
        dt_sim = self._step_data.get('dt', 0.0)
        failures = self._step_data.get('failures', 0)
        current_wall_time = datetime.datetime.now()

        # Calculate wall clock interval and speed
        wall_interval_sec = 0.0
        speed = 0.0
        if self._last_wall_time is not None:
            wall_interval = current_wall_time - self._last_wall_time
            wall_interval_sec = wall_interval.total_seconds()
            if wall_interval_sec > 0:
                speed = dt_sim / wall_interval_sec

        self._last_wall_time = current_wall_time
        self._step_data = {} # Clear step data after summary is generated

        # Format output
        return (
            f"PHGS Step {step:>5} | "
            f"Time: {sim_time:.4e} | "
            f"dt: {dt_sim:.4e} | "
            f"Wall Interval: {wall_interval_sec:.3f} s | "
            f"Speed: {speed:.2e} | "
            f"Failed Attempts: {failures}"
        )

    def _tool_specific_rules(self, line_str):
        line = line_str.strip()

        # 1. State check for SIMULATION TIME REPORT block
        if self._in_time_report:
            # If we hit the final dashed line, end the report state.
            if line.startswith('----------------'):
                self._in_time_report = False
            return True # Print all content inside the report block

        # 2. State trigger for SIMULATION TIME REPORT block
        if 'SIMULATION TIME REPORT' in line:
            self._in_time_report = True
            return True # Print the header

        # --- Transient Simulation Logic (only runs if NOT in time report) ---

        # 3. New Timestep (Start of state tracking)
        # Matches the dynamic separator line, e.g., "----------------------module1b Step: 4566------------------------"
        match_step = self._PHGS_STEP_RE.search(line)
        if match_step:
            # Initialize/reset state for the new step
            self._step_data = {'step': int(match_step.group(1)), 'failures': 0}
            return False # Suppress the raw step line

        # 4. Timestep Redo (Track failures)
        if 'redo timestep' in line:
            self._step_data['failures'] = self._step_data.get('failures', 0) + 1
            return False # Suppress the redo message

        # 5. Timestep Size (dt) Capture
        # Matches the line containing the simulation delta_t and "Accept timestep"
        match_dt = self._PHGS_DT_RE.search(line)
        if match_dt and 'Accept timestep' in line:
            self._step_data['dt'] = float(match_dt.group(1))
            return False # Suppress the raw dt line

        # 6. Accepted Time (End of state tracking and summary generation)
        match_time = self._PHGS_ACCEPTED_TIME_RE.search(line)
        if match_time:
            self._step_data['time'] = float(match_time.group(1))
            return self._format_summary() # Print custom summary string and suppress raw line

        # 7. Keep calculation lines (including steady-state)
        if 'flow solution' in line:
            return True

        # 8. Suppress everything else
        return False


class Hgs2VtuFilter(BaseOutputFilter):
    """Filtering rules specific to the HGS2VTU post-processor.

    This filter suppresses the banner and exit status is handled by the base class.
    """
    def _tool_specific_rules(self, line_str):
        # If a line is not empty and not suppressed by base rules, print it.
        # This keeps the general parsing/writing info.
        return bool(line_str.strip())

# --- Tool Runner Classes ---

class BaseRunner():
    """Base class for HGS-related runners (e.g., `HGSToolChainRun`).

    This class provides common utilities for managing a simulation environment,
    such as capturing pre-existing files, cleaning up simulation outputs, and
    parsing keep-file rules.
    """
    def __init__(self, exe, args=[], simdir='.', timeout=60.):
        """
        Parameters
        ----------
        exe : str
            Path to the executable to run.
        args : list of str, optional
            List of command line arguments for the executable, by default [].
        simdir : str, optional
            Directory of the simulation, by default '.'.
        timeout : float, optional
            Timeout in seconds for the subprocess execution, by default 60.
        """
        self.exe = exe
        self.args = args
        self.simdir, self.prefix, junk = parse_path_to_prefix(simdir)
        self.timeout = timeout
        self.captureRunPreexisting()

        logger.debug(f'BaseRunner with prefix {self.prefix}')

    def captureRunPreexisting(self):
        """Captures a snapshot of files existing in the simulation directory.

        Files captured are used later by `eraseRunOutputs` to distinguish
        pre-existing files from generated output files.
        """
        self.preexistingfiles = set(os.listdir(self.simdir))

    @staticmethod
    def read_keep_file(keep_file=None):
        """Reads rules from a text file to identify files/categories to keep.

        A KEEP_FILE contains lists of filenames or glob-style expressions
        (for pattern matching) of files to keep/retain when cleaning up
        (deleting files from) a simulation directory. Lines are e.g.,
        "somefile.txt", "somefile.*", or "*.[mf]props".

        Lines prefixed with 'Cat=' or 'Category=' match against a file's
        category instead of its name.

        Returns
        -------
        tuple of (list of str, list of str)
            - cats: List of category match patterns.
            - exprs: List of filename match patterns.

        Parameters
        ----------
        keep_file : str, optional
            Name of a file listing filenames or glob-style expressions (on
            separate lines) of additional files to keep. By default None.

        Raises
        ------
        ValueError
            If file reading or parsing fails due to unexpected content.
        """
        if keep_file is None or not Path(keep_file).is_file():
            return [], []

        cats = []
        exprs = []
        with open(keep_file,'r') as fin:
            for l in fin.readlines():

                l = l.strip()

                if not l or l.startswith('#'):
                    # ignore blank lines or comment-only lines
                    continue

                addto = exprs

                # split off a category prefix
                m=re.match('^cat(?:egory)?=', l, re.I)
                if m:
                    addto = cats
                    l = l[m.span()[1]:]

                # split off comments
                try:
                    l = shlex.split(l,comments=True)
                except ValueError as e:
                    logger.warning(f"Skipping line due to shlex error: {l}. Error: {e}")
                    continue

                if not l:
                    continue

                # grab just the first argument of the line --- users
                # beware that if you need a space in the expression
                # then it will need to be quoted.
                l = l[0]

                addto.append(l)

        return cats, exprs

    def eraseRunOutputs(self, keep=[], keep_file=None):
        """Erase simulation output files.

        This erases any files *in the simulation directory* or subdirectory that
        did not exist prior to the creation of this object or at the last call
        to `captureRunPreexisting()`.

        Parameters
        ----------
        keep : list-like, optional
            A list of file names to keep; these will not be deleted, even if
            they did not exist when this object was instantiated, by default [].
        keep_file : str, optional
            Path to a file with additional keep rules. See `read_keep_file`.
        """

        keep_set = set()

        def deglob(g):
            # check if the filename actually looks like a glob
            if re.search(r'[*[?]',g):
                return list( fn[len(self.simdir)+len(os.sep):]
                    for fn in glob.iglob(os.path.join(self.simdir,g)))
            else:
                return [g,]

        keep_set.update( chain.from_iterable(deglob(x) for x in keep) )

        if keep_file and os.path.isfile(keep_file):
            # Note: Assuming 'Runner' is an alias or meant to be 'BaseRunner'
            # Assuming 'Runner' is an imported name referring to BaseRunner or its child
            # For correctness inside the BaseRunner module, we use BaseRunner.read_keep_file
            try:
                cat_patterns, file_patterns = BaseRunner.read_keep_file(keep_file)
                keep_set.update( chain.from_iterable(deglob(x) for x in file_patterns) )
                # Category patterns are not currently processed for deletion and need file categorization logic.
                # Ignoring cat_patterns for now as they require `categorize_files`.
            except Exception as e:
                logger.warning(f"Error reading keep file {keep_file}: {e}")

        keep_set.update(self.preexistingfiles)

        # scan for items in the keep list that have a directory component.
        keep_dirs = list(os.path.normpath(d).split(os.sep)[0]
                for d in keep if os.path.dirname(d))
        keep_set.update( keep_dirs )

        fullsimdir = os.path.abspath(self.simdir)

        for f in sorted(set(os.listdir(self.simdir)) - keep_set):
            fullf = os.path.abspath(os.path.join(self.simdir,f))
            if os.path.commonprefix([fullf,fullsimdir,]) != fullsimdir:
                # file/directory not rooted below fsimdir
                pass
            elif os.path.isdir(fullf):
                shutil.rmtree(fullf)
            elif os.path.isfile(fullf):
                os.remove(fullf)

    def _cat_rec(self, d, files):
        """Recursive categorization function.

        Returns
        -------
        tuple of (dict, set)
            - dict: Dictionary of matching file entries by category.
            - set: Set of unmatched file entries.
        """

        ret = dict()

        cat_tbd = None
        if isinstance(d, dict):
            cat_tbd = d.items()
        elif isinstance(d, list):
            cat_tbd = iter(d) # assume that this is a list of 2-tuples
        else:
            raise ValueError(f'Unexpected type {type(d)}')

        for category, tbd in cat_tbd:

            if type(tbd) in (dict, list):
                subdict = tbd
                ret[category], files = self._cat_rec(subdict, files)

            else:
                match_func = tbd
                _matches = set(filter(match_func, files))
                files -= _matches
                if category in ret:
                    ret[category] |= _matches
                else:
                    ret[category] = _matches

        return ret, files


    def categorize_files(self, files=None):
        """Categorize files in the simulation directory.

        Returns
        -------
        tuple of (dict, set)
            - dict: Dictionary of categorized files.
            - set: Set of uncategorized files.
        """

        grokfn=self.prefix+'.grok'
        _input_pats = [
            'batch.pfx',
            grokfn,
            self.prefix+'*.control',
            'debug.control',
            'parallelindx.dat',
            'array_sizes.default',
        ]

        try:
            # Assuming grok_parse returns a dict-like object where keys starting
            # with 'files_' contain lists of filenames/patterns.
            grok = grok_parse(grokfn)
            _input_pats += list(
                chain.from_iterable( grok[a] for a in [
                    k for k in grok.keys() if k.startswith('files_') ])
            )
        except BaseException as e:
            logger.warning(f'Problem reading {grokfn}. Error:\n{e}')

        def _match(pats, p):
            for pat in pats:
                if re.match(pat, str(p)):
                    return True
            return False

        def _match_scratch(p):
            _pats = ['scratch_grok',
                'scratch_[cmfod]props',
                'scratch_etprops',
                'scratch_wellprops',
                'scratch_tileprops',
            ]
            return _match(_pats, p)


        def _match_input(p):
            '''Return True if Path p matches any known input file name pattern'''
            return _match(_input_pats, p)

        def _match_dbg(p):
            _pats = ['hs.dbg', 'grok.dbg','hsplot.dbg', ]
            return _match(_pats, p)


        if files is None:
            files = set(Path(self.simdir).iterdir())

        # convert Path objects to strings for categorization logic
        file_strings = {f.name for f in files}

        # final datastructure
        categories = {
            'input':(lambda p: _match_input(p)),
            'subdir':(lambda p: Path(self.simdir).joinpath(p).is_dir()),
            'output':{
                'scratch':_match_scratch,
                'dbg':_match_dbg,
                'intermediate':(lambda p:
                    re.match(self.prefix+r'o\..*\.\d+$', str(p))),
                'tecplot':(lambda p:
                    re.match(self.prefix+r'o\..*dat$', str(p))),
                'other':(lambda p:
                    re.match(self.prefix+r'o\..*', str(p))),
            },
        }

        cat, uncategorized_files = self._cat_rec(categories, file_strings)

        # logger.debug(pformat(cat, indent=2))

        return cat, uncategorized_files

    def write_excerpts(self, touch=False):
        """Write excerpts of o.eco and o.lst files.

        Parameters
        ----------
        touch : bool, optional
            If True, the files will be touched in the current directory. If False,
            the default, the excerpt files will be generated, by default False.

        Returns
        -------
        list of Path
            A list of the file paths of the generated excerpts.
        """
        ret = []

        totry = [
            Path(self.simdir) / (self.prefix+suff)
                for suff in 'o.eco o.lst'.split()
        ]

        for f in totry:
            pout = f.with_stem(self.prefix+'o.excerpt')

            if f.exists():
                ret.append(pout)
                if touch:
                    pout.touch()
                else:
                    # gather the original file creation and modification times
                    # and size
                    # write a header with this info

                    _fstat = f.stat()
                    fattr = {'size':None, 'ctime':None, 'mtime':None}
                    fattr['size'] = _fstat.st_size
                    fattr['mtime'] = datetime.datetime.fromtimestamp(_fstat.st_mtime)
                    try:
                        fattr['ctime'] = \
                            datetime.datetime.fromtimestamp(_fstat.st_birthtime)
                    except AttributeError:
                        fattr['ctime'] = \
                            str(datetime.datetime.fromtimestamp(_fstat.st_ctime)) \
                            + ' (possibly time of last metadata change, not original creation)'


                    s = f'Excerpt written by {__name__} on {datetime.datetime.now()}\n' \
                      + f'Original file: {f}, {fattr["size"]} bytes\n' \
                      + f'Creation time: {fattr["ctime"]}\n' \
                      + f'Modification time: {fattr["mtime"]}\n\n'

                    # do the excerpt
                    with open(pout, 'w', encoding='utf-8') as outfile:
                        outfile.write(s)
                        excerpt_large_file(str(f), outfile, 20, 21)

        return ret


class PyPowerShellRunner(object):
    """A placeholder class for running PowerShell scripts on Windows."""
    def __init__(self, script_path):
        """
        Parameters
        ----------
        script_path : str
            The path to the PowerShell script.
        """
        self.script_path = script_path
        # Implementation details for running PowerShell scripts go here...

def _list_pprocessing(
            fpfx='pre',
            dpfx='.',
            use_extension_hints=True,
        ):
    """Return a list of executables matching the name convention.

    Executables are sorted alphanumerically by file name.

    Parameters
    ----------
    fpfx : str, optional
        File prefix, either 'pre' or 'post', by default 'pre'.
    dpfx : str, optional
        Directory prefix to search in, by default '.'.
    use_extension_hints : bool, optional
        Whether to use file extension hints to determine the runner (e.g.,
        Python or PowerShell), by default True.

    Returns
    -------
    list of list of str
        A list of command line parameters suitable for `subprocess.run`.

    Raises
    ------
    ValueError
        If `fpfx` is not 'pre' or 'post'.
    """
    EXTENSION_RUNNER_HINTS = {
        '.ps1': [shutil.which('powershell'),'-NonInteractive','-File',],
        '.py': [shutil.which('python'),],
    }

    if fpfx not in ('pre','post'):
        raise ValueError('fpfx must be "pre" or "post"')

    names = sorted( ( os.path.basename(f)
              for f in glob.glob(os.path.join(dpfx,f'{fpfx}process*'))
              if os.access(f,os.X_OK)
        ) )

    ret = []

    if use_extension_hints:
        for n in names:
            ext = os.path.splitext(n)[1].lower()
            if ext in EXTENSION_RUNNER_HINTS:
                ret.append( EXTENSION_RUNNER_HINTS[ext] + [n,] )
            else:
                ret.append( [ n, ] )
    else:
        ret = [ [n,] for n in names ]

    return ret

class HGSToolChainRun(BaseRunner):
    """Run the Hydrogeosphere tool chain"""

    def __init__(self, simdir):
        """
        Parameters
        ----------
        simdir : str
            Directory containing simulation inputs path/to/sim/.

        Raises
        ------
        FileNotFoundError
            If 'batch.pfx' is not found in `simdir`.
        """
        # BaseRunner.__init__ handles captureRunPreexisting
        super().__init__(exe=None, simdir=simdir)

        owd = os.getcwd()
        os.chdir(self.simdir)

        self.invoke_time = datetime.datetime.now()

        # The prefix is already set by BaseRunner.__init__

        # put pre- and postprocessing commands into TOOL_CHAIN
        # Default core tools are grok, phgs, and hgs2vtu
        self._tool_chain = _list_pprocessing('pre') \
            + list(
                ( [shutil.which(exe),]
                  for exe in ('grok','phgs','hgs2vtu')) ) \
            + _list_pprocessing('post')


        os.chdir(owd)

        self._t_start = 0
        self._t_end = len(self._tool_chain)

    def set_tool_override(self, default_tool, optional_tool):
        """Overrides a default tool in the chain with an optional one.

        This method is designed to swap a tool like 'hgs2vtu' for 'hsplot'
        if the command-line flag is set.

        Parameters
        ----------
        default_tool : str
            The name of the tool to be replaced (e.g., 'hgs2vtu').
        optional_tool : str
            The name of the replacement tool (e.g., 'hsplot').

        Raises
        ------
        RuntimeError
            If the `optional_tool` executable cannot be found on the system path.
        """
        # Find the path to the optional tool first
        tool_path = shutil.which(optional_tool)
        if not tool_path:
            raise RuntimeError(f"Optional tool '{optional_tool}' not found in system path.")

        # Replace the tool in the chain
        for i, tool_list in enumerate(self._tool_chain):
            # The executable name is always the first item in the list
            base_exe = os.path.basename(tool_list[0])
            if base_exe == default_tool:
                # Replace the entire tool entry (path and arguments)
                self._tool_chain[i] = [tool_path]
                return # Assumes only one instance needs to be replaced

    def _check_write_access(self):
        """Tests that the simulation directory and key output files are writable.

        Creates and deletes a temporary file to check directory write access.
        Also checks write permission for existing key output files.

        Raises
        ------
        RuntimeError
            If write access fails for the directory or for key existing files.
        """
        check_files = [self.prefix+'o.eco', self.prefix+'o.lst']
        test_file = os.path.join(self.simdir, '__hgs_write_test__')

        # 1. Check directory writability by creating a test file
        try:
            with open(test_file, 'w') as f:
                f.write('test')
            os.remove(test_file)
        except Exception as e:
            raise RuntimeError(f"Directory {self.simdir} is not writable: {e}")

        # 2. Check writability of existing output files
        for f in check_files:
            file_path = os.path.join(self.simdir, f)
            if os.path.exists(file_path) and not os.access(file_path, os.W_OK):
                raise RuntimeError(f"Existing output file '{f}' is not writable. Check permissions.")


    def check_tools(self):
        """Check file and directory read/write/execute status.

        Includes a check for write access to the simulation directory and key output files.

        Returns
        -------
        tuple of (int, list of str)
            - returncode: 0 if everything is fine, 1 if there is a problem.
            - msgs: A list of messages detailing the status of each check.
        """
        owd = os.getcwd()
        os.chdir(self.simdir)

        msgs = []
        allOk = True
        returncode = 0

        # --- 1. Check Read/Write permissions (Comprehensive Check) ---
        try:
            self._check_write_access()
            msgs.append(f'OK. Directory {self.simdir} and key files are writable.')
        except RuntimeError as e:
            msgs.append(f'NOT OK. Write access failed: {e}')
            allOk = False
            returncode = 1

        if not os.access('.', os.R_OK):
            msgs.append(f'NOT OK. Cannot read files in {self.simdir}.')
            allOk = False
            returncode = 1
        else:
             msgs.append(f'OK. Can read files in {self.simdir}.')


        # --- 2. Check Tool Executability ---
        for n,p in enumerate(self._tool_chain):
            tmsg = 'SKIPPED'

            if self._t_start <= n < self._t_end:
                ok = os.access(p[0], os.X_OK)
                if not ok:
                    returncode = 1
                    allOk = False

                tmsg = f'{"OK" if ok else "NOT OK"}'

            itemstr = f'{n})'
            clstr = ' '.join(p)
            msgs.append(
                f'{itemstr:<4} {tmsg:7} {clstr}' )

        os.chdir(owd)

        return (returncode, msgs)

    def set_start(self, val):
        """Sets the starting tool in the chain.

        Parameters
        ----------
        val : str or int or None
            If int, interpret this as the tool index. If str, search for a tool
            that matches (e.g., 'grok'). If None, start from the beginning.

        Raises
        ------
        RuntimeError
            If the index is negative, if the tool string is not found, or if
            the start index is after the end index.
        """
        if val is None:
            return

        istart = 0
        try:
            istart=int(val)
        except ValueError:
            # try string search
            istart = -1
            for i,t in enumerate(self._tool_chain):
                tstr = ' '.join(t)
                if val in tstr:
                    istart = i
                    break
        else:
            if istart < 0:
                raise RuntimeError('Cannot have negative start index')

        if istart < 0:
            raise RuntimeError(f'Tool {val} not found in the chain.')

        if istart > self._t_end-1:
            raise RuntimeError('Cannot have starting tool > ending tool')

        self._t_start = istart


    def set_end(self, val):
        """Sets the ending tool in the chain.

        Parameters
        ----------
        val : str or int or None
            If int, interpret this as the tool index. If str, search for a tool
            that matches. If None, run to the end.

        Raises
        ------
        RuntimeError
            If the index is negative, if the tool string is not found, or if
            the end index is before the start index.
        """
        if val is None:
            return

        iend = 0
        try:
            iend=int(val)
        except ValueError:
            # try string search
            iend = -1
            for i,t in reversed(list(enumerate(self._tool_chain))):
                tstr = ' '.join(t)
                if val in tstr:
                    iend = i
                    break
        else:
            if iend < 0:
                raise RuntimeError('Cannot have negative end index')
            if iend >= len(self._tool_chain):
                # Allow one past the end as a valid index for convenience
                # but flag an error if it's much greater.
                if iend > len(self._tool_chain):
                    raise RuntimeError('End index out of range')


        if iend < 0:
            raise RuntimeError(f'Tool {val} not found in the chain.')

        if self._t_start > iend:
            raise RuntimeError('Cannot have ending tool < starting tool')

        # adopt end+1 convention
        self._t_end = iend+1

    def run(self):
        """Run the toolchain and return the exit code."""

        # The access check is done here, right before execution
        self._check_write_access()
        return asyncio.run(self._run())

    async def _run(self):
        """Run the toolchain core logic asynchronously.

        All console output is passed through a filtered stream, and a bounded
        buffer captures the raw output for error dumping on failure.

        Returns
        -------
        int
            The return code of the last executed tool (0 on success).
        """

        full_output_buffer = deque(maxlen=200)
        returncode = 0
        errmsgs = ''
        owd = os.getcwd()
        start_time = datetime.datetime.now()
        wall_time_format = '%Y-%m-%d %H:%M:%S'

        # Initialize the filters dynamically
        filters = {
            'grok': GrokFilter(filtering_on=True),
            'phgs': PhgsFilter(prefix=self.prefix, filtering_on=True),
            'hgs2vtu': Hgs2VtuFilter(filtering_on=True),
        }

        async def stream_output(proc_stdout, filter_instance, buffer):
            """Reads subprocess stdout line-by-line, buffers it, and applies filter."""
            try:
                # Iterate over lines read from the pipe
                async for line in proc_stdout:
                    line_str = line.decode('ascii', errors='ignore')
                    buffer.append(line_str)

                    # Apply filter logic
                    print_status = filter_instance.should_print(line_str)

                    if print_status is True:
                        sys.stdout.write(line_str)
                    elif isinstance(print_status, str):
                        # This is a custom summary line from the filter
                        sys.stdout.write(print_status + '\n')

                    sys.stdout.flush()

            except Exception as e:
                # Log streaming errors, but don't stop the whole process immediately
                logger.error(f"Error streaming output: {e}")


        def _dump_buffer(buffer, error_message):
            """Dumps the full buffered output to stderr on error/interrupt."""
            sys.stderr.write('\n\n' + 80*'#' + '\n')
            sys.stderr.write('### ERROR: Simulation Interrupted or Failed ###\n')
            sys.stderr.write(f'### Reason: {error_message} ###\n')
            sys.stderr.write('### Last 200 lines of buffered output: ###\n')
            sys.stderr.write(80*'#' + '\n')
            sys.stderr.writelines(list(buffer))
            sys.stderr.write(80*'#' + '\n\n')

        try:
            os.chdir(self.simdir)

            msg = '\n'+80*'='+f'\n{self!s}\n'
            sys.stdout.write(msg)
            sys.stdout.flush()

            for i,tool_args in enumerate(
                    self._tool_chain[self._t_start:self._t_end],
                    start=self._t_start):

                tool_name = os.path.basename(tool_args[0]).split('.')[0].lower()
                current_filter = filters.get(tool_name, BaseOutputFilter(filtering_on=False))

                # Tool starting message
                sys.stdout.write(
                    80*'='+f'\n{i}) {" ".join(tool_args)}\n'
                    + f'at {datetime.datetime.now().strftime(wall_time_format)}\n\n'
                )
                sys.stdout.flush()
                full_output_buffer.append(80*'='+f'\n{i}) {" ".join(tool_args)}\n')


                try:
                    proc = await asyncio.create_subprocess_exec(
                        *tool_args,
                        stdout=asyncio.subprocess.PIPE,
                        stderr=asyncio.subprocess.STDOUT, # Merge stderr into stdout pipe
                        limit=4096, # Increased limit for larger pipes
                        env=os.environ
                    )

                except Exception as e:
                    errmsgs = f"Failed to launch tool {tool_name}: {str(e)}"
                    returncode = 1 # signal launch failure
                    full_output_buffer.append(errmsgs + '\n')

                else:
                    # Run the streaming task concurrently with waiting for the process to finish
                    stream_task = asyncio.create_task(stream_output(
                        proc.stdout, current_filter, full_output_buffer))

                    # Wait for the process to finish
                    await proc.wait()
                    await stream_task # Ensure the streaming task finishes reading residual output

                    returncode = proc.returncode

                if returncode != 0:
                    break

        except (asyncio.CancelledError, KeyboardInterrupt) as e:
            # Catch Ctrl+C or external cancellation
            _dump_buffer(full_output_buffer, f"Process cancelled by user ({type(e).__name__})")
            return 1 # Return failure code on interrupt

        except Exception as e:
            # Catch unexpected internal Python errors
            _dump_buffer(full_output_buffer, f"Internal Runner Error ({type(e).__name__}): {e}")
            return 1

        finally:
            os.chdir(owd) # Always return to the original working directory

        if returncode != 0:
            msg = f'FAILED (code {returncode}): HGS tool chain for {self.getIdStr()}\n'
            if errmsgs:
                msg += ' ' + errmsgs
            sys.stderr.write(msg)

        return returncode

    def getIdStr(self):
        """Returns a formatted string identifying the simulation run."""
        return f'{self.prefix} in {self.simdir} @ {self.invoke_time.strftime("%Y-%m-%d %H:%M:%S")}'

    def __str__(self):
        """Returns a string representation of the toolchain configuration."""
        toolstrs = []
        for i,t in enumerate(self._tool_chain):
            skipstr = 'SKIPPED '
            if self._t_start <= i < self._t_end:
                skipstr = ''
            toolstrs.append( f'{i}) {skipstr}'+' '.join(t) )

        return 'Hydrogeosphere tool chain will run simulation ' \
            + f'{self.getIdStr()}\n' \
            + '\n'.join(toolstrs)
