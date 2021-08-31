"""Classes to support executing portions of the Hydrogeosphere toolchain
"""

import sys
import os
import glob
import shutil
#import subprocess
import asyncio
import datetime

from pyhgs._py_powershell_runner import PyPowerShellRunner

__all__ = [
    'PyPowerShellRunner',
    'HGSToolChainRun',
]

class _Tee(object):
    """https://stackoverflow.com/questions/616645/how-to-duplicate-sys-stdout-to-a-log-file
        http://mail.python.org/pipermail/python-list/2007-May/438106.html
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        self.encoding = sys.stdout.encoding
        sys.stdout = self
    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()
        self.stdout.flush()

EXTENSION_RUNNER_HINTS = {
    '.ps1': [shutil.which('powershell'),'-NonInteractive','-File',],
    '.py': [shutil.which('python'),],
}
"""Dictionary of file extensions with some runners and extra arguments needed by
subprocess.run to invoke these correctly.
"""


def _list_pprocessing(
            fpfx='pre',
            dpfx='.',
            use_extension_hints=True,
        ):
    """Return a list of executables matching the name convention.

        Executables are sorted alphanumerically by file name.

        The executable is actually a list of of command line parameters that are
        useful in calls to os.subprocess
    """
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

class HGSToolChainRun():
    """Run the Hydrogeosphere tool chain"""

    def __init__(self, sim_dir):

        owd = os.getcwd()
        os.chdir(sim_dir)

        self.sim_dir = sim_dir
        self.invoke_time = datetime.datetime.now()

        try:
            self.prefix = open('batch.pfx','r').read().strip()
        except:
            os.chdir(owd)
            raise

        # put pre- and postprocessing commands into TOOL_CHAIN
        self._tool_chain = _list_pprocessing('pre') \
            + list(
                ( [shutil.which(exe),]
                  for exe in ('grok','phgs','hsplot')) ) \
            + _list_pprocessing('post')


        os.chdir(owd)

        self._t_start = 0
        self._t_end = len(self._tool_chain)

    def check_tools(self):
        """Check file and directory read/write/execute status.

        Returns (int, [str,...], ) with ( all_tools_ok, and [ messages ] )
            with all_tools_ok = 0 if everything is fine, or
                all_tools_ok = 1 some problem.

        """
        owd = os.getcwd()
        os.chdir(self.sim_dir)

        msgs = []
        allOk = True

        if not os.access('.', os.R_OK):
            msgs.append(f'False. Cannot read files in {self.sim_dir}.')
            allOk &= False

        if not os.access('.', os.W_OK):
            msgs.append(f'False. Cannot write files in {self.sim_dir}.')
            allOk &= False

        for n,p in enumerate(self._tool_chain):
            tmsg = 'SKIPPED'

            if self._t_start <= n < self._t_end:
                ok = os.access(p[0], os.X_OK)
                allOk &= ok
                tmsg = f'{"OK" if ok else "NOT OK"}'

            itemstr = f'{n})'
            clstr = ' '.join(p)
            msgs.append(
                f'{itemstr:<4} {tmsg:7} {clstr}' )

        os.chdir(owd)

        if allOk:
            return (0, msgs)
        return (1, msgs)

    def set_start(self,val):
        """Sets the starting tool in the chain

        Arguments:
            val : str or int
                if int, interpret this as the index
                if str, search for a tool that matches
        """
        if not val:
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


    def set_end(self,val):
        """Sets the ending tool in the chain

        Arguments:
            val : str or int
                if int, interpret this as the index
                if str, search for a tool that matches
        """
        if not val:
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
            if iend >= len(self._tool_chain)-1:
                raise RuntimeError('Cannot have end index > number of tools')

        if iend < 0:
            raise RuntimeError(f'Tool {val} not found in the chain.')

        if self._t_start > iend:
            raise RuntimeError('Cannot have ending tool < starting tool')

        # adopt end+1 convention
        self._t_end = iend+1

    def run(self):
        """Run the toolchain and return the exit code"""

        return asyncio.run(self._run())

    async def _run(self):
        """Run the toolchain

        All console output is passed through.
        """

        consolefn = 'console.txt'
        returncode = 0
        errmsgs = ''

        with open(consolefn,'a') as fout:
            msg = '\n'+80*'='+f'\n{self!s}\n'
            print(msg,file=fout)
            print(msg)

        for i,tool_args in enumerate(
                self._tool_chain[self._t_start:self._t_end],
                start=self._t_start):

            outputtee = _Tee(consolefn,'a')

            # tool starting message
            outputtee.write(
                  80*'='+f'\n{i}) {" ".join(tool_args)}\n'
                  + f'at {datetime.datetime.now()}\n\n'
                  )

            try:
                proc = await asyncio.subprocess.create_subprocess_exec(
                    *tool_args,
                    #stdin=subprocess.PIPE,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.STDOUT,
                    limit=256,
                    #check=True,
                )

            except Exception as e:
                errmsgs = str(e)
                returncode = 1 # signal fail

            else:
                # poll the output streams and return code
                while True:
                    (pstdout, pstderr) = await proc.communicate()
                    # pstderr shoulde be '' from PIPE above
                    outputtee.write(pstdout.decode('ascii'))

                    if proc.returncode is not None:
                        break

                # handle any final output
                (pstdout, pstderr) = await proc.communicate()
                outputtee.write(pstdout.decode('ascii'))

                await proc.wait()

                returncode = proc.returncode
            
            del outputtee

            if returncode != 0:
                break

            with open(consolefn,'r') as console:
                pass
                # look for other ways that this may have failed. Sometimes
                # HGS has exit status of 0 (successs) but a failure message
                # in the console's text
                
        if returncode != 0:
            msg =  f'FAILED: HGS tool chain for {self.getIdStr()}'
            if errmsgs:
                msg += ' '+errmsgs
            print(msg, file=sys.stderr)
                
        return returncode

    def getIdStr(self):
        return f'{self.prefix} in {self.sim_dir} @ {self.invoke_time!s}'

    def __str__(self):
        toolstrs = []
        for i,t in enumerate(self._tool_chain):
            skipstr = 'SKIPPED '
            if self._t_start <= i < self._t_end:
                skipstr = ''
            toolstrs.append( f'{i}) {skipstr}'+' '.join(t) )

        return 'Hydrogeosphere tool chain will run simulation ' \
            + f'{self.getIdStr()}\n' \
            + '\n'.join(toolstrs)


