"""Classes to support executing portions of the Hydrogeosphere toolchain
"""

import os
import glob
import shutil

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




def _list_pprocessing(fpfx='pre',dpfx='.',add_python=True):
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

    if add_python:
        pyexe = shutil.which('python')
        for n in names:
            if n.lower().endswith('.py'):
                ret.append( [ pyexe, n, ] )
            else:
                ret.append( [ n, ] )

        return ret

    ret = [ [n,] for n in names ]

    return ret

class HGSToolChainRun():
    """Run the Hydrogeosphere tool chain"""

    def __init__(self, sim_dir):

        owd = os.getcwd()
        os.chdir(sim_dir)

        self.sim_dir = sim_dir
        with open('batch.pfx','r') as fin:
            self.prefix = fin.read().strip().split()[0:]

        # put pre- and postprocessing commands into TOOL_CHAIN
        self._tool_chain = _list_pprocessing('pre') \
            + list(
                ( [shutil.which(exe),]
                  for exe in ('grok','phgs','hsplot')) ) \
            + _list_pprocessing('post')


        os.chdir(owd)

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


    def __str__(self):
        toolstrs = []
        for i,t in enumerate(self._tool_chain):
            toolstrs.append( f'{i}) '+' '.join(t) )

        return 'Hydrogeosphere tool chain will run simulation ' +\
            f'{self.prefix} in directory {self.sim_dir}:\n' +\
            '\n'.join(toolstrs)

