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




def _find_pprocessing(dpfx='.',fpfx='pre'):
    """Return a sorted list of executable files matching the name convention"""
    if fpfx not in ['pre','post']:
        raise ValueError('fpfx must be "pre" or "post"')

    return list( os.path.basename(f)
            for f in glob.glob(os.path.join(dpfx,f'{fpfx}process*'))
            if os.access(f,os.X_OK)
        ).sort()

class HGSToolChainRun():
    """Run the Hydrogeosphere tool chain"""

    def __init__(self, sim_dir):
    
        owd = os.getcwd()
        os.chdir(sim_dir)

        self.sim_dir = sim_dir
        with open('batch.pfx','r') as fin:
            self.prefix = fin.read().strip().split()[0:]

        # put pre- and postprocessing commands into TOOL_CHAIN
        self.tool_chain = _find_pprocessing('','pre') \
            + [ shutil.which('grok'),
                shutil.which('phgs'),
                shutil.which('hsplot'), ] \
            + _find_pprocessing('','post')


        os.chdir(owd)

    def check_tools(self):
        """Check file and directory read/write/execute status.

        Returns (int, [str,...], ) with ( all_tools_ok, and [ messages ] )
            with all_tools_ok = 0 if everything is fine, or
                all_tools_ok = 1 some problem.

        """
        owd = os.getcwd()
        os.chdir(sim_dir)

        msgs = []
        allOk = True

        if not os.access('.', os.R_OK):
            msgs.append(f'False. Cannot read files in {self.sim_dir}.')
            allOk &= False

        if not os.access('.', os.W_OK):
            msgs.append(f'False. Cannot write files in {self.sim_dir}.')
            allOk &= False

        for n,p in enumerate(self.tool_chain):
            ok = os.access(p, os.X_OK)
            allOk &= ok
            itemstr = f'{n})'
            msgs.append( f'{itemstr:<4} {ok:5} {p}'

        os.chdir(owd)

        if allOk:
            return (0, msgs)
        return (1, msgs)


    def __str__(self):
        return 'Hydrogeosphere tool chain will run simulation '\
            f'{self.prefix} in directory {self.sim_dir}:'
            '\n'.join( f'{n}) {p}' for n,p in enumerate(self.tool_chain) )

