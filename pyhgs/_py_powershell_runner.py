"""HGS Runner for Window Powershell. PROVISIONAL

This implementation will be removed when the pyhgsrunner is operational
"""
import sys
import os
import re
import glob
import shlex
import subprocess
import shutil
import logging
from tempfile import SpooledTemporaryFile
from itertools import chain
from .runner import BaseRunner

logger = logging.getLogger('runner')

class PyPowerShellRunner(BaseRunner):                                              #{{{
    """Run HydroGeoSphere toolchain --- special case for PowerShell
        
    """
    def __init__(self, exe, args=[], simdir='.', timeout=60.):
        self.super(exe, args, simdir, timeout)

    def runSim(self):
        """Return (exitStatus, stdout messages, stderr messages)"""

        methodInvocationDir = os.getcwd()
        os.chdir(self.simdir)

        exitStatus = 128
        out = ''
        err = ''

        spewToConsole = logger.getEffectiveLevel() >= logging.INFO
            
        outObj = None
        errObj = None

        if not spewToConsole:
            outObj = SpooledTemporaryFile(max_size=2048,mode='w+')
            errObj = SpooledTemporaryFile(max_size=2048,mode='w+')

        try:
            cp = subprocess.run(
                [ 'powershell.exe', '-noprofile', self.exe, ] + self.args,
                encoding=sys.stdout.encoding,
                #stdout=outObj,
                #stderr=errObj,
                timeout=self.timeout,
                check=True,
                env=os.environ,
                )

        except subprocess.TimeoutExpired as e:
            exitStatus = 1
            out = str(e.stdout)
            #if e.stderr:
            #    err = f'{e.stderr}\n'
            err += f"Timeout occurred. "
            err += f"Process did not complete after {self.timeout}s"

        except subprocess.CalledProcessError  as e:
            exitStatus = 1
            if not spewToConsole:
                outObj.seek(0)
                out = outObj.read()
                errObj.seek(0)
                err = errObj.read()
            else:
                out = e.stdout
                err = e.stderr
            
        except Exception as e:
            raise Exception("PyPowerShellRunner subprocess failure") from e

        else:
            exitStatus = cp.returncode
            out = cp.stdout
            err = cp.stderr

        finally:
            if not spewToConsole:
                outObj.close()
                errObj.close()
            else:
                sys.stdout.flush()
                sys.stderr.flush()

            os.chdir(methodInvocationDir)


        return (exitStatus, str(out), str(err)) 

# }}} end of PyPowerShellRunner
