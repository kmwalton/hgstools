import sys
import os
import subprocess
from tempfile import SpooledTemporaryFile

__VERBOSE__=0

class PyPowerShellRunner:                                              #{{{
    """Run HydroGeoSphere toolchain --- special case for 
        
    """
    def __init__(self, exe, args=[], simdir='.', timeout=60.):

        self.exe = exe
        self.args = args
        self.simdir = simdir
        self.timeout = timeout
        self.preexistingfiles = set(os.listdir(simdir))


    def eraseRunOutputs(self):
        """Erase any new files in the simulation directory

        'New' means any file that did not exist when this object was created.
        """
        methodInvocationDir = os.getcwd()
        os.chdir(self.simdir)
        for f in set(os.listdir()) - self.preexistingfiles:
            os.remove(f)
        os.chdir(methodInvocationDir)
    

    def runSim(self):
        """Return (exitStatus, stdout messages, stderr messages)"""

        methodInvocationDir = os.getcwd()
        os.chdir(self.simdir)

        exitStatus = 128
        out = ''
        err = ''

        spewToConsole = __VERBOSE__ > 1
            
        outObj = None
        errObj = None

        if not spewToConsole:
            outObj = SpooledTemporaryFile(max_size=2048,mode='w+')
            errObj = SpooledTemporaryFile(max_size=2048,mode='w+')

        try:
            #import pdb ; pdb.set_trace()
            cp = subprocess.run(
                [ 'powershell.exe', self.exe, ] + self.args,
                encoding=sys.stdout.encoding,
                stdout=outObj,
                stderr=errObj,
                timeout=self.timeout,
                check=True
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
