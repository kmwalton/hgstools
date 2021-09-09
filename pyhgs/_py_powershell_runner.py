"""HGS Runner for Window Powershell. PROVISIONAL

This implementation will be removed when the pyhgsrunner is operational
"""
import sys
import os
import re
import glob
import subprocess
import shutil
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
        self.captureRunPreexisting()

    def captureRunPreexisting(self):
        self.preexistingfiles = set(os.listdir(self.simdir))

    def eraseRunOutputs(self, keep=[], keep_file=None):
        """Erase simulation output files.

        This erases any files *in the simulation directory* or subdirectory that
        did not exist prior to the creation of this object or at the last call
        to `captureRunPreexisting()`.

        Arguments:
            keep : list-like
                A list of file names to keep; these will not be deleted, even if
                they did not exist when this object was instantiated.

            keep_file : str
                Name of a file that lists filenames or glob-style expressions
                (on separate lines) of additional files to keep, or None (or a
                non-readable filename) if no keep files are to be
                specified in this manner. 
        """
        keep = set(keep)
        if keep_file and os.path.isfile(keep_file):
            with open(keep_file,'r') as fin:
                for l in fin.readlines():
                    l = l.strip()
                    # check if the filename actually looks like a glob
                    if re.search(r'[*[?]',l):
                        keep.update( fn[len(self.simdir):]
                                for fn in glob.iglob(os.path.join(simdir,l)))
                    else:
                        keep.update(l)

        keep.update(self.preexistingfiles)

        # scan for items in the keep list that have a directory component.
        keep.update( d.split(os.sep)[0] for d in keep if os.path.dirname(d) )

        fullsimdir = os.path.abspath(simdir)

        for f in sorted(set(os.listdir(simdir)) - keep):
            fullf = os.path.abspath(os.path.join(self.simdir,f))
            if os.path.commonprefix(fullf,fullsimdir) != fullsimdir:
                # file/directory not rooted below fsimdir
                pass
            elif os.path.isdir(ff):
                shutil.rmtree(ff)
            elif os.path.isfile(ff):
                os.remove(ff)
    

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
                #stdout=outObj,
                #stderr=errObj,
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
