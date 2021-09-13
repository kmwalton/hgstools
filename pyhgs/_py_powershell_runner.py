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
from tempfile import SpooledTemporaryFile
from itertools import chain

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

    @staticmethod
    def read_keep_file(keep_file):
        """Return a list (string or unexpanded glob filenames) from keep_file.

        Arguments:

            keep_file : str
                Name of a file that lists filenames or glob-style expressions
                (on separate lines; one per line) of additional files to keep.
                Lines may be blank. Python-style comments are allowed.
                File names or expressions with spaces must be quoted; file names
                with a '#' character must be escaped by '\#', else the balance
                of the line will be treated as a comment.

                If no keep files are to be specified in this manner, pass None
                (default) or a string representing a non-existant or non-readable
                filename.
        """
        exprs = []
        with open(keep_file,'r') as fin:
            for l in fin.readlines():

                # split off comments
                l = shlex.split(l,comments=True)

                if not l:
                    # ignore blank lines or comment-only lines
                    continue
                else:
                    # grab just the first argument of the line --- users
                    # bewrare that if you need a space in the expression
                    # then it will need to be quoted.
                    l = l[0]

                exprs.append(l)

        return exprs

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
                See `PyPowerShellRunner.read_keep_file`.
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
            keep_set.update( chain.from_iterable(deglob(x)
                for x in PyPowerShellRunner.read_keep_file(keep_file)) )

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
