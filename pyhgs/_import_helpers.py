"""Help find the ofrac package!

Due to "legacy issues" the ofracs.py module was usually in the top-level of
segments of the PYTHONPATH environment variable. A switch to the normal
conditions, where ofracs.py is accessed from its package, like 'import
ofrac.ofracs' has left a lot of old code (or old environment setups) not
working. As a stop-gap, this is a helper module to transparently find things
within the ofrac.ofracs module and give warnings when the environment needs to
be updated.
"""
import sys
import importlib
import warnings

def _get_ofracs_module():
    """Return the ofracs module from wherever its found"""

    ofracs_module = None

    if 'ofrac.ofracs' in sys.modules:
        ofracs_module = sys.modules['ofrac.ofracs']
    elif 'ofracs' in sys.modules:
        ofracs_module = sys.modules['ofracs']
    else:
        try:
            ofracs_module = importlib.import_module('ofrac.ofracs')
        except ModuleNotFoundError as e:
            try:
                ofracs_module = importlib.import_module('ofracs')
            finally:
                warnings.warn('This environment includes a one-too-deep ' \
                    'path, path/to/package_libary/ofrac.')

    return ofracs_module

def _get_from_ofracs(*names):
    """Return a dict of names->class_function_or_object"""

    ofracs_module = _get_ofracs_module()

    retdict = {}

    for n in names:
        retdict[n] = getattr(ofracs_module, n)

    return retdict

def _find_OFracGrid():
    """Find/Return the class object 'OFracGrid'

    To mitigate issues in PYTHONPATH'ing of ofracs...
    """

    _ofracgrid_modules_search = [ 'ofrac.ofracs', 'ofracs', ]
    for _mod_pth in _ofracgrid_modules_search:
        try:
            _mod = importlib.import_module(_mod_pth)
        except ImportError:
            pass
        else:
            OFracGrid = getattr(_mod, 'OFracGrid')
            break

    if not 'OFracGrid' in dir():
        raise ImportError('Could not import OFracGrid from any in '
            + str(_ofracgrid_modules_search))

    return OFracGrid

