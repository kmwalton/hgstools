"""Classes and methods for various HGS post-processing calculations


Two loggers are available for inspecting this module's operations:

    import logging
    import pyhgs.calcs
    logger1 = logging.getLogger('pyhgs.calcs')
    logger2 = logging.getLogger('pyhgs.calcs_perf')

Stream handlers with custom formats may be attached to these loggers:

    fmtr = logging.Formatter('PYHGSCALC %(levelname)s\t%(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmtr)
    logger2.addHandler(ch)


"""

__docformat__ = 'numpy'

from .aabbox import (AABBox,)
__all__ = ['AABBox',]


# guard against these modules being unavailable
try:
    from .avconc import (AvCalc,)
    __all__.extend(['AvCalc',])
except ImportError:
    pass

try:
    from .avreggrid import (AvRegGrid,)
    __all__.extend(['AvRegGrid',])
except ImportError:
    pass

try:
    from .discharge import (DischargeCalc,)
    __all__.extend(['DischargeCalc',])
except ImportError:
    pass

