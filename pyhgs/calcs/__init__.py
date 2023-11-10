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

from .avconc import (AvCalc, )
from .avreggrid import (AvRegGrid, )

__all__ = ['AvCalc', 'AvRegGrid',]
