
import time
from collections import deque

class _PerfLogStack:
    """`logging.Logger` object with FIFO internal datastructure

    Geared towards measuring performance of nested levels of tasks.
    """
    def __init__(self, logfunc):
        self.perfstk = deque()
        self.logfunc = logfunc
    def push(self, msg=None):
        if msg is not None:
            self.logfunc(msg)
        self.perfstk.append(time.perf_counter())
    def pop(self, msg=None):
        if not msg:
            msg = 'Last operation'
        self.logfunc(
            msg
            +f' took {time.perf_counter()-self.perfstk.pop():.2f}s')
        return self

