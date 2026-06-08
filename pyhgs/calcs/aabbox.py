"""Deprecated location for `AABBox`.

`AABBox` was promoted to the top-level `pyhgs.aabbox` module so it can serve as
the canonical *blockspec* type across `pyhgs`. This module re-exports it for
backward compatibility; import from `pyhgs.aabbox` in new code.
"""

from ..aabbox import AABBox

__all__ = ['AABBox']
