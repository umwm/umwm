"""Redirection of the standard output file descriptor.

The model writes its messages through the Fortran runtime, bypassing Python's
``sys.stdout``, so they can only be captured or silenced at the level of the
file descriptor. This works in scripts, under ``mpiexec``, and in Jupyter
kernels, whose own standard output is a separate stream.
"""

from __future__ import annotations

import contextlib
import os
import sys
from pathlib import Path

__all__ = ["redirect_stdout_fd"]


@contextlib.contextmanager
def redirect_stdout_fd(target, flush=None):
    """Send everything written to the standard output file descriptor during
    the block to ``target``: a path (the file is appended to), an open file
    object, or ``None`` to discard it. ``flush``, if given, is called before
    the descriptor is restored, to flush buffers of other runtimes."""
    sys.stdout.flush()
    if target is None:
        handle = open(os.devnull, "ab")
        close = True
    elif isinstance(target, (str, os.PathLike)):
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        handle = open(target, "ab")
        close = True
    else:
        handle = target
        close = False
        handle.flush()
    saved = os.dup(1)
    try:
        os.dup2(handle.fileno(), 1)
        yield
    finally:
        if flush is not None:
            flush()
        try:
            sys.stdout.flush()
        except (OSError, ValueError):
            pass
        os.dup2(saved, 1)
        os.close(saved)
        if close:
            handle.close()
