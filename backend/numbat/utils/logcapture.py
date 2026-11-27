"""Capture stdout and stderr, including output from compiled (Fortran/C) code."""

import contextlib
import ctypes
import os
import sys


def _flush_all():
    for s in (sys.stdout, sys.stderr, sys.__stdout__, sys.__stderr__):
        try:
            s.flush()
        except Exception:
            pass
    try:
        ctypes.CDLL(None).fflush(None)  # C stdio buffers
    except Exception:
        pass


@contextlib.contextmanager
def redirect_output_to_file(log_prefix=None):
    """Redirect stdout and stderr (at file descriptor level) to <log_prefix>_log.txt.

    Because the redirection is done on fds 1 and 2, output written by Fortran
    and C code is captured as well as Python output. Original streams are
    restored on exit. If log_prefix is None, nothing is redirected.
    """
    if log_prefix is None:
        yield None
        return

    logname = f"{log_prefix}_log.txt"
    _flush_all()
    try:
        log_stat = os.stat(logname)
    except FileNotFoundError:
        log_stat = None
    if log_stat is not None and all(
        os.path.samestat(log_stat, os.fstat(fd)) for fd in (1, 2)
    ):
        yield logname
        return

    with contextlib.ExitStack() as stack:
        saved_out = os.dup(1)
        stack.callback(os.close, saved_out)
        saved_err = os.dup(2)
        stack.callback(os.close, saved_err)
        with open(logname, "w") as logf:
            try:
                os.dup2(logf.fileno(), 1)
                os.dup2(logf.fileno(), 2)
                yield logname
            finally:
                _flush_all()
                os.dup2(saved_out, 1)
                os.dup2(saved_err, 2)
