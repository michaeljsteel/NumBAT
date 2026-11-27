"""Control stdout around calls to the compiled Fortran backend."""

import contextlib
import os
from threading import RLock

from fortran import nb_fortran
from numtools.output import Verbosity, get_verbosity

from .logcapture import _flush_all, output_capture_active

_fortran_output_lock = RLock()


def _flush_fortran_stdout():
    try:
        flush = nb_fortran.flush_fortran_stdout
    except AttributeError as exc:
        raise RuntimeError(
            "The Fortran backend must be rebuilt to support quiet output."
        ) from exc
    status = flush()
    if status:
        raise OSError(f"Could not flush Fortran stdout (status {status}).")


@contextlib.contextmanager
def fortran_output():
    """Suppress process stdout in quiet mode unless a NumBAT log is active."""
    with _fortran_output_lock:
        if get_verbosity() != Verbosity.QUIET.value:
            yield
            return

        if output_capture_active():
            try:
                yield
            finally:
                _flush_fortran_stdout()
            return

        _flush_all()
        _flush_fortran_stdout()
        saved_stdout = os.dup(1)
        try:
            with open(os.devnull, "w") as null:
                os.dup2(null.fileno(), 1)
                try:
                    yield
                finally:
                    try:
                        _flush_fortran_stdout()
                        _flush_all()
                    finally:
                        os.dup2(saved_stdout, 1)
        finally:
            os.close(saved_stdout)


def call_fortran(function, *args, **kwargs):
    """Call a Fortran entry point with the current output policy."""
    with fortran_output():
        return function(*args, **kwargs)
