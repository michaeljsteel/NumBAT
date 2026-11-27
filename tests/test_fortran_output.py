"""Verify quiet Fortran output and restoration in isolated processes."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"


def _run(code):
    env = dict(os.environ, PYTHONPATH=str(BACKEND_DIR))
    return subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        cwd=BACKEND_DIR,
        env=env,
        timeout=120,
        check=False,
    )


def test_quiet_output_capture_nesting_and_restoration(tmp_path):
    result = _run(
        f"""
import ctypes
import os
from types import SimpleNamespace
from numbat.utils import fortran_output as fo
from numbat.utils.logcapture import redirect_output_to_file, output_capture_active
from numtools.output import set_verbosity

libc = ctypes.CDLL(None)
fo.nb_fortran = SimpleNamespace(flush_fortran_stdout=lambda: 0)
def emit():
    os.write(1, b"native stdout\\n")
    os.write(2, b"native stderr\\n")
    libc.printf(b"buffered stdout\\n")
    return 17

set_verbosity("normal")
assert fo.call_fortran(emit) == 17
set_verbosity("quiet")
assert fo.call_fortran(emit) == 17
try:
    def fail():
        os.write(1, b"hidden exception output\\n")
        raise RuntimeError("expected")
    fo.call_fortran(fail)
except RuntimeError as exc:
    assert str(exc) == "expected"
else:
    raise AssertionError("exception not propagated")
os.write(1, b"stdout restored\\n")
with redirect_output_to_file({str(tmp_path / "outer")!r}):
    assert output_capture_active()
    fo.call_fortran(emit)
    with redirect_output_to_file({str(tmp_path / "inner")!r}):
        assert output_capture_active()
        fo.call_fortran(emit)
    assert output_capture_active()
    fo.call_fortran(emit)
assert not output_capture_active()
fo.call_fortran(emit)
os.write(1, b"done\\n")
"""
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.count("native stdout") == 1
    assert result.stdout.count("buffered stdout") == 1
    assert "hidden exception output" not in result.stdout
    assert "stdout restored" in result.stdout
    assert "done" in result.stdout
    assert result.stderr.count("native stderr") == 3
    outer = (tmp_path / "outer_log.txt").read_text()
    inner = (tmp_path / "inner_log.txt").read_text()
    assert outer.count("native stdout") == 2
    assert outer.count("buffered stdout") == 2
    assert outer.count("native stderr") == 2
    assert inner.count("native stdout") == 1
    assert inner.count("buffered stdout") == 1
    assert inner.count("native stderr") == 1


def test_missing_flush_entry_point_and_flush_failure():
    result = _run(
        """
import os
from types import SimpleNamespace
from numbat.utils import fortran_output as fo
from numtools.output import set_verbosity

set_verbosity("quiet")
fo.nb_fortran = SimpleNamespace()
try:
    fo.call_fortran(lambda: None)
except RuntimeError as exc:
    assert "rebuilt" in str(exc)
else:
    raise AssertionError("missing flush not reported")
fo.nb_fortran = SimpleNamespace(flush_fortran_stdout=lambda: 5)
try:
    fo.call_fortran(lambda: None)
except OSError as exc:
    assert "status 5" in str(exc)
else:
    raise AssertionError("flush failure not reported")
os.write(1, b"restored\\n")
"""
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == "restored\n"


def test_compiled_fortran_fatal_stdout_is_quiet():
    probe = _run(
        "from fortran import nb_fortran; "
        "print(int(hasattr(nb_fortran, 'flush_fortran_stdout')))"
    )
    assert probe.returncode == 0, probe.stderr
    if probe.stdout.strip() != "1":
        pytest.skip("rebuild the Fortran backend to verify native output")

    code = """
import numpy as np
from fortran import nb_fortran
from numbat.utils.fortran_output import call_fortran
from numtools.output import set_verbosity

set_verbosity(VERBOSITY)
nodes = np.arange(1, 7, dtype=np.int64).reshape(6, 1)
coords = np.zeros((2, 6))
beta = np.ones(1, dtype=np.complex128)
field = np.zeros((3, 13, 1, 1), dtype=np.complex128)
call_fortran(nb_fortran.h_mode_field_ez, 1., 1, 1, 6, nodes, coords, beta, field)
print("unexpected return")
"""
    normal = _run(code.replace("VERBOSITY", "'normal'"))
    quiet = _run(code.replace("VERBOSITY", "'quiet'"))
    assert "H_mode_field_Ez" in normal.stdout
    assert "H_mode_field_Ez" not in quiet.stdout
    assert "unexpected return" not in normal.stdout
    assert "unexpected return" not in quiet.stdout
