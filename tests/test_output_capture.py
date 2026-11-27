"""Exercise process-wide output capture in isolated interpreters."""

import os
import subprocess
import sys
from pathlib import Path

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"


def run_script(code):
    env = dict(os.environ, PYTHONPATH=str(BACKEND_DIR))
    return subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=env,
        cwd=BACKEND_DIR,
        timeout=120,
        check=False,
    )


def test_capture_and_restore_with_nested_capture(tmp_path):
    result = run_script(
        f"""
import os
import ctypes
import subprocess
import sys
import numbat
from numbat.utils.logcapture import redirect_output_to_file

app = numbat.NumBATApp("run", outdir={str(tmp_path)!r}, ignore_user_settings=True)
print("before")
assert app.start_output_capture() == {str(tmp_path / "run_log.txt")!r}
print("python stdout")
print("python stderr", file=sys.stderr)
os.write(1, b"native stdout\\n")
os.write(2, b"native stderr\\n")
ctypes.CDLL(None).printf(b"buffered C stdout\\n")
subprocess.run([sys.executable, "-c", "print('child output')"], check=True)
try:
    app.start_output_capture()
except RuntimeError:
    pass
else:
    raise AssertionError("starting capture twice should fail")
with redirect_output_to_file(app.outpath()):
    print("nested output")
print("after nested")
with redirect_output_to_file(app.outpath("separate")):
    print("separate output")
print("after separate")
app.stop_output_capture()
app.stop_output_capture()
print("after")
print("restored stderr", file=sys.stderr)
"""
    )
    assert result.returncode == 0, result.stderr
    log = (tmp_path / "run_log.txt").read_text()
    for text in (
        "python stdout",
        "python stderr",
        "native stdout",
        "native stderr",
        "buffered C stdout",
        "child output",
        "nested output",
        "after nested",
        "after separate",
    ):
        assert text in log
        assert text not in result.stdout
        assert text not in result.stderr
    assert "before" in result.stdout
    assert "after\n" in result.stdout
    assert "restored stderr" in result.stderr
    assert "restored stderr" not in log
    assert "separate output" not in log
    assert (tmp_path / "separate_log.txt").read_text() == "separate output\n"


def test_restart_custom_prefix_and_exit_flush(tmp_path):
    result = run_script(
        f"""
import sys
import numbat

app = numbat.NumBATApp("run", outdir={str(tmp_path)!r}, ignore_user_settings=True)
app.start_output_capture()
print("first run")
app.stop_output_capture()
app.start_output_capture()
print("second run")
app.stop_output_capture()
assert app.start_output_capture("custom") == {str(tmp_path / "custom_log.txt")!r}
print("buffered at exit", end="")
print("stderr at exit", file=sys.stderr)
"""
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "run_log.txt").read_text() == "second run\n"
    log = (tmp_path / "custom_log.txt").read_text()
    assert "buffered at exit" in log
    assert "stderr at exit" in log
    assert "buffered at exit" not in result.stdout


def test_failed_start_does_not_activate_capture(tmp_path):
    result = run_script(
        f"""
import os
import numbat

app = numbat.NumBATApp("run", outdir={str(tmp_path)!r}, ignore_user_settings=True)
try:
    app.start_output_capture("missing/path")
except FileNotFoundError:
    pass
else:
    raise AssertionError("missing directory should fail")
print("still on terminal")
original_dup2 = os.dup2
def fail_stderr_once(source, target):
    if target == 2:
        os.dup2 = original_dup2
        raise OSError("simulated redirection failure")
    return original_dup2(source, target)
os.dup2 = fail_stderr_once
try:
    app.start_output_capture()
except OSError:
    pass
else:
    raise AssertionError("redirection failure should propagate")
print("restored after failed redirection")
app.start_output_capture()
print("working capture")
app.stop_output_capture()
"""
    )
    assert result.returncode == 0, result.stderr
    assert "still on terminal" in result.stdout
    assert "restored after failed redirection" in result.stdout
    assert (tmp_path / "run_log.txt").read_text() == "working capture\n"
