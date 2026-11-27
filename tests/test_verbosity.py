"""Tests for NumBAT's configurable Python output."""

import os
import subprocess
import sys
from pathlib import Path

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"


def _run(args):
    env = dict(os.environ, PYTHONPATH=str(BACKEND_DIR))
    return subprocess.run(
        [sys.executable, *args],
        capture_output=True,
        text=True,
        cwd=BACKEND_DIR,
        env=env,
        timeout=120,
        check=False,
    )


def test_verbosity_and_print_semantics(tmp_path):
    code = f"""
import contextlib
import io
import numbat
from numbat import reporting
from numbat.plotting.ipython import img_single
from numbat.mesh.fem import FemMesh
from nbanalytic.elasticstack import ElasticStack
from types import SimpleNamespace

out = io.StringIO()
with contextlib.redirect_stdout(out):
    app = numbat.NumBATApp(outdir={str(tmp_path)!r},
                           ignore_user_settings=True, verbosity="quiet")
    assert app.verbosity() == "quiet"
    numbat.nbprint("hidden")
    print("user output")
    reporting.register_warning("visible warning")
    img_single("nonexistent-image.png")
assert out.getvalue() == ("user output\\nWarning: visible warning\\n"
                          "Can't find image file nonexistent-image.png.\\n")

mesh = SimpleNamespace(num_nodes=6, num_elements=1)
props = SimpleNamespace(n_mats_em=1, v_refindexn=[1.5])
material = SimpleNamespace(material_name="Glass")
guide = SimpleNamespace(optical_props=props, get_optical_materials=lambda: [material])
out = io.StringIO()
with contextlib.redirect_stdout(out):
    FemMesh._report_em_mesh_properties(mesh, guide)
    assert ElasticStack._filter_determinant_minima(None, [(1., 1.)], 0.1, 1., True) == []
assert out.getvalue() == ""

app.set_verbosity("normal")
out = io.StringIO()
with contextlib.redirect_stdout(out):
    FemMesh._report_em_mesh_properties(mesh, guide)
    ElasticStack._filter_determinant_minima(None, [(1., 1.)], 0.1, 1., True)
    numbat.nbprint("one", "two", sep=":", end="!", flush=True)
    numbat.nbprint("detailed hidden", level="detailed")
assert "final EM sim mesh has 6 nodes" in out.getvalue()
assert "Glass" in out.getvalue()
assert "Dropping 1 modes" in out.getvalue()
assert out.getvalue().endswith("one:two!")

app.set_verbosity("detailed")
out = io.StringIO()
numbat.nbprint("details", file=out, level="detailed")
assert out.getvalue() == "details\\n"
for invalid in ("silent", "", None):
    try:
        app.set_verbosity(invalid)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid verbosity accepted")
assert app.verbosity() == "detailed"
try:
    numbat.nbprint("bad", level="silent")
except ValueError:
    pass
else:
    raise AssertionError("invalid output level accepted")
"""
    result = _run(["-c", code])
    assert result.returncode == 0, result.stderr


def test_standalone_analytic_import_does_not_load_numbat():
    result = _run([
        "-c",
        "import sys; import nbanalytic.elasticcylinder, numtools.optimize; "
        "assert 'numbat' not in sys.modules",
    ])
    assert result.returncode == 0, result.stderr


def test_temporary_verbosity_context(tmp_path):
    code = f"""
import contextlib
import io
import numbat

app = numbat.NumBATApp(outdir={str(tmp_path)!r}, ignore_user_settings=True)
assert app.verbosity() == "normal"
scope = app.verbosity_at()
assert app.verbosity() == "normal"
out = io.StringIO()
with contextlib.redirect_stdout(out):
    numbat.nbprint("before")
    with scope:
        assert app.verbosity() == "quiet"
        numbat.nbprint("hidden")
        numbat.nbprint("warning", level=numbat.OutputLevel.ALWAYS)
        with app.verbosity_at(numbat.Verbosity.DETAILED):
            assert app.verbosity() == "detailed"
            numbat.nbprint("detail", level=numbat.OutputLevel.DETAILED)
        assert app.verbosity() == "quiet"
        with app.verbosity_at("normal"):
            numbat.nbprint("nested")
        assert app.verbosity() == "quiet"
    numbat.nbprint("after")
assert out.getvalue() == "before\\nwarning\\ndetail\\nnested\\nafter\\n"
assert app.verbosity() == "normal"

try:
    with app.verbosity_at("detailed"):
        app.set_verbosity("quiet")
        raise RuntimeError("expected")
except RuntimeError as exc:
    assert str(exc) == "expected"
else:
    raise AssertionError("exception not propagated")
assert app.verbosity() == "normal"

try:
    with app.verbosity_at("invalid"):
        raise AssertionError("invalid scope entered")
except ValueError:
    pass
else:
    raise AssertionError("invalid verbosity accepted")
assert app.verbosity() == "normal"
"""
    result = _run(["-c", code])
    assert result.returncode == 0, result.stderr


def test_output_enums_and_string_compatibility(tmp_path):
    code = f"""
import contextlib
import io
import numbat
from numtools.output import OutputLevel, Verbosity, set_verbosity

assert numbat.OutputLevel is OutputLevel
assert numbat.Verbosity is Verbosity
app = numbat.NumBATApp(outdir={str(tmp_path)!r}, ignore_user_settings=True,
                       verbosity=Verbosity.QUIET)
assert app.verbosity() == "quiet"

expected = {{
    Verbosity.QUIET: {{OutputLevel.ALWAYS}},
    Verbosity.NORMAL: {{OutputLevel.ALWAYS, OutputLevel.NORMAL}},
    Verbosity.DETAILED: set(OutputLevel),
}}
for verbosity in Verbosity:
    for selected_verbosity in (verbosity, verbosity.value):
        app.set_verbosity(selected_verbosity)
        assert app.verbosity() == verbosity.value
        for level in OutputLevel:
            for selected_level in (level, level.value):
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    numbat.nbprint("message", level=selected_level)
                assert out.getvalue() == ("message\\n" if level in expected[verbosity] else "")

for invalid in (None, "bad", Verbosity.QUIET, Verbosity.NORMAL):
    try:
        numbat.nbprint("bad", level=invalid)
    except ValueError:
        pass
    else:
        raise AssertionError(f"invalid output level accepted: {{invalid!r}}")
for invalid in (None, "bad", OutputLevel.NORMAL, OutputLevel.DETAILED):
    try:
        set_verbosity(invalid)
    except ValueError:
        pass
    else:
        raise AssertionError(f"invalid verbosity accepted: {{invalid!r}}")
assert app.verbosity() == "detailed"
"""
    result = _run(["-c", code])
    assert result.returncode == 0, result.stderr


def test_reporting_warnings_and_fatal_errors_in_quiet_mode():
    code = """
import contextlib
import io
import logging
from numbat import reporting

reporting.set_verbosity("quiet")
records = []
class RecordHandler(logging.Handler):
    def emit(self, record):
        records.append(record)
logger = logging.getLogger("NumBAT log")
logger.addHandler(RecordHandler())
logger.setLevel(logging.WARNING)

stdout, stderr = io.StringIO(), io.StringIO()
with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
    reporting.report("incomplete parameters")
    reporting.register_warning("check the modes")
    reporting.assertion(False, "nonfatal", severity="warning")
    try:
        reporting.report_and_exit("invalid input")
    except SystemExit as exc:
        assert exc.code == 1
    else:
        raise AssertionError("fatal error did not exit")

assert stdout.getvalue() == "Warning: check the modes\\n"
assert stderr.getvalue() == ("Warning: incomplete parameters\\n"
                             "Warning: NumBAT assertion failed: nonfatal\\n"
                             "\\nFatal error: \\n invalid input\\n\\n\\n")
assert "check the modes" in reporting.report_warnings()
assert "Warning: check the modes" not in reporting.report_warnings()
assert records == []
"""
    result = _run(["-c", code])
    assert result.returncode == 0, result.stderr


def test_verbosity_reaches_spawned_workers(tmp_path):
    script = tmp_path / "worker.py"
    script.write_text(
        """
import multiprocessing
import numbat
from numbat.utils.threading import launch_worker_processes_and_wait

def task(value):
    numbat.nbprint("worker normal output")
    numbat.nbprint("worker warning", level="always")
    return value * 2

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    app = numbat.NumBATApp(ignore_user_settings=True, verbosity="quiet")
    work = multiprocessing.JoinableQueue()
    result = multiprocessing.Queue()
    work.put(3)
    launch_worker_processes_and_wait(1, task, result, work)
    assert result.get(timeout=5) == 6
    print("user output")
"""
    )
    result = _run([str(script)])
    assert result.returncode == 0, result.stderr
    assert "worker normal output" not in result.stdout
    assert "worker warning" in result.stdout
    assert "user output" in result.stdout
    assert "Assigning" not in result.stdout
    assert "exhausted" not in result.stdout
