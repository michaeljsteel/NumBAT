"""Tests for scripts/nb_convert_imports.py."""

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "nb_convert_imports.py"

spec = importlib.util.spec_from_file_location("nb_convert_imports", SCRIPT)
conv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(conv)


def test_line_rules():
    cl = lambda s: conv.convert_line(s)[0]
    assert cl("import materials") == "from numbat import materials"
    assert cl("  import integration  # x") == "  from numbat.solvers import sbsgain as integration  # x"
    assert cl("from integration import GainProps") == "from numbat.solvers import GainProps"
    assert cl("import materials as m") == "from numbat import materials as m"
    assert cl("import plotgain") == "from numbat import plotting as plotgain"
    assert cl("from nbtypes import SI_GHz, twopi") == "from numbat.types import SI_GHz, twopi"
    assert cl("from nbthread import launch_worker_processes_and_wait").startswith("from numbat.parallel")
    assert cl("import numbat") == "import numbat"
    assert cl("import numpy as np") == "import numpy as np"
    assert conv.convert_line("import numbat, materials")[1] is True


def test_idempotent():
    src = "import numbat\nimport materials\nfrom nbtypes import SI_GHz\n"
    once, _ = conv.convert_source(src)
    twice, _ = conv.convert_source(once)
    assert once == twice


def test_notebook_and_cli(tmp_path):
    nb = {"cells": [{"cell_type": "code", "metadata": {}, "outputs": [], "execution_count": None,
                     "source": ["import numbat\n", "import materials\n"]}],
          "metadata": {}, "nbformat": 4, "nbformat_minor": 5}
    f = tmp_path / "a.ipynb"
    f.write_text(json.dumps(nb, indent=1) + "\n")
    res = subprocess.run([sys.executable, str(SCRIPT), "--write", str(f)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert json.loads(f.read_text())["cells"][0]["source"][1] == "from numbat import materials\n"
    assert (tmp_path / "a.ipynb.orig").exists()


def test_converted_imports_resolve():
    code = """
import numbat
from numbat import materials, plotting as plotgain, mesh as nbgmsh
from numbat.solvers import integration, modecalcs, modes, sbsgain, simresult
from numbat.geometry import structure
from numbat.types import SI_GHz
"""
    env = dict(os.environ, PYTHONPATH=str(ROOT / "backend"))
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, cwd=ROOT / "backend")
    assert res.returncode == 0, res.stderr
