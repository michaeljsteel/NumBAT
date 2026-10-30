"""Phase 2 checks: leaf modules moved into ``numbat_core`` keep their old names.

Each flat module at the backend root is a facade over its ``numbat_core``
home. Every public name in the new module must resolve through the facade to
the same object, so existing ``import nbtypes`` style code keeps working.
Runs in a fresh interpreter with PYTHONPATH=backend.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"

FACADE_PAIRS = [
    ("nbtypes", "numbat_core.nbtypes"),
    ("reporting", "numbat_core.reporting"),
    ("nbversion", "numbat_core.nbversion"),
    ("nbplatform", "numbat_core.nbplatform"),
    ("numbattools", "numbat_core.utils.numeric"),
    ("crystalsyms", "numbat_core.materials.crystalsyms"),
    ("materialprops", "numbat_core.materials.materialprops"),
    ("bulkprops", "numbat_core.materials.bulkprops"),
    ("voigt", "numbat_core.materials.elastic_tensors"),
]


def _run(code: str) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(BACKEND_DIR)
    return subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=env,
        cwd=BACKEND_DIR,
        timeout=120,
    )


@pytest.mark.parametrize("old, new", FACADE_PAIRS)
def test_facade_exposes_new_module_public_names_by_identity(old, new):
    code = f"""
import importlib
old = importlib.import_module({old!r})
new = importlib.import_module({new!r})
public = [n for n in dir(new) if not n.startswith('_')]
missing = [n for n in public if not hasattr(old, n)]
assert not missing, f'facade missing names: {{missing}}'
mismatch = [n for n in public if getattr(old, n) is not getattr(new, n)]
assert not mismatch, f'facade names not identical: {{mismatch}}'
print('ok', len(public))
"""
    res = _run(code)
    assert res.returncode == 0, res.stderr
    assert res.stdout.startswith("ok")


@pytest.mark.parametrize("_old, new", FACADE_PAIRS)
def test_numbat_core_module_imports_standalone(_old, new):
    res = _run(f"import importlib; importlib.import_module({new!r})")
    assert res.returncode == 0, res.stderr
