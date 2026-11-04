"""Import checks for the supported backend namespace.

Each check runs in a fresh interpreter with PYTHONPATH=backend. Public and
internal implementation modules must import successfully.
"""

import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"

TUTORIAL_IMPORT_ORDER = [
    "numbat",
    "numbat.materials",
    "numbat.geometry",
    "numbat.mesh",
    "numbat.solvers",
    "numbat.plotting",
]

IMPORTABLE_AFTER_NUMBAT = [
    "numbat",
    "numbat_core.mesh",
    "numbat_core.mesh.fem",
    "numbat_core.mesh.mail",
    "numbat_core.mesh.templates",
    "numbat_core.mesh.builtin",
    "numbat_core.plotting",
    "numbat_core.plotting.axes",
    "numbat_core.plotting.fields",
    "numbat_core.plotting.gain",
    "numbat_core.plotting.materials",
    "numbat_core.plotting.mesh",
    "numbat_core.plotting.modes",
    "numbat_core.plotting.preferences",
    "numbat_core.plotting.profiles",
    "numbat_core.plotting.raw_modes",
    "numbat_core.plotting.tools",
    "numbat_core.geometry",
    "numbat_core.geometry.structure",
    "numbat_core.geometry.usermesh",
    "numbat_core.solvers",
    "numbat_core.solvers.modes",
    "numbat_core.solvers.modecalcs",
    "numbat_core.solvers.sbsgain",
    "numbat_core.solvers.simresult",
    "numbat_core.api",
    "numbat_core.api.application",
]

IMPORTABLE_STANDALONE = [
    "numbat",
    "numbat_core.mesh",
    "numbat_core.geometry",
    "numbat_core.geometry.structure",
    "numbat_core.geometry.usermesh",
    "numbat_core.solvers",
    "numbat_core.solvers.modes",
    "numbat_core.solvers.modecalcs",
    "numbat_core.solvers.sbsgain",
    "numbat_core.solvers.simresult",
    "numbat_core.api",
    "numbat_core.api.application",
    "numbat_core.mesh.fem",
    "numbat_core.mesh.mail",
    "numbat_core.mesh.templates",
    "numbat_core.mesh.builtin",
    "numbat_core.plotting",
    "numbat_core.plotting.axes",
    "numbat_core.plotting.fields",
    "numbat_core.plotting.gain",
    "numbat_core.plotting.materials",
    "numbat_core.plotting.mesh",
    "numbat_core.plotting.modes",
    "numbat_core.plotting.preferences",
    "numbat_core.plotting.profiles",
    "numbat_core.plotting.raw_modes",
    "numbat_core.plotting.tools",
]

# Entry points that may import the application root directly.
APP_ROOT_MODULES = {"numbat.py", "nb_install_tester.py"}


def _run_imports(*modules, cwd):
    env = {**os.environ, "PYTHONPATH": str(BACKEND_DIR)}
    code = "\n".join(f"import {module}" for module in modules)
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )


def _assert_imports_succeed(result):
    assert result.returncode == 0, result.stderr


def test_tutorial_import_sequence_succeeds_in_one_interpreter(tmp_path):
    _assert_imports_succeed(_run_imports(*TUTORIAL_IMPORT_ORDER, cwd=tmp_path))


@pytest.mark.parametrize("module", IMPORTABLE_AFTER_NUMBAT)
def test_module_imports_after_numbat(tmp_path, module):
    _assert_imports_succeed(_run_imports("numbat", module, cwd=tmp_path))


@pytest.mark.parametrize("module", IMPORTABLE_STANDALONE)
def test_module_imports_standalone(tmp_path, module):
    _assert_imports_succeed(_run_imports(module, cwd=tmp_path))


def _imports_application_root(path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any(alias.name == "numbat" for alias in node.names):
                return True
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module == "numbat":
            return True
    return False


def test_low_level_modules_do_not_import_application_root():
    offenders = [
        str(path.relative_to(BACKEND_DIR))
        for path in sorted(BACKEND_DIR.rglob("*.py"))
        if path.name not in APP_ROOT_MODULES
        and _imports_application_root(path)
    ]
    assert offenders == [], f"low-level modules import numbat: {offenders}"


def test_mesh_modules_do_not_import_internal_plotting():
    mesh_dir = BACKEND_DIR / "numbat_core" / "mesh"
    offenders = []
    for path in sorted(mesh_dir.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith(
                "numbat_core.plotting"
            ):
                offenders.append(str(path.relative_to(BACKEND_DIR)))
                break
            if isinstance(node, ast.Import) and any(
                alias.name.startswith("numbat_core.plotting") for alias in node.names
            ):
                offenders.append(str(path.relative_to(BACKEND_DIR)))
                break
    assert offenders == [], f"mesh modules import plotting: {offenders}"
