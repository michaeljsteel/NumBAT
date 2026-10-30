"""Generated mesh files go to a scratch folder in the working directory.

The scratch folder is created under the current working directory, not inside
the backend source tree, and is reused across calls.
"""

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat_core import runtime  # noqa: E402


def test_scratch_dir_is_created_in_working_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    scratch = runtime.path_scratch()

    assert scratch == tmp_path / "numbat_tmp"
    assert scratch.is_dir()


def test_scratch_dir_is_reused(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    marker = runtime.path_scratch() / "marker.txt"
    marker.write_text("keep")

    assert runtime.path_scratch() == tmp_path / "numbat_tmp"
    assert marker.read_text() == "keep"


def test_structure_does_not_write_into_source_tree():
    structure_src = (Path(__file__).resolve().parents[1] / "backend" / "structure.py").read_text()
    assert "/ 'build'" not in structure_src
    assert "runtime.path_scratch()" in structure_src
