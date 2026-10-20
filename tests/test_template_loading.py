import json
import pickle
import sys
from pathlib import Path

import pytest

sys.path.append("backend")

from meshing import template_loader
from meshing.templates import _load_waveguide_templates


def _unload_module(module_path):
    resolved_path = Path(module_path).resolve()
    module_name = template_loader._module_name_for_path(resolved_path)
    template_loader._loaded_modules.pop(resolved_path, None)
    sys.modules.pop(module_name, None)


@pytest.fixture(autouse=True)
def clean_loaded_template_modules():
    original_modules = dict(template_loader._loaded_modules)
    yield
    for module_path in list(template_loader._loaded_modules):
        if module_path not in original_modules:
            _unload_module(module_path)


def test_load_waveguide_module_caches_resolved_path_and_registers_module(tmp_path):
    module_path = tmp_path / "template.py"
    module_path.write_text("execution_count = 1\n")

    first_module = template_loader.load_waveguide_module(module_path)
    second_module = template_loader.load_waveguide_module(module_path.resolve())

    assert first_module is second_module
    assert first_module.execution_count == 1
    assert first_module.__name__ == template_loader._module_name_for_path(module_path.resolve())
    assert sys.modules[first_module.__name__] is first_module


def test_load_waveguide_module_distinguishes_same_named_files_in_different_directories(tmp_path):
    first_path = tmp_path / "first" / "template.py"
    second_path = tmp_path / "second" / "template.py"
    first_path.parent.mkdir()
    second_path.parent.mkdir()
    first_path.write_text("value = 'first'\n")
    second_path.write_text("value = 'second'\n")

    first_module = template_loader.load_waveguide_module(first_path)
    second_module = template_loader.load_waveguide_module(second_path)

    assert first_module is not second_module
    assert first_module.__name__ != second_module.__name__
    assert first_module.value == "first"
    assert second_module.value == "second"


def test_load_waveguide_module_removes_failed_import_from_module_registry(tmp_path):
    module_path = tmp_path / "broken_template.py"
    module_path.write_text("raise RuntimeError('broken template')\n")
    module_name = template_loader._module_name_for_path(module_path.resolve())

    with pytest.raises(RuntimeError, match="broken template"):
        template_loader.load_waveguide_module(module_path)

    assert module_name not in sys.modules
    assert module_path.resolve() not in template_loader._loaded_modules


def test_waveguide_index_reuses_one_module_for_classes_from_same_implementation(tmp_path):
    module_path = tmp_path / "templates.py"
    module_path.write_text("class FirstTemplate:\n    pass\n\nclass SecondTemplate:\n    pass\n")
    index_path = tmp_path / "waveguides.json"
    index_path.write_text(
        json.dumps(
            {
                "wguides": [
                    {"wg_impl": "templates.py", "wg_class": "FirstTemplate", "inc_shape": "first"},
                    {
                        "wg_impl": "templates.py",
                        "wg_class": "SecondTemplate",
                        "inc_shape": "second",
                    },
                ]
            }
        )
    )

    waveguides = _load_waveguide_templates(tmp_path, index_path)

    first_class = waveguides[0]["wg_template_cls"]
    second_class = waveguides[1]["wg_template_cls"]
    assert first_class.__module__ == second_class.__module__
    assert first_class.__module__ == template_loader._module_name_for_path(module_path.resolve())


def test_pickle_rebuilds_file_based_geometry_with_its_saved_state(tmp_path):
    module_path = tmp_path / "geometry.py"
    module_path.write_text(
        "from usermesh import UserGeometryBase\n\n"
        "class ExampleGeometry(UserGeometryBase):\n"
        "    pass\n"
    )
    module = template_loader.load_waveguide_module(module_path)
    geometry = module.ExampleGeometry({"size": 3}, None)
    geometry.extra_state = "preserved"
    serialized_geometry = pickle.dumps(geometry)

    _unload_module(module_path)
    rebuilt_geometry = pickle.loads(serialized_geometry)

    assert type(rebuilt_geometry).__name__ == "ExampleGeometry"
    assert type(rebuilt_geometry).__module__ == template_loader._module_name_for_path(
        module_path.resolve()
    )
    assert rebuilt_geometry._d_params == {"size": 3}
    assert rebuilt_geometry.extra_state == "preserved"
