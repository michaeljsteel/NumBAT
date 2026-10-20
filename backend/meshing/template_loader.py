# Copyright (C) 2017-2025  Michael Steel.
#
# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Load file-based waveguide template modules with stable identities."""

import hashlib
import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType

_loaded_modules = {}


def _module_name_for_path(module_path: Path) -> str:
    """Return a deterministic import name for a resolved template file."""
    safe_stem = re.sub(r"\W|^(?=\d)", "_", module_path.stem)
    path_digest = hashlib.sha256(str(module_path).encode("utf-8")).hexdigest()[:16]
    return f"_numbat_waveguide_{safe_stem}_{path_digest}"


def load_waveguide_module(module_file: str | Path) -> ModuleType:
    """Load a waveguide module once for each resolved implementation path."""
    module_path = Path(module_file).resolve()
    module_name = _module_name_for_path(module_path)
    cached_module = _loaded_modules.get(module_path)

    if cached_module is not None:
        registered_module = sys.modules.get(module_name)
        if registered_module is None:
            sys.modules[module_name] = cached_module
        elif registered_module is not cached_module:
            raise ImportError(
                f"Waveguide template module name collision for implementation file: {module_path}"
            )
        return cached_module

    registered_module = sys.modules.get(module_name)
    if registered_module is not None:
        registered_file = getattr(registered_module, "__file__", None)
        if registered_file is None or Path(registered_file).resolve() != module_path:
            raise ImportError(
                f"Waveguide template module name collision for implementation file: {module_path}"
            )
        _loaded_modules[module_path] = registered_module
        return registered_module

    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load waveguide template module: {module_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        if sys.modules.get(module_name) is module:
            del sys.modules[module_name]
        raise

    _loaded_modules[module_path] = module
    return module
