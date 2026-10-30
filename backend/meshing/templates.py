# Copyright (C) 2017-2025  Michael Steel.
# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import json
import sys
from importlib import import_module
from pathlib import Path

from numbat_core import reporting
import traceback

g_waveguide_templates = {}   # module level dictionary of waveguide template name to implementing class
_registered_waveguide_directories = []
_waveguide_manifest_filename = "numbat_waveguides.json"


def _read_waveguide_manifest(manifest_path):
    try:
        with manifest_path.open() as manifest_file:
            manifest = json.load(manifest_file)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid waveguide manifest: {manifest_path}") from exc

    if not isinstance(manifest, dict) or not isinstance(manifest.get("wguides"), list):
        raise ValueError(
            f"Waveguide manifest {manifest_path} must contain a 'wguides' list."
        )
    for index, waveguide in enumerate(manifest["wguides"]):
        if not isinstance(waveguide, dict):
            raise ValueError(
                f"Waveguide manifest {manifest_path} entry {index} must be an object."
            )
        missing_fields = {"module", "wg_class", "inc_shape"} - waveguide.keys()
        if missing_fields:
            missing = ", ".join(sorted(missing_fields))
            raise ValueError(
                f"Waveguide manifest {manifest_path} entry {index} is missing: {missing}."
            )
        if not all(
            isinstance(waveguide[field], str) and waveguide[field]
            for field in ("module", "wg_class", "inc_shape")
        ):
            raise ValueError(
                f"Waveguide manifest {manifest_path} entry {index} must give non-empty "
                "'module', 'wg_class', and 'inc_shape' strings."
            )
    return manifest


def register_waveguide_directory(directory):
    """Register a directory containing a NumBAT waveguide manifest.

    The directory is made available as a Python import root. Its active
    waveguides are loaded when ``NumBATApp`` initializes.
    """
    waveguide_directory = Path(directory).expanduser().resolve()
    if not waveguide_directory.is_dir():
        raise NotADirectoryError(
            f"Waveguide directory does not exist or is not a directory: {waveguide_directory}"
        )

    manifest_path = waveguide_directory / _waveguide_manifest_filename
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Waveguide directory {waveguide_directory} does not contain "
            f"{_waveguide_manifest_filename}"
        )

    _read_waveguide_manifest(manifest_path)

    if waveguide_directory not in _registered_waveguide_directories:
        _registered_waveguide_directories.append(waveguide_directory)
    if str(waveguide_directory) not in sys.path:
        sys.path.insert(0, str(waveguide_directory))


def _load_waveguide_templates(pth_wgtemplate_dir, pth_wgtemplate_index):
    """Loads and instantiates waveguide templates specified in .json file {pth_wgtemplate_index} in directory {pth_wgtemplate_dir}.

    Returns an index of waveguide inc_shape to class."""

    wg_index = _read_waveguide_manifest(Path(pth_wgtemplate_index))["wguides"]

    for wg in wg_index: # eg rectangular, circular, Onion, etc

        if not wg.get('active', 1):  # this mesh is turned off for now
            continue

        wgclsnm = wg['wg_class']  # Python class implementing this template
        module_name = wg["module"]
        try:
            py_mod = import_module(module_name)
        except Exception as ex:
            reporting.report_and_exit(
                f"Python couldn't import waveguide template module '{module_name}'."
                "\n Your code likely contains a syntax error or an attempt to load another "
                "module that was not successful."
                f"\n\n The Python traceback contained the error message:\n   {str(ex)}."
                f"\n\n\n The full traceback was {traceback.format_exc()}"
            )
        template_source = f"module: {module_name}"

        # Now extract the class that implements this waveguide template
        if not hasattr(py_mod, wgclsnm):
            reporting.report_and_exit(
                f"Can't find waveguide template class {wgclsnm} in {template_source}"
            )
        wgcls = getattr(py_mod, wgclsnm)

        wg['wg_template_cls'] = wgcls

    return wg_index


g_waveguide_templates = {}   # module level dictionary of waveguide template name to implementing class


def _add_waveguide_templates(waveguides, source):
    existing_shapes = {waveguide["inc_shape"] for waveguide in g_waveguide_templates}
    for waveguide in waveguides:
        inc_shape = waveguide["inc_shape"]
        if inc_shape in existing_shapes:
            reporting.report_and_exit(
                f"Waveguide template shape '{inc_shape}' from {source} is already registered."
            )
        existing_shapes.add(inc_shape)
    g_waveguide_templates.extend(waveguides)


# called at startup from NumBATApp.__init__
def initialise_waveguide_templates(nbapp):
    global g_waveguide_templates
    g_waveguide_templates = []
    pmsh_dir = nbapp.path_mesh_templates()

    pmsh_index_builtin = pmsh_dir / 'builtin_waveguides.json'
    if pmsh_index_builtin.exists():
        _add_waveguide_templates(
            _load_waveguide_templates(pmsh_dir, pmsh_index_builtin), pmsh_index_builtin
        )
    else:
        reporting.register_warning(
            f"Couldn't find builtin waveguide template index file: {pmsh_index_builtin}")

    for waveguide_directory in _registered_waveguide_directories:
        manifest_path = waveguide_directory / _waveguide_manifest_filename
        _add_waveguide_templates(
            _load_waveguide_templates(waveguide_directory, manifest_path), manifest_path
        )



def get_waveguide_template_class(inc_shape):

    for wg in g_waveguide_templates:
        if inc_shape in wg['inc_shape']:  # is the desired shape supported by this template class?

            wg_cls = wg['wg_template_cls']
            return wg_cls

    else:  # didn't find required wg
        reporting.report_and_exit(f"Selected inc_shape = '{inc_shape}' "
                                      'is not currently implemented. \nPlease make a mesh with gmsh and '
                                      'consider contributing this to NumBAT via github.')



def print_waveguide_help(inc_shape):
    for wg in g_waveguide_templates:
        if inc_shape in wg['inc_shape']:  # is the desired shape supported by this template class?
            #found = True

            # Instantiate the class that defines this waveguide model
            wg_geom = wg['wg_template_cls'](None, None)
            wg_geom.init_geometry()
            print(wg_geom.get_parameter_help_summary())