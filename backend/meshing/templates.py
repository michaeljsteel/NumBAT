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
import importlib

import reporting
import traceback

g_waveguide_templates = {}   # module level dictionary of waveguide template name to implementing class


def _load_waveguide_templates(pth_wgtemplate_dir, pth_wgtemplate_index):
    """Loads and instantiates waveguide templates specified in .json file {pth_wgtemplate_index} in directory {pth_wgtemplate_dir}.

    Returns an index of waveguide inc_shape to class."""

    try:
        with open(pth_wgtemplate_index) as fin:
            wg_index = json.load(fin)['wguides']
    except Exception as ex:
        raise Exception('JSON parse error in reading user mesh index file' + str(ex)) from ex

    for wg in wg_index: # eg rectangular, circular, Onion, etc

        if not wg.get('active', 1):  # this mesh is turned off for now
            continue

        wgnm = wg['inc_shape']    # Name of inc_shape requested by user
        wgclsnm = wg['wg_class']  # Python class implementing this template
        wgpy = wg['wg_impl']      # Python file defining that class

        # Find and load the module containing the class that implements this waveguide template
        pth_mod = pth_wgtemplate_dir / wgpy

        if not pth_mod.exists():
            reporting.report_and_exit(f'Missing waveguide template implementation file: {str(wgpy)} named in str({pth_wgtemplate_index}).'
                            + f'\nFile was expected in directory {str(pth_wgtemplate_dir)}')

        # wgnm is just a convenient name
        spec = importlib.util.spec_from_file_location(wgnm, pth_mod)

        try:
            py_mod = spec.loader.load_module()
        except Exception as ex:
            reporting.report_and_exit(f"Python couldn't load the user module '{str(wgpy)}'."+
                                      "\n Your code likely contains a syntax error or an attempt to load another module that was not successful." + f'\n\n The Python traceback contained the error message:\n   {str(ex)}.'+
                                      f'\n\n\n The full traceback was {traceback.format_exc()}')

        # Now extract the class that implements this waveguide template
        if not hasattr(py_mod, wgclsnm):
            reporting.report_and_exit(f"Can't find waveguide template class {wgclsnm} in implementation file: {wgpy}")
        wgcls = getattr(py_mod, wgclsnm)

        wg['wg_template_cls'] = wgcls

    return wg_index


g_waveguide_templates = {}   # module level dictionary of waveguide template name to implementing class

# called at startup from NumBATApp.__init__
def initialise_waveguide_templates(nbapp):
    global g_waveguide_templates
    pmsh_dir = nbapp.path_mesh_templates()

    pmsh_index_builtin = pmsh_dir / 'builtin_waveguides.json'
    pmsh_index_user = pmsh_dir / 'user_waveguides.json'

    if pmsh_index_builtin.exists():
        g_waveguide_templates = _load_waveguide_templates(pmsh_dir, pmsh_index_builtin)
    else:
        reporting.register_warning(
            f"Couldn't find builtin waveguide template index file: {pmsh_index_builtin}")

    if pmsh_index_user.exists():
        user_wg_templates =  _load_waveguide_templates(pmsh_dir, pmsh_index_user)
        g_waveguide_templates.extend(user_wg_templates)
    else:
        reporting.register_warning(
            f"Couldn't find user waveguide template index file: {pmsh_index_user}")




def get_waveguide_template_class(inc_shape):

    for wg in g_waveguide_templates:
        if inc_shape in wg['inc_shape']:  # is the desired shape supported by this template class?

            wg_cls = wg['wg_template_cls']
            return wg_cls

    else:  # didn't find required wg
        reporting.report_and_exit(f"Selected inc_shape = '{self.inc_shape}' "
                                      'is not currently implemented. \nPlease make a mesh with gmsh and '
                                      'consider contributing this to NumBAT via github.')