# objects.py is a subroutine of NumBAT. It contains the Struct
# objects that represent the structure being simulated.

# Copyright (C) 2017  Bjorn Sturmberg, Kokou Dossou.

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


# import time
import os
import subprocess
import copy
import traceback


import tempfile
from pathlib import Path
import json
import importlib

import numpy as np



from mode_calcs import Simulation
import numbat
from fortran import NumBAT
import reporting
import numbattools as nbtools
import materials


# def dec_float_str(dec_float):
#     ''' Convert float with decimal point into string with '_' in place of '.' '''
#     try:
#         s = '%8.1f' % dec_float
#     except Exception:
#         s = ''
#     # if type(dec_float) is float or type(dec_float) is int:
#     #    s = '%8.1f' % dec_float
#     # else:
#     #    s = ''
#     # s = s.replace('.', '_')
#     s = s.replace('.', ',')
#     s = s.replace(' ', '')
#     return s


def _load_waveguide_templates(p_msh_dir, p_msh_index):
    '''Loads and instantiates waveguide templates returning an index of waveguide inc_shape to class.'''

    try:
        with open(p_msh_index) as fin:
            msh_index = json.load(fin)['wguides']
    except Exception as ex:
        raise Exception('JSON parse error in reading user mesh index file' + str(ex)) from ex

    for msh in msh_index:
        if not msh.get('active', 1):  # this mesh is turned off for now
            continue

        mshnm = msh['inc_shape']
        mshpy = msh['wg_impl']
        mshclsnm = msh['wg_class']

        pth_mod = Path(p_msh_dir, mshpy)

        if not pth_mod.exists():
        #    raise Exception
            reporting.report_and_exit(f'Missing waveguide template implementation file: {str(mshpy)} named in str({p_msh_index}).'
                            + f'\nFile was expected in directory {str(p_msh_dir)}')

        spec = importlib.util.spec_from_file_location(mshnm, pth_mod)  # TODO: is mshnm the right name for this?
        try:
            py_mod = spec.loader.load_module()
        except Exception as ex:
            reporting.report_and_exit(f"Python couldn't load the user module '{str(mshpy)}'."+
                                      "\n Your code likely contains a syntax error or an attempt to load another module that was not successful." + f'\n\n The Python traceback contained the error message:\n   {str(ex)}.'+
                                      f'\n\n\n The full traceback was {traceback.format_exc()}')
        if not hasattr(py_mod, mshclsnm):
            reporting.report_and_exit(f"Can't find waveguide template class {mshclsnm} in implementation file: {mshpy}")
        mshcls = getattr(py_mod, mshclsnm)

        msh['mesh_cls'] = mshcls

    return msh_index


class Structure(object):
    ''' Represents the geometry and  material properties (elastic and electromagnetic) of a waveguide structure.

        Args:
            unitcell_x  (float):
                The horizontal period of the unit cell in nanometers.

            inc_a_x  (float):
                The horizontal diameter of the primary inclusion in nm.

        Keyword Args:
            unitcell_y  (float):
                The vertical period of the unit cell in nanometers. If None, unitcell_y = unitcell_x.

            inc_a_y  (float):
                The vertical diameter of the primary inclusion in nm.

            inc_shape  (str):
                Shape of inclusions that have template mesh, currently:
                 ``circular`` ``rectangular`` ``slot`` ``rib`` ``slot_coated`` ``rib_coated``
                 ``rib_double_coated`` ``pedestal`` ``onion`` ``onion2`` ``onion3``.
                 rectangular is default.

            slab_a_x  (float):
                The horizontal diameter in nm of the slab directly below the inclusion.

            slab_a_y  (float):
                The vertical diameter in nm of the slab directly below the inclusion.

            slab_b_x  (float):
                The horizontal diameter in nm of the slab separated from the inclusion by slab_a.

            slab_b_y  (float):
                The vertical diameter in nm of the slab separated from the inclusion by slab_a.

            two_inc_sep  (float):
                Separation between edges of two inclusions in nm.

            incs_y_offset  (float):
                 Vertical offset between centers of two inclusions in nm.

            coat_x  (float): The width of the first coat layer around
                the inclusion.

            coat_y  (float):
                The thickness of the first coating layer around the inclusion.

            coat2_x  (float):
                The width of the second coating layer around the inclusion.

            coat2_y  (float):
                The thickness of the second coating layer around the inclusion.

            symmetry_flag  (bool):
                True if materials all have sufficient symmetry that their tensors contain only 3 unique values.
                If False must specify full [3,3,3,3] tensors.

            material_bkg  (Material):
                The outer background material.

            material_a  (Material):
                The primary inclusion material.

            material_b-r (Material):
                Materials of additional layers.

            loss  (bool): If False, Im(n) = 0, if True n as in
                ``Material`` instance.

            make_mesh_now  (bool): If True, program creates a FEM mesh with
                provided :Struct: parameters. If False, must provide
                mesh_file name of existing .mail that will be run despite
                :Struct: parameters.

            force_mesh  (bool): If True, a new mesh is created despite
                existence of mesh with same parameter. This is used to make
                mesh with equal period etc. but different lc refinement.

            mesh_file  (str): If using a set pre-made mesh give its name
                including .mail. It must be located in backend/fortran/msh/
                Note: len(mesh_file) < 100.

            plt_mesh  (bool): Plot a png of the geometry and mesh files.

            check_mesh  (bool): Inspect the geometry and mesh files in gmsh.

            lc_bkg  (float): Length constant of meshing of background medium
                (smaller = finer mesh)

            lc_refine_1  (float): factor by which lc_bkg will be reduced on inclusion
                surfaces; lc_surface = lc_bkg / lc_refine_1. Larger lc_refine_1 = finer mesh.

            lc_refine_2-6'  (float): factor by which lc_bkg will be reduced on chosen
                surfaces; lc_surface = lc_bkg / lc_refine_2. see relevant .geo files.

            plotting_fields  (bool): Unless set to true field data deleted.
                Also plots modes (ie. FEM solutions) in gmsh format.
                Plots epsilon*|E|^2 & choice of real/imag/abs of
                x,y,z components & field vectors. Fields are saved as gmsh
                files, but can be converted by running the .geo file found in
                Bloch_fields/PNG/

            plot_real  (bool): Choose to plot real part of modal fields.

            plot_imag  (bool): Choose to plot imaginary part of modal fields.

            plot_abs  (bool): Choose to plot absolute value of modal fields.
    '''


    _user_mesh_templates = {}

    # called at startup from NumBATApp.__init__
    @classmethod
    def initialise_waveguide_templates(cls, nbapp):
        pmsh_dir = nbapp.path_mesh_templates()

        pmsh_index_builtin = Path(pmsh_dir, 'builtin_waveguides.json')
        pmsh_index_user = Path(pmsh_dir, 'user_waveguides.json')

        if not pmsh_index_builtin.exists():
            reporting.register_warning(
                f"Couldn't find builtin waveguide template index file: {pmsh_index_builtin}")
        else:
            cls._user_mesh_templates = _load_waveguide_templates(pmsh_dir, pmsh_index_builtin)
        
        if not pmsh_index_user.exists():
            reporting.register_warning(
                f"Couldn't find user waveguide template index file: {pmsh_index_user}")
        else:
            user_wg_templates =  _load_waveguide_templates(pmsh_dir, pmsh_index_user)
            cls._user_mesh_templates.extend(user_wg_templates)



    def __init__(self, unitcell_x, inc_a_x,
                 unitcell_y=None, inc_a_y=None, inc_shape='rectangular',
                 slab_a_x=None, slab_a_y=None, slab_b_x=None, slab_b_y=None,
                 coat_x=None, coat_y=None, coat2_x=None, coat2_y=None,
                 inc_b_x=None, inc_b_y=None,
                 two_inc_sep=None, incs_y_offset=None,
                 pillar_x=None, pillar_y=None,
                 inc_c_x=None, inc_d_x=None, inc_e_x=None, inc_f_x=None,
                 inc_g_x=None, inc_h_x=None, inc_i_x=None, inc_j_x=None,
                 inc_k_x=None, inc_l_x=None, inc_m_x=None, inc_n_x=None,
                 inc_o_x=None,
                 material_bkg=None,
                 material_a=None, material_b=None,
                 material_c=None, material_d=None,
                 material_e=None, material_f=None,
                 material_g=None, material_h=None,
                 material_i=None, material_j=None,
                 material_k=None, material_l=None,
                 material_m=None, material_n=None,
                 material_o=None, material_p=None,
                 material_q=None, material_r=None,
                 loss=True, symmetry_flag=True,
                 make_mesh_now=True, force_mesh=True,
                 mesh_file='NEED_FILE.mail', check_mesh=False, plt_mesh=False,
                 lc_bkg=0.09, lc_refine_1=1.0, lc_refine_2=1.0, lc_refine_3=1.0, lc_refine_4=1.0, lc_refine_5=1.0,
                 plotting_fields=False, plot_real=1, plot_imag=0, plot_abs=0, plot_field_conc=False,
                 direct_call=True):

        if plt_mesh:
            print("\n Warning: Option 'plt_mesh' is deprecated. Call method .plot_mesh() on your Struct object.\n\n")

        if check_mesh:
            print("\n Warning: Option 'check_mesh' is deprecated. Call method .check_mesh() on your Struct object.\n\n")

        numbat.assert_numbat_object_created()
        if direct_call:
            reporting.register_warning(
                'Calling objects.Structure directly is deprecated. Please switch to calling nbapp.make_structure()')

        self.shift_em_x = 0  # user requested offsets to coord-system
        self.shift_em_y = 0
        self.shift_ac_x = 0  # user requested offsets to coord-system
        self.shift_ac_y = 0


        # Structures material properties - need to check geometry definition
        # to ensure connecting material type with correct surface of geometry
        mat_vac = materials.make_material('Vacuum')

        self.d_materials = {}
        mat_pairs = {'bkg': material_bkg,
                     'a': material_a, 'b': material_b, 'c': material_c, 'd': material_d, 'e': material_e,
                     'f': material_f, 'g': material_g, 'h': material_h, 'i': material_i, 'j': material_j,
                     'k': material_k, 'l': material_l, 'm': material_m, 'n': material_n, 'o': material_o,
                     'p': material_p, 'q': material_q, 'r': material_r}
        for (tag, mat) in mat_pairs.items():
            self.d_materials[tag] = mat if (mat is not None) else copy.deepcopy(mat_vac)

        # Structures geometric shapes
        self.inc_shape = inc_shape
        unitcell_x = float(unitcell_x)
        self.inc_a_x = inc_a_x
        self.unitcell_x = unitcell_x

        if unitcell_y is None:
            unitcell_y = float(unitcell_x)
        else:
            unitcell_y = float(unitcell_y)

        if inc_a_y is None:
            inc_a_y = float(inc_a_x)
        else:
            inc_a_y = float(inc_a_y)

        if inc_b_x is not None:
            if inc_b_y is None:
                inc_b_y = float(inc_b_x)
            else:
                inc_b_y = float(inc_b_y)

        self.all_params = {'lc': lc_bkg}
        for p in ['unitcell_x', 'inc_a_x', 'unitcell_y', 'inc_a_y', 'inc_shape', 'slab_a_x',
                  'slab_a_y', 'slab_b_x', 'slab_b_y', 'coat_x', 'coat_y', 'coat2_x', 'coat2_y',
                  'inc_b_x', 'inc_b_y', 'two_inc_sep', 'incs_y_offset', 'pillar_x', 'pillar_y',
                  'inc_c_x', 'inc_d_x', 'inc_e_x', 'inc_f_x', 'inc_g_x', 'inc_h_x', 'inc_i_x',
                  'inc_j_x', 'inc_k_x', 'inc_l_x', 'inc_m_x', 'inc_n_x', 'inc_o_x',
                  'lc_refine_1', 'lc_refine_2', 'lc_refine_3', 'lc_refine_4', 'lc_refine_5']:
            self.all_params[p] = eval(p)

        self.loss = loss       # TODO: Boolean. Needs better name
        self.force_mesh = force_mesh

        if make_mesh_now:
            self.build_waveguide_geometry(self.d_materials)
        else:  # TODO: this seems to be broken. But also not really worth supporting? Mesh construction is not hard
            print(f"Using mesh from existing file '{mesh_file}'.")
            self.mesh_file = mesh_file
            with open(self.mesh_file) as f:
                self.mail_data = f.readlines()

        if plotting_fields:  # TODO: This is for internal fortran plotting. Should have a less appealing name
            reporting.register_warning('''Calling plotting_fields in objects.Structure.
                                       This is deprecated. Please report to the github page.''')
            self.plotting_fields = 1
            if not os.path.exists('Bloch_fields'):
                os.mkdir('Bloch_fields')
            if not os.path.exists('Bloch_fields/PDF'):
                os.mkdir('Bloch_fields/PDF')
            if not os.path.exists('AC_fields'):
                os.mkdir('AC_fields')
        else:
            self.plotting_fields = 0

        self.plot_real = plot_real
        self.plot_imag = plot_imag
        self.plot_abs = plot_abs
        self.plot_field_conc = plot_field_conc

        self.symmetry_flag = symmetry_flag

        # print('symflag', symmetry_flag)
        # el_conv_table = {}
        # i = 1; j = 1
        # for matter in acoustic_props:
        #     if matter != None:
        #         el_conv_table[i] = j
        #         j += 1
        #     i += 1
        # self.typ_el_AC = el_conv_table
        # print el_conv_table

        self._build_elastic_tensors(symmetry_flag)

    def _build_elastic_tensors(self, symmetry_flag):

        # construct list of materials with nonzero density, ie with acoustic properties likely defined
        # Any material not given acoustic_props assumed to be vacuum.
        acoustic_props = [
            m for m in self.d_materials.values() if m.has_elastic_properties()]

        # Number of different acoustic materials
        self.n_typ_el_AC = len(acoustic_props)

        rho = np.zeros(self.n_typ_el_AC)

        # stiffness tensor in 6x6 Voigt notation
        c_tensor = np.zeros((6, 6, self.n_typ_el_AC))

        # stiffness tensor as rank 4 ijkz tensor
        c_tensor_z = np.zeros((3, 3, 3, self.n_typ_el_AC))

        # photelastic tensor as rank 4 ijkl tensor
        p_tensor = np.zeros((3, 3, 3, 3, self.n_typ_el_AC))

        # eta tensor as rank 4 ijkl tensor
        eta_tensor = np.zeros((3, 3, 3, 3, self.n_typ_el_AC))

        voigt_map = {(0, 0): 1, (1, 1): 2, (2, 2): 3, (2, 1): 4,
                     (2, 0): 5, (0, 1): 6, (1, 2): 4, (0, 2): 5, (1, 0): 6}

        for k_typ in range(self.n_typ_el_AC):
            if acoustic_props[k_typ]:
                t_ac = acoustic_props[k_typ]
                t_ac_c = t_ac.c_tensor
                t_ac_p = t_ac.p_tensor
                t_ac_eta = t_ac.eta_tensor

                rho[k_typ] = t_ac.rho

                if symmetry_flag:  # is it actually worth making this saving?
                    c_tensor[0, 0, k_typ] = t_ac_c[1, 1]
                    c_tensor[1, 1, k_typ] = t_ac_c[1, 1]
                    c_tensor[2, 2, k_typ] = t_ac_c[1, 1]
                    c_tensor[0, 1, k_typ] = t_ac_c[1, 2]
                    c_tensor[0, 2, k_typ] = t_ac_c[1, 2]
                    c_tensor[1, 0, k_typ] = t_ac_c[1, 2]
                    c_tensor[1, 2, k_typ] = t_ac_c[1, 2]
                    c_tensor[2, 0, k_typ] = t_ac_c[1, 2]
                    c_tensor[2, 1, k_typ] = t_ac_c[1, 2]
                    c_tensor[3, 3, k_typ] = t_ac_c[4, 4]
                    c_tensor[4, 4, k_typ] = t_ac_c[4, 4]
                    c_tensor[5, 5, k_typ] = t_ac_c[4, 4]

                    c_tensor_z[2, 2, 2, k_typ] = t_ac_c[1, 1]
                    c_tensor_z[2, 0, 0, k_typ] = t_ac_c[1, 2]
                    c_tensor_z[2, 1, 1, k_typ] = t_ac_c[1, 2]
                    c_tensor_z[1, 1, 2, k_typ] = t_ac_c[4, 4]
                    c_tensor_z[1, 2, 1, k_typ] = t_ac_c[4, 4]
                    c_tensor_z[0, 0, 2, k_typ] = t_ac_c[4, 4]
                    c_tensor_z[0, 2, 0, k_typ] = t_ac_c[4, 4]

                    p_tensor[0, 0, 0, 0, k_typ] = t_ac_p[1, 1]
                    p_tensor[1, 1, 1, 1, k_typ] = t_ac_p[1, 1]
                    p_tensor[2, 2, 2, 2, k_typ] = t_ac_p[1, 1]
                    p_tensor[0, 0, 1, 1, k_typ] = t_ac_p[1, 2]
                    p_tensor[0, 0, 2, 2, k_typ] = t_ac_p[1, 2]
                    p_tensor[1, 1, 0, 0, k_typ] = t_ac_p[1, 2]
                    p_tensor[1, 1, 2, 2, k_typ] = t_ac_p[1, 2]
                    p_tensor[2, 2, 0, 0, k_typ] = t_ac_p[1, 2]
                    p_tensor[2, 2, 1, 1, k_typ] = t_ac_p[1, 2]
                    p_tensor[1, 2, 1, 2, k_typ] = t_ac_p[4, 4]
                    p_tensor[1, 2, 2, 1, k_typ] = t_ac_p[4, 4]
                    p_tensor[2, 1, 1, 2, k_typ] = t_ac_p[4, 4]
                    p_tensor[2, 1, 2, 1, k_typ] = t_ac_p[4, 4]
                    p_tensor[0, 2, 0, 2, k_typ] = t_ac_p[4, 4]
                    p_tensor[0, 2, 2, 0, k_typ] = t_ac_p[4, 4]
                    p_tensor[2, 0, 0, 2, k_typ] = t_ac_p[4, 4]
                    p_tensor[2, 0, 2, 0, k_typ] = t_ac_p[4, 4]
                    p_tensor[0, 1, 0, 1, k_typ] = t_ac_p[4, 4]
                    p_tensor[0, 1, 1, 0, k_typ] = t_ac_p[4, 4]
                    p_tensor[1, 0, 0, 1, k_typ] = t_ac_p[4, 4]
                    p_tensor[1, 0, 1, 0, k_typ] = t_ac_p[4, 4]

                    eta_tensor[0, 0, 0, 0, k_typ] = t_ac_eta[1, 1]
                    eta_tensor[1, 1, 1, 1, k_typ] = t_ac_eta[1, 1]
                    eta_tensor[2, 2, 2, 2, k_typ] = t_ac_eta[1, 1]
                    eta_tensor[0, 0, 1, 1, k_typ] = t_ac_eta[1, 2]
                    eta_tensor[0, 0, 2, 2, k_typ] = t_ac_eta[1, 2]
                    eta_tensor[1, 1, 0, 0, k_typ] = t_ac_eta[1, 2]
                    eta_tensor[1, 1, 2, 2, k_typ] = t_ac_eta[1, 2]
                    eta_tensor[2, 2, 0, 0, k_typ] = t_ac_eta[1, 2]
                    eta_tensor[2, 2, 1, 1, k_typ] = t_ac_eta[1, 2]
                    eta_tensor[1, 2, 1, 2, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[1, 2, 2, 1, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[2, 1, 1, 2, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[2, 1, 2, 1, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[0, 2, 0, 2, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[0, 2, 2, 0, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[2, 0, 0, 2, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[2, 0, 2, 0, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[0, 1, 0, 1, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[0, 1, 1, 0, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[1, 0, 0, 1, k_typ] = t_ac_eta[4, 4]
                    eta_tensor[1, 0, 1, 0, k_typ] = t_ac_eta[4, 4]

                else:
                    for i in range(6):
                        for j in range(6):
                            c_tensor[i, j, k_typ] = t_ac_c[i+1, j+1]

                    for i in [0, 1, 2]:
                        for j in [0, 1, 2]:
                            I = voigt_map[(i, j)]
                            for k in [0, 1, 2]:
                                Jz = voigt_map[(k, 2)]
                                c_tensor_z[i, j, k, k_typ] = t_ac_c[I, Jz]
                                for l in [0, 1, 2]:
                                    J = voigt_map[(k, l)]
                                    p_tensor[i, j, k, l, k_typ] = t_ac_p[I, J]
                                    eta_tensor[i, j, k, l,
                                               k_typ] = t_ac_eta[I, J]

        self.rho = rho
        self.c_tensor = c_tensor
        self.c_tensor_z = c_tensor_z
        self.p_tensor = p_tensor
        self.eta_tensor = eta_tensor


    def get_material(self, k):
        return self.d_materials[k]

    def set_xyshift_em(self, x, y):
        # Sets shift in grid from user perspective in nm
        self.shift_em_x = x*1e-9
        self.shift_em_y = y*1e-9

    def set_xyshift_ac(self, x, y):
        self.shift_ac_x = x*1e-9
        self.shift_ac_y = y*1e-9

    def _new_mesh_required(self):  # TODO: msh_name ?
        return self.force_mesh or not os.path.exists(self.msh_location_in + msh_name + '.mail')




    def get_mail_data(self):
        '''Returns list of lines in the NumBAT mesh format .mail file'''
        return self.mail_data

    def build_waveguide_geometry(self, d_materials):
        ''' Take the parameters specified in python and make a Gmsh FEM mesh.
            Creates a .geo and .msh file from the .geo template,
            then uses Fortran conv_gmsh routine
            to convert .msh into .mail, which is used in NumBAT FEM routine.
        '''

        print('Building mesh')

        # full path to backend directory that this code file is in
        this_directory = os.path.dirname(os.path.realpath(__file__))

        # msh_location_in = Path(this_directory, 'fortran', 'msh', '')  # msh directory inside backend
        # msh directory inside backend

        self.msh_location_in = Path(this_directory, 'msh')
        self.msh_location_out = Path(self.msh_location_in, 'build')

        self.wg_geom = None

        if not self.msh_location_out.is_dir():
            self.msh_location_out.mkdir()

        found = False 
        for msh in Structure._user_mesh_templates: 
            if self.inc_shape in msh['inc_shape']:
                found = True
                wg_geom = msh['mesh_cls'](self.all_params, d_materials)
                wg_geom.init_geometry()

                self.n_typ_el = wg_geom.num_type_elements()

                break

        if not found:
            raise NotImplementedError(f"\n Selected inc_shape = '{self.inc_shape}' "
                                      'is not currently implemented. \nPlease make a mesh with gmsh and '
                                      'consider contributing this to NumBAT via github.')

        self.linear_element_shapes = []
        self.curvilinear_element_shapes = []

        if wg_geom is not None:
            if wg_geom.is_curvilinear():
                self.linear_element_shapes.append(wg_geom.geom_name())
            else:
                self.curvilinear_element_shapes.append(wg_geom.geom_name())

        self.wg_geom = wg_geom

        self._build_mesh()


    def _build_mesh(self):
        '''Instantiates template gmsh file to specific gmsh then runs conv_gmsh to generate NumBAT .mail file.'''

        msh_template = self.wg_geom.gmsh_template_filename()
        msh_fname = self.wg_geom.get_instance_filename()

        if self._new_mesh_required():

            geo = self.wg_geom.make_geometry(numbat.NumBATApp().path_mesh_templates())

            fname = Path(self.msh_location_out, msh_fname)

            with open(str(fname) + '.geo', 'w') as fout:
                fout.write(geo)

            # Convert our Gmsh .geo file into Gmsh .msh and then NumBAT .mail
            assertions_on = False
            
            err_no, err_msg = NumBAT.conv_gmsh(str(fname), assertions_on)
            if err_no != 0:

                s = f'Terminating after Fortran error in processing .geo file "{fname}%s.geo".'
                if len(err_msg):
                    s += f'\nMessage was:\n {err_msg}'
                s += f'''

                Is the mesh template file "backend/fortran/msh/{msh_template}_msh_template.geo" designed correctly?
                To help diagnose the problem, try viewing the generated mesh file in gmsh by running:
                  gmsh {fname}.geo'''
                reporting.report_and_exit(s)

        self.mesh_file = str(fname) + '.mail'
        # TODO: curently used onyl to generate filenames for plot_mesh. Needed?
        self.msh_name = msh_fname

        # read in first line giving number of msh points and elements
        with open(self.mesh_file) as f:
            self.mail_data = f.readlines()



    def plot_mesh(self, outpref):
        '''Visualise mesh with gmsh and save to a file.'''

        # Manipulate scripts in backend/fortran/build
        # Writes final png file to user directory

        nbapp = numbat.NumBATApp()
        gmsh_exe = nbapp.path_gmsh()


        outprefix = Path(numbat.NumBATApp().outdir(), outpref)
        tdir = tempfile.TemporaryDirectory()
        tmpoutpref = str(Path(tdir.name, outpref))

        conv_tmp = open(Path(self.msh_location_in, 'geo2png.scr'), 'r').read()
        conv = conv_tmp.replace('tmp', str(tmpoutpref) + '-mesh_geom')

        fn_scr = Path(self.msh_location_out, self.msh_name + '_geo2png.scr')

        cmd = [gmsh_exe, self.msh_name + '.geo',
               self.msh_name + '_geo2png.scr']

        with open(fn_scr, 'w') as fout:
            fout.write(conv)
        subprocess.Popen(cmd, cwd=self.msh_location_out)

        os.wait()

        conv_tmp = open(Path(self.msh_location_in, 'msh2png.scr'), 'r').read()
        conv = conv_tmp.replace('tmp', str(tmpoutpref) + '-mesh_nodes')

        fn_scr = Path(self.msh_location_out, self.msh_name + '_msh2png.scr')

        cmd = [gmsh_exe, self.msh_name + '.msh',
               self.msh_name + '_msh2png.scr']
        with open(fn_scr, 'w') as fout:
            fout.write(conv)

        subprocess.Popen(cmd, cwd=self.msh_location_out)

        os.wait()

        nbtools.join_figs(tmpoutpref+'-mesh_geom.png',
                          tmpoutpref+'-mesh_nodes.png',
                          str(outprefix)+'-mesh.png')

    def check_mesh(self):
        '''Visualise geometry and mesh with gmsh.'''


        print('checking mesh')
        nbapp = numbat.NumBATApp()
        gmsh_exe = str(nbapp.path_gmsh())

        gmsh_cmd = gmsh_exe + " " + str(self.msh_location_out) + '/'+self.msh_name + '.geo'
        print(gmsh_cmd)
        os.system(gmsh_cmd)

        gmsh_cmd = gmsh_exe + " " + str(self.msh_location_out) + '/'+self.msh_name + '.msh'
        os.system(gmsh_cmd)




    def calc_EM_modes(self, num_modes, wl_nm, n_eff, Stokes=False, debug=False,
                      **args):
        ''' Run a simulation to find the Struct's EM modes.

            Args:
                num_modes  (int): Number of EM modes to solve for.

                wl_nm  (float): Wavelength of EM wave in vacuum in nanometres.

                n_eff  (float): Guesstimated effective index of
                    fundamental mode, will be origin of FEM search.

            Returns:
                ``Simulation`` object
        '''
        sim = Simulation(self, num_modes=num_modes, wl_nm=wl_nm,
                         n_eff=n_eff, Stokes=Stokes, debug=debug, **args)

        print('Calculating EM modes:')
        sim.calc_EM_modes()
        return sim

    def calc_AC_modes(self, num_modes, q_AC,
                      shift_Hz=None, EM_sim=None, bcs=None, debug=False, **args):
        ''' Run a simulation to find the Struct's acoustic modes.

            Args:
                num_modes  (int): Number of AC modes to solve for.

            Keyword Args:
                q_AC  (float): Wavevector of AC modes.

                shift_Hz  (float): Guesstimated frequency of modes,
                    will be origin of FEM search. NumBAT will make
                    an educated guess if shift_Hz=None.
                    (Technically the shift and invert parameter).

                EM_sim  (``Simulation`` object): Typically an acoustic
                    simulation follows on from an optical one.
                    Supply the EM ``Simulation`` object so the AC FEM mesh
                    can be constructed from this.
                    This is done by removing vacuum regions.

            Returns:
                ``Simulation`` object
        '''

        sim = Simulation(self, num_modes=num_modes, q_AC=q_AC,
                         shift_Hz=shift_Hz, EM_sim=EM_sim, debug=debug, **args)

        print('\n\nCalculating AC modes')

        sim.calc_AC_modes(bcs)
        return sim


# called at startup from NumBATApp.__init__
def initialise_waveguide_templates(numbatapp):
    Structure.initialise_waveguide_templates(numbatapp)

