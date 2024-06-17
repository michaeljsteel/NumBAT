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



#TODO: reduce number of imports

import os
import subprocess
import copy
import traceback

import tempfile
from pathlib import Path
import json
import importlib
import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate



import numbat
import reporting
from nbtypes import SI_nm

import materials
from fortran import nb_fortran
import plotting
import plottools
from mode_calcs import EMSimulation, ACSimulation
import nbgmsh

def _load_waveguide_templates(p_wgtemplate_dir, p_wgtemplate_index):
    '''Loads and instantiates waveguide templates specified in .json file {p_wgtemplate_index} in directory {p_wgtemplate_dir}.

    Returns an index of waveguide inc_shape to class.'''

    try:
        with open(p_wgtemplate_index) as fin:
            wg_index = json.load(fin)['wguides']
    except Exception as ex:
        raise Exception('JSON parse error in reading user mesh index file' + str(ex)) from ex

    for wg in wg_index:
        if not wg.get('active', 1):  # this mesh is turned off for now
            continue

        wgnm = wg['inc_shape']    # Name of inc_shape requested by user
        wgclsnm = wg['wg_class']  # Python class implementing this template
        wgpy = wg['wg_impl']      # Python file definind that class

        pth_mod = Path(p_wgtemplate_dir, wgpy)

        if not pth_mod.exists():
        #    raise Exception
            reporting.report_and_exit(f'Missing waveguide template implementation file: {str(wgpy)} named in str({p_wgtemplate_index}).'
                            + f'\nFile was expected in directory {str(p_wgtemplate_dir)}')

        spec = importlib.util.spec_from_file_location(wgnm, pth_mod)  # TODO: is wgnm the right name for this?
        try:
            py_mod = spec.loader.load_module()
        except Exception as ex:
            reporting.report_and_exit(f"Python couldn't load the user module '{str(wgpy)}'."+
                                      "\n Your code likely contains a syntax error or an attempt to load another module that was not successful." + f'\n\n The Python traceback contained the error message:\n   {str(ex)}.'+
                                      f'\n\n\n The full traceback was {traceback.format_exc()}')
        if not hasattr(py_mod, wgclsnm):
            reporting.report_and_exit(f"Can't find waveguide template class {wgclsnm} in implementation file: {wgpy}")
        wgcls = getattr(py_mod, wgclsnm)

        wg['wg_template_cls'] = wgcls

    return wg_index

class ElectromagneticProps:
    '''EM properties in unit-indexed forms suitable for fortran'''

    def __init__(self, v_mats_em, n_mats_em, loss):


        self.n_mats_em = n_mats_em         # number of materials with em properties (basically all used by structure)

        matvals = v_mats_em[:n_mats_em]    # Material objects of those

        # Set up mapping tables for refractive indices
        # (Why we need el_conv_table_n is mystery)
        # el_conv_table_n maps the number of the material to the position in the nonzero v_refindexn
        # el_conv_table_n[ith material] = index into v_refindexn  of non-zero refractive indices
        # Except for zero index materials,
        #  it will always be {1:1, 2:2, 3:3, .., num_mats:num_mats}

        # trivial identity map of index to active material
        self.el_conv_table_n = {i:i for i in range(1, self.n_mats_em+1)}   #_n for refractive index n


        # Array[0:n_mats_em] - refractive index of each active material

        self.v_refindexn =np.array([m.refindex_n for m in matvals])
        if not loss:
            self.v_refindexn = self.v_refindexn.real




class ElasticProps:
    '''Elastic tensors in unit-indexed forms suitable for fortran'''

    def __init__(self, v_acoustic_mats, symmetry_flag):

        self.n_mats_ac = len(v_acoustic_mats)
        self.v_active_mats = v_acoustic_mats

        # density  shape = [n_mats_ac]
        self.rho = np.zeros(self.n_mats_ac)

        # stiffness tensor in 6x6 Voigt notation [3 x 3  x  n_mats_ac]
        self.c_IJ = np.zeros((6, 6, self.n_mats_ac))

        # stiffness tensor as rank 4 ijkz tensor [3x3x3x1  x  n_mats_ac]
        self.c_ijkz = np.zeros((3, 3, 3, self.n_mats_ac))

        # photelastic tensor as rank 4 ijkl tensor  # [3x3x3x3  x  n_mats_ac]
        self.p_ijkl = np.zeros((3, 3, 3, 3, self.n_mats_ac))

        # eta tensor as rank 4 ijkl tensor [3x3x3x3  x  n_mats_ac]
        self.eta_ijkl = np.zeros((3, 3, 3, 3, self.n_mats_ac))

        self.fill_tensors(v_acoustic_mats, symmetry_flag)


    def extract_elastic_mats(self, structure, opt_props):
        el_conv_table = {}
        oldloc = 1
        newloc = 1
        d_mats_AC = {}

            #No need to examine any materials beyond the max in the EM simulation (they are all vacuum anyway)
        for mat in list(structure.d_materials.values())[:opt_props.n_mats_em]:
            if mat.has_elastic_properties():
                el_conv_table[oldloc] = newloc
                newloc += 1
                d_mats_AC[mat.material_name] = mat
            oldloc += 1

        self.typ_el_AC = {}
        for k, v in el_conv_table.items():
            # now keeps its own rather than take from simres_EM which might not exist
            self.typ_el_AC[opt_props.el_conv_table_n[k]] = v

        #TODO: are these two in any way different?
        print('building elastic lists', el_conv_table,  self.typ_el_AC)

    def is_elastic_material_index(self, idx):
        return idx in self.typ_el_AC

    def active_material_index(self, idx):
        return self.typ_el_AC[idx]

    def fill_tensors(self, v_acoustic_mats, symmetry_flag):

        # map a zero-indexed 3x3 elt to unit indexed 6x1 form.  eg x,x == 0,0 == 1
        # TODO: use a zero-indexed form of toVoigt map
        voigt_map = {(0, 0): 1, (1, 1): 2, (2, 2): 3, (2, 1): 4,
                     (2, 0): 5, (0, 1): 6, (1, 2): 4, (0, 2): 5, (1, 0): 6}


        # Build zero-based material tensors from unit-based
        for k_typ in range(self.n_mats_ac):
            if v_acoustic_mats[k_typ]:
                t_ac = v_acoustic_mats[k_typ]
                t_ac_c_IJ = t_ac.stiffness_c_IJ
                t_ac_p_IJ = t_ac.photoel_p_IJ
                t_ac_eta_IJ = t_ac.viscosity_eta_IJ

                self.rho[k_typ] = t_ac.rho

                if symmetry_flag:  # is it actually worth making this saving?
                    print('Surprise: using symmetry_flag tensor buildings.')
                    self.c_IJ[0, 0, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_IJ[1, 1, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_IJ[2, 2, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_IJ[0, 1, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[0, 2, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[1, 0, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[1, 2, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[2, 0, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[2, 1, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_IJ[3, 3, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_IJ[4, 4, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_IJ[5, 5, k_typ] = t_ac_c_IJ[4, 4]

                    self.c_ijkz[2, 2, 2, k_typ] = t_ac_c_IJ[1, 1]
                    self.c_ijkz[2, 0, 0, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_ijkz[2, 1, 1, k_typ] = t_ac_c_IJ[1, 2]
                    self.c_ijkz[1, 1, 2, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_ijkz[1, 2, 1, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_ijkz[0, 0, 2, k_typ] = t_ac_c_IJ[4, 4]
                    self.c_ijkz[0, 2, 0, k_typ] = t_ac_c_IJ[4, 4]

                    self.p_ijkl[0, 0, 0, 0, k_typ] = t_ac_p_IJ[1, 1]
                    self.p_ijkl[1, 1, 1, 1, k_typ] = t_ac_p_IJ[1, 1]
                    self.p_ijkl[2, 2, 2, 2, k_typ] = t_ac_p_IJ[1, 1]
                    self.p_ijkl[0, 0, 1, 1, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[0, 0, 2, 2, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[1, 1, 0, 0, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[1, 1, 2, 2, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[2, 2, 0, 0, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[2, 2, 1, 1, k_typ] = t_ac_p_IJ[1, 2]
                    self.p_ijkl[1, 2, 1, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[1, 2, 2, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 1, 1, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 1, 2, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 2, 0, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 2, 2, 0, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 0, 0, 2, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[2, 0, 2, 0, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 1, 0, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[0, 1, 1, 0, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[1, 0, 0, 1, k_typ] = t_ac_p_IJ[4, 4]
                    self.p_ijkl[1, 0, 1, 0, k_typ] = t_ac_p_IJ[4, 4]

                    self.eta_ijkl[0, 0, 0, 0, k_typ] = t_ac_eta_IJ[1, 1]
                    self.eta_ijkl[1, 1, 1, 1, k_typ] = t_ac_eta_IJ[1, 1]
                    self.eta_ijkl[2, 2, 2, 2, k_typ] = t_ac_eta_IJ[1, 1]
                    self.eta_ijkl[0, 0, 1, 1, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[0, 0, 2, 2, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[1, 1, 0, 0, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[1, 1, 2, 2, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[2, 2, 0, 0, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[2, 2, 1, 1, k_typ] = t_ac_eta_IJ[1, 2]
                    self.eta_ijkl[1, 2, 1, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[1, 2, 2, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 1, 1, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 1, 2, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 2, 0, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 2, 2, 0, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 0, 0, 2, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[2, 0, 2, 0, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 1, 0, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[0, 1, 1, 0, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[1, 0, 0, 1, k_typ] = t_ac_eta_IJ[4, 4]
                    self.eta_ijkl[1, 0, 1, 0, k_typ] = t_ac_eta_IJ[4, 4]

                else:
                    for i in range(6):
                        for j in range(6):
                            self.c_IJ[i, j, k_typ] = t_ac_c_IJ[i+1, j+1]  # TODO: replace with Voigt.value() ?

                    for i in [0, 1, 2]:
                        for j in [0, 1, 2]:
                            I = voigt_map[(i, j)]
                            for k in [0, 1, 2]:
                                Jz = voigt_map[(k, 2)]
                                self.c_ijkz[i, j, k, k_typ] = t_ac_c_IJ[I, Jz]
                                for l in [0, 1, 2]:
                                    J = voigt_map[(k, l)]
                                    self.p_ijkl[i, j, k, l, k_typ] = t_ac_p_IJ[I, J]
                                    self.eta_ijkl[i, j, k, l, k_typ] = t_ac_eta_IJ[I, J]


class Structure:
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


    _waveguide_templates = {}   # class-level dictionary of waveguide template name to implementing class

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
            cls._waveguide_templates = _load_waveguide_templates(pmsh_dir, pmsh_index_builtin)

        if not pmsh_index_user.exists():
            reporting.register_warning(
                f"Couldn't find user waveguide template index file: {pmsh_index_user}")
        else:
            user_wg_templates =  _load_waveguide_templates(pmsh_dir, pmsh_index_user)
            cls._waveguide_templates.extend(user_wg_templates)



    def __init__(self, inc_shape=None, unitcell_x=None, unitcell_y=None,
                 inc_a_x=None, inc_a_y=None, inc_b_x=None, inc_b_y=None,
                 inc_c_x=None, inc_d_x=None, inc_e_x=None, inc_f_x=None,
                 inc_g_x=None, inc_h_x=None, inc_i_x=None, inc_j_x=None,
                 inc_k_x=None, inc_l_x=None, inc_m_x=None, inc_n_x=None,
                 inc_o_x=None,
                 slab_a_x=None, slab_a_y=None, slab_b_x=None, slab_b_y=None,
                 coat_x=None, coat_y=None, coat2_x=None, coat2_y=None,
                 two_inc_sep=None, incs_y_offset=None, pillar_x=None, pillar_y=None,
                 material_bkg=None,
                 material_a=None, material_b=None, material_c=None, material_d=None,
                 material_e=None, material_f=None, material_g=None, material_h=None,
                 material_i=None, material_j=None, material_k=None, material_l=None,
                 material_m=None, material_n=None, material_o=None, material_p=None,
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

        #self.n_mats_em = 0     # total number of materials _declared_ to be used in the structure. May not always be accurage

        # material identities seem to depend on stable ordering of the material dictionary when converted to a list

        self.optical_props = None
        self.elastic_props = None

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

        # Fill up dictionary with all mats 'bkg' to 'r' with vacuum if needed.
        # TODO: Why do this. Better to leave undefined and fail if needed materials are missing.

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

        n_mats_em = 0
        if make_mesh_now:
            n_mats_em = self.build_waveguide_geometry(self.d_materials)
        else:  # TODO: this seems to be broken. But also not really worth supporting? Mesh construction is not hard
            print(f"Using mesh from existing file '{mesh_file}'.")
            self.mesh_mail_fname = mesh_file
            with open(self.mesh_mail_fname) as f:
                self.mail_data = f.readlines()

        # TODO: much of this plotting stuffis no longer needed
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
        # for matter in v_acoustic_mats:
        #     if matter != None:
        #         el_conv_table[i] = j
        #         j += 1
        #     i += 1
        # self.typ_el_AC = el_conv_table
        # print el_conv_table

        self.optical_props = ElectromagneticProps(list(self.d_materials.values()), n_mats_em, self.loss)

        self._build_elastic_tensors(symmetry_flag)


    def _build_elastic_tensors(self, symmetry_flag):
        '''Put all elastic tensor properties into single 3 x 3 x num_material arrays, zero-indexed'''

        # construct list of materials with nonzero density, ie with acoustic properties likely defined
        # Any material not given v_acoustic_mats assumed to be vacuum.

        v_acoustic_mats = [m for m in self.d_materials.values() if m.has_elastic_properties()]

        self.elastic_props = ElasticProps(v_acoustic_mats, symmetry_flag)

    def get_material(self, k):
        return self.d_materials[k]

    def set_xyshift_em(self, x, y):
        # Sets shift in grid from user perspective in nm
        self.shift_em_x = x * SI_nm
        self.shift_em_y = y * SI_nm

    def set_xyshift_ac(self, x, y):
        self.shift_ac_x = x * SI_nm
        self.shift_ac_y = y * SI_nm

    def _new_mesh_required(self):  # TODO: msh_name ? REMOVE ME
        msh_name='missing_mesh_name'
        return self.force_mesh or not os.path.exists(self.msh_location_in + msh_name + '.mail')




    def get_mail_mesh_data(self):
        '''Returns Object representing mesh in  .mail file'''
        return self.mail_data

    def build_waveguide_geometry(self, d_materials):
        ''' Take the parameters specified in python and make a Gmsh FEM mesh.
            Creates a .geo and .msh file from the .geo template,
            then uses Fortran conv_gmsh routine
            to convert .msh into .mail, which is used in NumBAT FEM routine.
        '''

        print('Building waveguide')

        # full path to backend directory that this code file is in
        this_directory = os.path.dirname(os.path.realpath(__file__))

        # locations of gmsh input and output files
        self.msh_location_in = Path(this_directory, 'msh')
        self.msh_location_out = Path(self.msh_location_in, 'build')


        if not self.msh_location_out.is_dir():
            self.msh_location_out.mkdir()

        self.wg_geom = None

        found = False
        for wg in Structure._waveguide_templates:
            if self.inc_shape in wg['inc_shape']:  # is the desired shape supported by this template class?
                found = True

                # Instantiate the class that defines this waveguide model
                wg_geom = wg['wg_template_cls'](self.all_params, d_materials)

                wg_geom.init_geometry()

                n_mats_em = wg_geom.num_type_materials()  # This is number of distinct materials == element types materials declared by the template

                break

        if not found:
            raise NotImplementedError(f"\n Selected inc_shape = '{self.inc_shape}' "
                                      'is not currently implemented. \nPlease make a mesh with gmsh and '
                                      'consider contributing this to NumBAT via github.')

        self.linear_element_shapes = []
        self.curvilinear_element_shapes = []

        if wg_geom is not None:
            if wg_geom.is_curvilinear():  #TODO backwards!
                self.linear_element_shapes.append(wg_geom.geom_name())
            else:
                self.curvilinear_element_shapes.append(wg_geom.geom_name())

        self.wg_geom = wg_geom

        self._build_mesh()

        return n_mats_em

    def using_linear_elements(self):
        return self.inc_shape in self.linear_element_shapes

    def using_curvilinear_elements(self):
        return self.inc_shape in self.curvilinear_element_shapes



    def _build_mesh(self):
        '''Instantiates generic template gmsh file to aspecific gmsh then runs conv_gmsh
        to generate the NumBAT .mail file.'''

        print('  Processing Gmsh')

        msh_template = self.wg_geom.gmsh_template_filename()
        msh_fname = self.wg_geom.get_instance_filename()

        if self._new_mesh_required():

            # create string in .geo format from the template file with all parameters adjusted for our design
            geo = self.wg_geom.make_geometry(numbat.NumBATApp().path_mesh_templates())

            fname = Path(self.msh_location_out, msh_fname)

            with open(str(fname) + '.geo', 'w') as fout:
                fout.write(geo)

            # Convert our Gmsh .geo file into Gmsh .msh and then NumBAT .mail
            assertions_on = False

            err_no, err_msg = nb_fortran.conv_gmsh(str(fname), assertions_on)
            if err_no != 0:

                s = f'Terminating after Fortran error in processing .geo file "{fname}%s.geo".'
                if len(err_msg):
                    s += f'\nMessage was:\n {err_msg}'
                s += f'''

                Is the mesh template file "backend/fortran/msh/{msh_template}_msh_template.geo" designed correctly?
                To help diagnose the problem, try viewing the generated mesh file in gmsh by running:
                  gmsh {fname}.geo'''
                reporting.report_and_exit(s)

        self.mesh_mail_fname = str(fname) + '.mail'
        # TODO: curently used onyl to generate filenames for plot_mesh. Needed? Fix the filenames.
        self.msh_name = msh_fname

        # read in first line giving number of msh points and elements
        with open(self.mesh_mail_fname) as f:
            mail_lines= f.readlines()

        self.mail_data = nbgmsh.MailData(mail_lines)  #keep track of the Mail format


    def plot_mesh(self, outpref):
        '''Visualise mesh with gmsh and save to a file.'''

        # Manipulate scripts in backend/fortran/build
        # Writes final png file to user directory

        nbapp = numbat.NumBATApp()
        gmsh_exe = nbapp.path_gmsh()


        outprefix = Path(numbat.NumBATApp().outdir(), outpref)
        tdir = tempfile.TemporaryDirectory()
        tmpoutpref = str(Path(tdir.name, outpref))

        with open(Path(self.msh_location_in, 'geo2png.scr'), 'r') as fin:
            conv_tmp = fin.read()
        conv = conv_tmp.replace('tmp', str(tmpoutpref) + '-mesh_geom')

        fn_scr = Path(self.msh_location_out, self.msh_name + '_geo2png.scr')

        cmd = [gmsh_exe, self.msh_name + '.geo',
               self.msh_name + '_geo2png.scr']

        with open(fn_scr, 'w') as fout:
            fout.write(conv)
        subprocess.Popen(cmd, cwd=self.msh_location_out)

        os.wait()

        with open(Path(self.msh_location_in, 'msh2png.scr'), 'r') as fin:
            conv_tmp = fin.read()
        conv = conv_tmp.replace('tmp', str(tmpoutpref) + '-mesh_nodes')

        fn_scr = Path(self.msh_location_out, self.msh_name + '_msh2png.scr')

        cmd = [gmsh_exe, self.msh_name + '.msh',
               self.msh_name + '_msh2png.scr']
        with open(fn_scr, 'w') as fout:
            fout.write(conv)

        subprocess.Popen(cmd, cwd=self.msh_location_out)

        os.wait()

        plottools.join_figs([tmpoutpref+'-mesh_geom.png',
                          tmpoutpref+'-mesh_nodes.png',],
                          str(outprefix)+'-mesh.png',
                          clip=(60,50,60,50))

    def check_mesh(self):
        '''Visualise geometry and mesh with gmsh.'''

        nbapp = numbat.NumBATApp()
        gmsh_exe = str(nbapp.path_gmsh())

        gmsh_cmd = gmsh_exe + " " + str(self.msh_location_out) + '/'+self.msh_name + '.geo'
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
        sim = EMSimulation(self, num_modes=num_modes, wl_nm=wl_nm,
                         n_eff_target=n_eff, Stokes=Stokes, debug=debug, **args)

        sim.calc_modes()

        return sim.get_sim_result()

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

                EM_sim  (``EMSimResult`` object): Typically an acoustic
                    simulation follows on from an optical one.
                    Supply the ``EMSimResult`` object so the AC FEM mesh
                    can be constructed from this.
                    This is done by removing vacuum regions.

            Returns:
                ``Simulation`` object
        '''

        sim = ACSimulation(self, num_modes=num_modes, q_AC=q_AC,
                         shift_Hz=shift_Hz, simres_EM=EM_sim, debug=debug, **args)

        sim.calc_modes(bcs)

        return sim.get_sim_result()

    def plot_refractive_index_profile(self, prefix, as_epsilon=False):
        print('\n\nPlotting ref index')


        mail = self.mail_data
        v_x, v_y = mail.v_centx, mail.v_centy
        v_elt_indices = mail.v_elts[:,-1]  # the elt number column
        v_refindex = 0*v_x

        print('Mesh props', mail.n_msh_pts, mail.n_msh_elts)

        uniq_elts = set(list(v_elt_indices))


        # find ref index at each centroid
        v_elt_refindex = np.zeros(len(uniq_elts))

        for i in range(len(v_elt_refindex)):
            v_elt_refindex[i] = np.real(list(self.d_materials.values())[i].refindex_n)

        for i,elt in enumerate(v_elt_indices):
            v_refindex[i] = v_elt_refindex[elt-1]


        # Now have an irregular x,y,n array to interpolate onto.

        x_min = np.min(v_x)
        x_max = np.max(v_x)
        y_min = np.min(v_y)
        y_max = np.max(v_y)

        n_points = 40
        area = abs((x_max-x_min)*(y_max-y_min))
        n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
        n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))


        v_regx = np.linspace(x_min, x_max, n_pts_x)
        v_regy = np.linspace(y_min, y_max, n_pts_y)
        m_regx, m_regy = np.meshgrid(v_regx, v_regy)

        xy_in = np.array([v_x, v_y]).T
        xy_out = np.vstack([m_regx.ravel(), m_regy.ravel()]).T


        v_regindex = scipy.interpolate.griddata(xy_in, v_refindex, xy_out).reshape([n_pts_y, n_pts_x])
        fig, ax = plt.subplots()

        if as_epsilon:
            v_regindex = v_regindex**2
            fig.suptitle('Dielectric constant')
        else:
            fig.suptitle('Refractive index')

        cmap='cool'
        cf=ax.contourf(m_regx, m_regy, v_regindex, cmap=cmap, vmin=1.0)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        fig.colorbar(cf)


        plotting.save_and_close_figure(fig, prefix+'refn.png')


    def plot_phase_velocity_z_profile(self, prefix):
        print('plotting phasevel index')


# called at startup from NumBATApp.__init__
def initialise_waveguide_templates(numbatapp):
    Structure.initialise_waveguide_templates(numbatapp)

