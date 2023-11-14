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


#import time
import os
import subprocess
import copy


#import numbers

import uuid
import tempfile
import shutil
import platform
from pathlib import Path

import numpy as np


import matplotlib.patches as mplpatches
#import matplotlib.collections as mplcollections
#import matplotlib.colors as mplcolors


from mode_calcs import Simulation
import numbat
from fortran import NumBAT
import reporting
import numbattools as nbtools
import materials


g_onion_layer_thicknesses = ()


def is_real_number(x):
    # return isinstance(x, float) or isinstance(x, int)
    # return isinstance(x, numbers.Number)
    try:
        xx = float(x)
        return True
    except Exception:
        return False


def dec_float_str(dec_float):
    ''' Convert float with decimal point into string with '_' in place of '.' '''
    try:
        s = '%8.1f' % dec_float
    except Exception:
        s = ''
    # if type(dec_float) is float or type(dec_float) is int:
    #    s = '%8.1f' % dec_float
    # else:
    #    s = ''
    # s = s.replace('.', '_')
    s = s.replace('.', ',')
    s = s.replace(' ', '')
    return s


def _frm_onion(ax, layers):
    nmtoum = 1.e-3   # radii are in nm but plots are in microns

    r = 0
    for l in layers:
        if l is not None:
            r += l
            circ = mplpatches.Circle((0, 0), r*nmtoum, facecolor=None, fill=False, edgecolor='gray',
                                     linewidth=.75)
            ax.add_patch(circ)


def _my_frm_onion(ax):
    global g_onion_layer_thicknesses
    return _frm_onion(ax, g_onion_layer_thicknesses[:self.n_typ_el-1])


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
                 ``circular`` ``rectangular`` ``slot`` ``rib`` ``slot_coated`` ``rib_coated`` ``rib_double_coated`` ``pedestal`` ``onion`` ``onion2`` ``onion3``.
                 ectangular is default.

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
                 plotting_fields=False, plot_real=1, plot_imag=0, plot_abs=0,
                 plot_field_conc=False):

        numbat.assert_numbat_object_created()

        self.shift_em_x = 0  # user requested offsets to coord-system
        self.shift_em_y = 0
        self.shift_ac_x = 0  # user requested offsets to coord-system
        self.shift_ac_y = 0

        self.mpl_wg_frame_drawer = None

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
            self.d_materials[tag] = mat if (
                mat is not None) else copy.deepcopy(mat_vac)

        # Structures geometric shapes
        self.unitcell_x = float(unitcell_x)
        self.inc_a_x = inc_a_x

        if unitcell_y is None:
            self.unitcell_y = float(unitcell_x)
        else:
            self.unitcell_y = float(unitcell_y)

        if inc_a_y is None:
            self.inc_a_y = float(inc_a_x)
        else:
            self.inc_a_y = float(inc_a_y)

        self.inc_b_x = inc_b_x
        if inc_b_x is not None:
            if inc_b_y is None:
                self.inc_b_y = float(inc_b_x)
            else:
                self.inc_b_y = float(inc_b_y)

        # structure shape and dimensions
        self.inc_shape = inc_shape
        self.slab_a_x = slab_a_x
        self.slab_a_y = slab_a_y
        self.slab_b_x = slab_b_x
        self.slab_b_y = slab_b_y
        self.pillar_x = pillar_x
        self.pillar_y = pillar_y
        self.inc_c_x = inc_c_x
        self.inc_d_x = inc_d_x
        self.inc_e_x = inc_e_x
        self.inc_f_x = inc_f_x
        self.inc_g_x = inc_g_x
        self.inc_h_x = inc_h_x
        self.inc_i_x = inc_i_x
        self.inc_j_x = inc_j_x
        self.inc_k_x = inc_k_x
        self.inc_l_x = inc_l_x
        self.inc_m_x = inc_m_x
        self.inc_n_x = inc_n_x
        self.inc_o_x = inc_o_x
        self.coat_x = coat_x
        self.coat_y = coat_y
        self.coat2_x = coat2_x
        self.coat2_y = coat2_y
        self.two_inc_sep = two_inc_sep
        self.incs_y_offset = incs_y_offset

        self.loss = loss       # TODO: Boolean. Needs better name

        if plt_mesh:
            print("\n Warning: Option 'plt_mesh' is deprecated. Call method .plot_mesh() on your Struct object.\n\n")

        if check_mesh:
            print("\n Warning: Option 'check_mesh' is deprecated. Call method .check_mesh() on your Struct object.\n\n")

        # grid resolution controls
        self.lc = lc_bkg
        self.lc_refine_1 = lc_refine_1
        self.lc_refine_2 = lc_refine_2
        self.lc_refine_3 = lc_refine_3
        self.lc_refine_4 = lc_refine_4
        self.lc_refine_5 = lc_refine_5
        self.force_mesh = force_mesh

        if make_mesh_now:
            self.make_mesh(self.d_materials)
        else:  # TODO: this seems to be broken. But also not really worth supporting? Mesh construction is not hard
            print(f"Using mesh from existing file '{mesh_file}'.")
            self.mesh_file = mesh_file
            # read in first line giving number of msh points and elements
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



        self.rho, self.c_tensor, self.c_tensor_z, self.p_tensor, self.eta_tensor = self._build_elastic_tensors(symmetry_flag)

        self.linear_element_shapes = ['rectangular', 'slot', 'slot_coated', 'rib',
                                      'rib_coated', 'rib_double_coated', 'pedestal']
        self.curvilinear_element_shapes = ['circular', 'onion', 'onion2', 'onion3',
                                           'circ_onion', 'circ_onion2', 'circ_onion3']

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
                                    eta_tensor[i, j, k, l, k_typ] = t_ac_eta[I, J]

        return  rho,  c_tensor,  c_tensor_z,  p_tensor,  eta_tensor

#def set_threadsafe(self):
#        pass

    def get_material(self, k):
        return self.d_materials[k]

    def set_xyshift_em(self, x, y):
        # Sets shift in grid from user perspective in nm
        self.shift_em_x = x*1e-9
        self.shift_em_y = y*1e-9

    def set_xyshift_ac(self, x, y):
        self.shift_ac_x = x*1e-9
        self.shift_ac_y = y*1e-9

    def _build_new_mesh(self):  #TODO: msh_name ?
        return self.force_mesh or not os.path.exists(self.msh_location_in + msh_name + '.mail')

    def _load_mesh_template(self, msh_template):
        geo = open(self.msh_location_in + f'{msh_template}_msh_template.geo' , 'r').read()
        return geo

    # TODO: in a plug-in system for .geo templates, this function doesn't need to be in the class
    #      turn the _make_geo_* functions into non-member functions as call-backs

    def _make_mesh_name(self, msh_template, l_dims):
        '''Make name for the concrete instantiation of a given mesh .geo template'''
        msh_name = msh_template

        # made crazy long names, not helping
        # for v in l_dims:
        #    if is_real_number(v): msh_name += '_%s' % dec_float_str(v)

        msh_name += '--pid-%d' % os.getpid()

        # need to make name unique to support parallel processing
        msh_name += '--'+str(uuid.uuid4())

        return msh_name

    def _apply_geo_subs(self, msh_template, subs):
        '''Turn the template into a concrete instantiation of a given mesh .geo template'''
        geo = self._load_mesh_template(msh_template)
        for (olds, news, val) in subs:
            if val is None:  # unset value not overridden or dropped
                continue
            elif is_real_number(val):
                geo = geo.replace(olds, news % val)
            else:
                geo = geo.replace(olds, news)
        return geo

        # TODO: the 1_on_slab, 2_on_slab, 1_on_2slab etc seem never to have been hooked up.
        #      Are they useful?

    def _make_geo_circular_rectangular(self, inc_shape, d_materials):
        subs = ''

        if self.slab_b_x is not None:
            raise ValueError(
                f"NumBAT doesn't understand your geometry: with shape {self.inc_shape}, I did not expect values for slab_b.")
        elif self.slab_a_x is not None:
            raise ValueError(
                f"NumBAT doesn't understand your geometry: with shape {self.inc_shape}, I did not expect values for slab_a.")
        elif self.inc_a_x is not None:
            if self.coat_y is None and self.inc_b_x is None:  # One inclusion, no coating
                msh_template = 'oneincl'  # used to be just '1'
                self.n_typ_el = 2         # bkg, core (mat_a)

                msh_name = self._make_mesh_name(msh_template,
                                                (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y))

            elif self.coat_y is None and self.inc_b_x is not None:  # Two inclusions, no coating
                msh_template = 'twoincl'  # used to be just '2'
                self.n_typ_el = 3         # bkg, core 1 (mat_a), core 2 (mat_b)

                msh_name = self._make_mesh_name(msh_template,
                                                (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y, self.inc_a_x, self.inc_a_y))

            # Two inclusions, with coating # TODO:implement
            elif self.coat_y is not None and self.inc_b_x is not None:
                raise NotImplementedError(
                    'Have not implemented 2 coated inclusions.')

            elif self.coat_y is not None and self.inc_b_x is None:  # One inclusion, with coating # TODO:implement
                raise NotImplementedError(
                    'Have not implemented 1 coated inclusions.')

            else:
                raise ValueError("NumBAT doesn't understand your geometry.")
        else:
            raise ValueError('must have at least one nonzero inclusion.')

        if inc_shape == 'circular':
            msh_name += '-c'

        # TODO: these are crazy small defaults
        subs = [('d_in_nm = 100;', 'd_in_nm = %f;', self.unitcell_x)]
        subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.unitcell_y))
        subs.append(('a1 = 20;', 'a1 = %f;', self.inc_a_x))
        subs.append(('a1y = 10;', 'a1y = %f;', self.inc_a_y))

        if self.inc_shape == 'circular':
            subs.append(('rect = 1;', 'rect = 0;', ''))

        subs.append(('lc = 0;', 'lc = %f;', self.lc))
        subs.append(
            ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
        subs.append(
            ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))

        if msh_template in ['twoincl', '2', '2_on_s', '2_on_2s']:
            subs.append(('a2 = 10;', 'a2 = %f;', self.inc_b_x))
            subs.append(('a2y = 20;', 'a2y = %f;', self.inc_b_y))
            subs.append(('sep = 10;', 'sep = %f;', self.two_inc_sep))

            # geo = geo.replace('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;' % self.lc_refine_3)
        if msh_template == '2':
            subs.append(('yoff = -5;', 'yoff = %f;', self.incs_y_offset))

        if msh_template in ['1_on_slab', '1_on_2slabs', '1_on_slab', '2_on_2slabs']:
            subs.append(('slab_width = d_in_nm;',
                        'slab_width = %f;', self.slab_a_x))
            subs.append(
                ('slab_height = 10;', 'slab_height = %f;', self.slab_a_y))
            subs.append(
                ('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;', self.lc_refine_3))
            subs.append(
                ('lc_refine_4 = lc/1;', 'lc_refine_4 = lc/%f;', self.lc_refine_4))

        if msh_template in ['1_on_2slabs', '2_on_2slabs']:
            subs.append(('slab2_width = d_in_nm;',
                        'slab2_width = %f;', self.slab_b_x))
            subs.append(
                ('slab2_height = 5;', 'slab2_height = %f;', self.slab_b_y))
            subs.append(
                ('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;', self.lc_refine_3))
            subs.append(
                ('lc_refine_4 = lc/1;', 'lc_refine_4 = lc/%f;', self.lc_refine_4))

        return msh_template, msh_name, subs

    def _make_geo_slot(self, inc_shape, d_materials):
        subs = ''

        if self.coat_y is not None:
            msh_template = 'slot_coated'
            self.n_typ_el = 5

            msh_name = self._make_mesh_name(msh_template,
                                            (self.unitcell_x, self.unitcell_y, self.inc_a_x,
                                             self.inc_a_y, self.inc_b_x, self.slab_a_x,
                                             self.slab_a_y, self.coat_y))

            subs = [('d_in_nm = 100;', 'd_in_nm = %f;', self.unitcell_x)]
            subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.unitcell_y))
            subs.append(('a1 = 20;', 'a1 = %f;', self.inc_a_x))
            subs.append(('a1y = 10;', 'a1y = %f;', self.inc_a_y))
            subs.append(('a2 = 20;', 'a2 = %f;', self.inc_b_x))
            subs.append(('slabx = 80;', 'slabx = %f;', self.slab_a_x))
            subs.append(('slaby = 10;', 'slaby = %f;', self.slab_a_y))
            subs.append(('c1y = 10;', 'c1y = %f;', self.coat_y))
            subs.append(('lc = 0;', 'lc = %f;', self.lc))
            subs.append(
                ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
            subs.append(
                ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))

            # msh_name = 'slot_c_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(cc)s_%(ccc)s' % {
            msh_name = self._make_mesh_name(msh_template,
                                            (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y, self.inc_b_x, self.slab_a_y, self.coat_y))
        else:
            msh_template = 'slot'
            self.n_typ_el = 4
            # msh_name = 'slot_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(cc)s_%(ccc)s' % {
            msh_name = self._make_mesh_name(msh_template,
                                            (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y, self.inc_b_x, self.slab_a_x, self.slab_a_y))

            subs = [('d_in_nm = 100;', 'd_in_nm = %f;', self.unitcell_x)]
            subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.unitcell_y))
            subs.append(('a1 = 20;', 'a1 = %f;', self.inc_a_x))
            subs.append(('a1y = 10;', 'a1y = %f;', self.inc_a_y))
            subs.append(('a2 = 20;', 'a2 = %f;', self.inc_b_x))
            subs.append(('slabx = 80;', 'slabx = %f;', self.slab_a_x))
            subs.append(('slaby = 10;', 'slaby = %f;', self.slab_a_y))
            subs.append(('lc = 0;', 'lc = %f;', self.lc))
            subs.append(
                ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
            subs.append(
                ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))
            subs.append(
                ('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;', self.lc_refine_3))

        return msh_template, msh_name, subs

    def _make_geo_rib(self, inc_shape, d_materials):
        subs = ''
        msh_template = 'rib'

        # bkg, rib (mat_a), slab (mat_b), [outer ribs (mat_c)]
        if d_materials['c'].is_vacuum():  # TODO: perhaps a better test is whether bkg = mat_c
            self.n_typ_el = 3
        else:
            self.n_typ_el = 4

        msh_name = self._make_mesh_name(msh_template,
                                        (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y, self.slab_a_x, self.slab_a_y))

        subs = [('d_in_nm = 100;', 'd_in_nm = %f;', self.unitcell_x)]
        subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.unitcell_y))
        subs.append(('a1 = 20;', 'a1 = %f;', self.inc_a_x))
        subs.append(('a1y = 10;', 'a1y = %f;', self.inc_a_y))
        subs.append(('slabx = 80;', 'slabx = %f;', self.slab_a_x))
        subs.append(('slaby = 10;', 'slaby = %f;', self.slab_a_y))
        subs.append(('lc = 0;', 'lc = %f;', self.lc))
        subs.append(
            ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
        subs.append(
            ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))

        return msh_template, msh_name, subs

    def _make_geo_rib_coated(self, inc_shape, d_materials):
        subs = ''

        msh_template = 'rib_coated'
        self.n_typ_el = 4  # bkg, rib (mat_a), slab (mat_b), coating (mat_c)
        msh_name = self._make_mesh_name(msh_template,
                                        (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y, self.coat_x, self.coat_y, self.slab_a_x, self.slab_a_y))

        subs = [('d_in_nm = 100;', 'd_in_nm = %f;', self.unitcell_x)]
        subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.unitcell_y))
        subs.append(('a1 = 20;', 'a1 = %f;', self.inc_a_x))
        subs.append(('a1y = 10;', 'a1y = %f;', self.inc_a_y))
        subs.append(('slabx = 80;', 'slabx = %f;', self.slab_a_x))
        subs.append(('slaby = 10;', 'slaby = %f;', self.slab_a_y))
        subs.append(('coatx = 2;', 'coatx = %f;', self.coat_x))
        subs.append(('coaty = 2;', 'coaty = %f;', self.coat_y))
        subs.append(('lc = 0;', 'lc = %f;', self.lc))
        subs.append(
            ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
        subs.append(
            ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))
        subs.append(
            ('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;', self.lc_refine_3))

        return msh_template, msh_name, subs

    def _make_geo_rib_double_coated(self, inc_shape, d_materials):
        subs = ''
        msh_template = 'rib_double_coated'
        # bkg, rib (mat_a), slab (mat_b), top inner coating (mat_c), bottom coating (mat_d), top outer-coating (mat_e)
        self.n_typ_el = 6

        msh_name = self._make_mesh_name(msh_template,
                                        (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y,
                                         self.coat_x, self.coat_y, self.coat2_y, self.slab_a_x, self.slab_a_y, self.slab_b_y))

        subs = [('d_in_nm = 100;', 'd_in_nm = %f;', self.unitcell_x)]
        subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.unitcell_y))
        subs.append(('a1 = 20;', 'a1 = %f;', self.inc_a_x))
        subs.append(('a1y = 10;', 'a1y = %f;', self.inc_a_y))
        subs.append(('slabx = 80;', 'slabx = %f;', self.slab_a_x))
        subs.append(('slaby = 10;', 'slaby = %f;', self.slab_a_y))
        subs.append(('slab2y = 5;', 'slab2y = %f;', self.slab_b_y))
        subs.append(('coatx = 2;', 'coatx = %f;', self.coat_x))
        subs.append(('coaty = 2;', 'coaty = %f;', self.coat_y))
        subs.append(('coat2x = 4;', 'coat2x = %f;', self.coat2_x))
        subs.append(('coat2y = 4;', 'coat2y = %f;', self.coat2_y))
        subs.append(('lc = 0;', 'lc = %f;', self.lc))
        subs.append(
            ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
        subs.append(
            ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))
        subs.append(
            ('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;', self.lc_refine_3))
        subs.append(
            ('lc_refine_4 = lc/1;', 'lc_refine_4 = lc/%f;', self.lc_refine_4))
        subs.append(
            ('lc_refine_5 = lc/1;', 'lc_refine_5 = lc/%f;', self.lc_refine_5))

        return msh_template, msh_name, subs

    def _make_geo_pedestal(self, inc_shape, d_materials):
        subs = ''
        msh_template = 'pedestal'
        self.n_typ_el = 4
        msh_name = self._make_mesh_name(msh_template,
                                        (self.unitcell_x, self.unitcell_y, self.inc_a_x, self.inc_a_y, self.pillar_x, self.pillar_y, self.slab_a_x, self.slab_a_y))

        subs = [('d_in_nm = 100;', 'd_in_nm = %f;', self.unitcell_x)]
        subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.unitcell_y))
        subs.append(('a1 = 20;', 'a1 = %f;', self.inc_a_x))
        subs.append(('a1y = 10;', 'a1y = %f;', self.inc_a_y))
        subs.append(('a1top = 15;', 'a1top = %f;', self.inc_b_x))
        subs.append(('slabx = 80;', 'slabx = %f;', self.slab_a_x))
        subs.append(('slaby = 10;', 'slaby = %f;', self.slab_a_y))
        subs.append(('slabxtop = 60;', 'slabxtop = %f;', self.slab_b_x))
        subs.append(('px = 2;', 'px = %f;', self.pillar_x))
        subs.append(('py = 5;', 'py = %f;', self.pillar_y))
        subs.append(('lc = 0;', 'lc = %f;', self.lc))
        subs.append(
            ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
        subs.append(
            ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))

        return msh_template, msh_name, subs

    # This global passing is horrible.
    # It would be so much nicer to make this a lambda inside _make_onion_geo,
    # but then it is hard to pickle for multiprocessing

    def _make_geo_onion(self, inc_shape, d_materials):
        # TODO: Check:  ax+2*(bx+cx+...) <=unitcell_x

        subs = ''
        msh_template = self.inc_shape
        if inc_shape in ['onion', 'circ_onion']:
            self.n_typ_el = 16
        elif inc_shape in ['onion1', 'circ_onion1']:
            self.n_typ_el = 2
        elif inc_shape in ['onion2', 'circ_onion2']:
            self.n_typ_el = 3
        elif inc_shape in ['onion3', 'circ_onion3']:
            self.n_typ_el = 4

        msh_name = self._make_mesh_name(msh_template,
                                        (self.unitcell_x, self.inc_a_x, self.inc_b_x, self.inc_c_x, self.inc_d_x,
                                         self.inc_e_x, self.inc_f_x, self.inc_g_x, self.inc_h_x, self.inc_i_x, self.inc_j_x, self.inc_k_x, self.inc_l_x, self.inc_m_x, self.inc_n_x,
                                         self.inc_o_x))

        subs = [('d_in_nm = 2000;', 'd_in_nm = %f;', self.unitcell_x)]
        subs.append(('dy_in_nm = 2000;', 'dy_in_nm = %f;', self.unitcell_y))
        subs.append(('a1 = 100;', 'a1 = %f;', self.inc_a_x))
        subs.append(('a2 = 100;', 'a2 = %f;', self.inc_b_x))
        subs.append(('a3 = 100;', 'a3 = %f;', self.inc_c_x))
        subs.append(('a4 = 100;', 'a4 = %f;', self.inc_d_x))
        subs.append(('a5 = 100;', 'a5 = %f;', self.inc_e_x))
        subs.append(('a6 = 100;', 'a6 = %f;', self.inc_f_x))
        subs.append(('a7 = 100;', 'a7 = %f;', self.inc_g_x))
        subs.append(('a8 = 100;', 'a8 = %f;', self.inc_h_x))
        subs.append(('a9 = 100;', 'a9 = %f;', self.inc_i_x))
        subs.append(('a10 = 100;', 'a10 = %f;', self.inc_j_x))
        subs.append(('a11 = 100;', 'a11 = %f;', self.inc_k_x))
        subs.append(('a12 = 100;', 'a12 = %f;', self.inc_l_x))
        subs.append(('a13 = 100;', 'a13 = %f;', self.inc_m_x))
        subs.append(('a14 = 100;', 'a14 = %f;', self.inc_n_x))
        subs.append(('a15 = 100;', 'a15 = %f;', self.inc_o_x))
        subs.append(('lc = 0;', 'lc = %f;', self.lc))
        subs.append(
            ('lc_refine_1 = lc/1;', 'lc_refine_1 = lc/%f;', self.lc_refine_1))
        subs.append(
            ('lc_refine_2 = lc/1;', 'lc_refine_2 = lc/%f;', self.lc_refine_2))

        # make wg frame for plots. factor 2 because first element is the diameter,
        # remainder are the annular thicknesses
        global g_onion_layer_thicknesses
        g_onion_layer_thicknesses = (self.inc_a_x/2, self.inc_b_x, self.inc_c_x, self.inc_d_x, self.inc_e_x,
                                     self.inc_f_x, self.inc_g_x, self.inc_h_x, self.inc_i_x, self.inc_j_x,
                                     self.inc_k_x, self.inc_l_x, self.inc_m_x, self.inc_n_x, self.inc_o_x)

        # TODO: FIX ME: Bring this funcitonality back
        # self.mpl_wg_frame_drawer = _my_frm_onion

        return msh_template, msh_name, subs

    def _make_geo_trapezoidal_rib(self, inc_shape, d_materials):
        # inc_a_x - width of the top of the rib
        # inc_a_y - height of the top of the rib
        # slab_a_x - width of the middle of the rib
        # slab_a_y - height of the buried part of the rib

        msh_template = 'trapezoidal_rib'
        self.n_typ_el = 4
        msh_name = self._make_mesh_name(msh_template,
                                        (self.unitcell_x, self.inc_a_x, self.inc_a_y, self.slab_a_x, self.slab_a_y))

        subs = [('top_rib_width = 600.0;', "top_rib_width = %f;", self.inc_a_x)]
        subs.append(('mid_rib_width = 900.0;',
                    "mid_rib_width = %f;", self.slab_a_x))
        subs.append(('rib_height = 500.0;', "rib_height = %f;", self.inc_a_y))
        subs.append(('slab_thickness = 300.0;',
                    "slab_thickness = %f;", self.slab_a_y))
        subs.append(('lc = 0.020000;', "lc = %f;", self.lc))
        subs.append(('lc_refine_1 = lc/10.0;',
                    "lc_refine_1 = lc/%f;", self.lc_refine_1))
        subs.append(('lc_refine_2 = lc/5.0;',
                    "lc_refine_2 = lc/%f;", self.lc_refine_2))

        return msh_template, msh_name, subs

    def get_mail_data(self):
        '''Returns list of lines in the NumBAT mesh format .mail file'''
        return self.mail_data

    def make_mesh(self, d_materials):
        ''' Take the parameters specified in python and make a Gmsh FEM mesh.
            Creates a .geo and .msh file from the .geo template,
            then uses Fortran conv_gmsh routine
            to convert .msh into .mail, which is used in NumBAT FEM routine.
        '''

        print('Building mesh')

        # full path to backend directory that this code file is in
        this_directory = os.path.dirname(os.path.realpath(__file__))
        # msh_location_in = os.path.join(this_directory, 'fortran', 'msh', '')  # msh directory inside backend
        # msh directory inside backend
        msh_location_in = os.path.join(this_directory, 'msh', '')

        self.msh_location_in = msh_location_in
        self.msh_location_out = os.path.join(msh_location_in, 'build', '')

        if not os.path.exists(self.msh_location_out):
            os.mkdir(self.msh_location_out)

        if self.inc_shape in ['circular', 'rectangular', 'twoincl']:
            msh_template, msh_name, subs = self._make_geo_circular_rectangular(
                self.inc_shape, d_materials)

        elif self.inc_shape in ['rib']:
            msh_template, msh_name, subs = self._make_geo_rib(
                self.inc_shape, d_materials)

        elif self.inc_shape in ['rib_coated']:
            msh_template, msh_name, subs = self._make_geo_rib_coated(
                self.inc_shape, d_materials)

        elif self.inc_shape in ['rib_double_coated']:
            msh_template, msh_name, subs = self._make_geo_rib_double_coated(
                self.inc_shape, d_materials)

        elif self.inc_shape in ['slot', 'slot_coated']:
            msh_template, msh_name, subs = self._make_geo_slot(
                self.inc_shape, d_materials)

        elif self.inc_shape in ['pedestal']:
            msh_template, msh_name, subs = self._make_geo_pedestal(
                self.inc_shape, d_materials)

        elif self.inc_shape in ['onion', 'onion1', 'onion2', 'onion3',
                                'circ_onion', 'circ_onion1', 'circ_onion2', 'circ_onion3']:
            msh_template, msh_name, subs = self._make_geo_onion(
                self.inc_shape, d_materials)

        elif self.inc_shape in ['trapezoidal_rib']:
            msh_template, msh_name, subs = self._make_geo_trapezoidal_rib(
                self.inc_shape, d_materials)

        else:
            raise NotImplementedError(f"\n Selected inc_shape = '{self.inc_shape}' "
                                      'is not currently implemented. \nPlease make a mesh with gmsh and '
                                      'consider contributing this to NumBAT via github.')

        # check for gmsh install

        if platform.system() == 'Darwin':  # TODO: do better testing for MAC
            if not Path('/Applications/Gmsh.app/Contents/MacOS/gmsh').is_file():
                reporting.report_and_exit(
                    'Gmsh does not seem to be installed at /Applications/Gmsh.app/Contents/MacOS/gmsh.')

        else:
            if shutil.which('gmsh') is None:
                reporting.report_and_exit(
                    'Gmsh does not seem to be installed.')

        if self._build_new_mesh():
            geo = self._apply_geo_subs(msh_template, subs)

            fname = self.msh_location_out + msh_name

            with open(fname + '.geo', 'w') as fout:
                fout.write(geo)

            # Convert our Gmsh .geo file into Gmsh .msh and then NumBAT .mail
            err_no, err_msg = NumBAT.conv_gmsh(fname)
            if err_no != 0:

                s = 'Terminating after Fortran error in processing .geo file "%s.geo".' % fname
                if len(err_msg):
                    s += '\nMessage was:\n %s' % err_msg
                s += f'''

                Is the mesh template file "backend/fortran/msh/{msh_template}.geo" designed correctly?
                To help diagnose the problem, try viewing the generated mesh file in gmsh by running:
                  gmsh {fname}.geo'''
                reporting.report_and_exit(s)

        self.mesh_file = fname + '.mail'
        # TODO: curently used onyl to generate filenames for plot_mesh. Needed?
        self.msh_name = msh_name

        # read in first line giving number of msh points and elements
        with open(self.mesh_file) as f:
            self.mail_data = f.readlines()

    def plot_mesh(self, outpref):
        '''Visualise mesh with gmsh and save to a file.'''

        # Manipulate scripts in backend/fortran/build
        # Writes final png file to user directory

        nbapp = numbat.get_NumBATApp()

        gmsh_exe = nbapp.path_gmsh()
        print('gmshexe', gmsh_exe)

        cwd = os.getcwd()
        outprefix = os.path.join(cwd, outpref)
        tdir = tempfile.TemporaryDirectory()
        tmpoutpref = os.path.join(tdir.name, outpref)

        conv_tmp = open(self.msh_location_in + 'geo2png.scr', 'r').read()
        conv = conv_tmp.replace('tmp', tmpoutpref + '-mesh_geom')

        fn_scr = self.msh_location_out + self.msh_name + '_geo2png.scr'

        cmd = [gmsh_exe, self.msh_name + '.geo', self.msh_name + '_geo2png.scr']

        with open(fn_scr, 'w') as fout:
            fout.write(conv)
        subprocess.Popen(cmd, cwd=self.msh_location_out)

        os.wait()

        conv_tmp = open(self.msh_location_in + 'msh2png.scr', 'r').read()
        conv = conv_tmp.replace('tmp', tmpoutpref + '-mesh_nodes')

        fn_scr = self.msh_location_out + self.msh_name + '_msh2png.scr'
        # cmd = ' '.join(['gmsh', f' -log {outprefix}-gmshlog.log',
        #                self.msh_name + '.msh',
        #                self.msh_name + '_msh2png.scr' ])
        cmd = [gmsh_exe, self.msh_name + '.msh', self.msh_name + '_msh2png.scr']
        with open(fn_scr, 'w') as fout:
            fout.write(conv)

        subprocess.Popen(cmd, cwd=self.msh_location_out)

        (pid, status) = os.wait()

        nbtools.join_figs(tmpoutpref+'-mesh_geom.png',
                          tmpoutpref+'-mesh_nodes.png',
                          outprefix+'-mesh.png')

    def check_mesh(self):
        '''Visualise geometry and mesh with gmsh.'''

        gmsh_cmd = 'gmsh ' + self.msh_location_out + self.msh_name + '.geo'
        os.system(gmsh_cmd)

        gmsh_cmd = 'gmsh ' + self.msh_location_out + self.msh_name + '.msh'
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
