# objects.py is a subroutine of NumBAT. It contains the Struct
# objects that represent the structure being simulated.

# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

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

import copy

import numbat
import reporting
from nbtypes import SI_nm
from numbattools import run_subprocess

import meshing.templates as mshtemplates

import materials
from modecalcs import em_mode_calculation, ac_mode_calculation
from materialprops import OpticalProps, ElasticProps

import plotting.gmsh as pltgmsh
import plotting.profiles as pltprof

from fortran import nb_fortran


class Structure:
    """ Represents the geometry and  material properties (elastic and optical) of a waveguide structure.

        Args:
            domain_x  (float):
                The horizontal period of the unit cell in nanometers.

            domain_y  (float):
                The vertical period of the unit cell in nanometers. If None, domain_y = domain_x.

            inc_shape  (str):
                Shape of inclusions that have template mesh, currently:
                 ``circular`` ``rectangular`` ``slot`` ``rib`` ``slot_coated`` ``rib_coated``
                 ``rib_double_coated`` ``pedestal`` ``onion`` ``onion2`` ``onion3``.
                 rectangular is default.

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

            lc_bkg  (float): Length constant of meshing of background medium
                (smaller = finer mesh)

            lc_refine_1  (float): factor by which lc_bkg will be reduced on inclusion
                surfaces; lc_surface = lc_bkg / lc_refine_1. Larger lc_refine_1 = finer mesh.

            lc_refine_2-6'  (float): factor by which lc_bkg will be reduced on chosen
                surfaces; lc_surface = lc_bkg / lc_refine_2. see relevant .geo files.

    """

    def _clean_and_handle_args_and_parameters(self, *largs, **kwargs):

        if kwargs.get('direct_call', 0):
            reporting.register_warning(
                'Calling objects.Structure directly is deprecated. Please switch to calling nbapp.make_structure()')

        if 'inc_shape' not in kwargs and len(largs) == 0:
            reporting.report_and_exit('Must provide an inc_shape argument to make_structure()')

        self.all_params = copy.deepcopy(kwargs)
        self.all_params['lc'] = self.all_params['lc_bkg']   # aim to move to pure lc_bkg and delete this

        del self.all_params['direct_call']

        if 'inc_shape' in self.all_params: del self.all_params['inc_shape']
        self.inc_shape = largs[0] if len(largs)>0 else kwargs['inc_shape']

        # Support old calling convention
        if 'domain_x' not in kwargs:
            self.all_params['domain_x'] = largs[1] if len(largs)>1 else float(kwargs['unitcell_x'])
        if 'domain_y' not in kwargs:
            self.all_params['domain_y'] = largs[2] if len(largs)>2 else float(kwargs['unitcell_y'])

        self.domain_x = self.all_params['domain_x']
        self.domain_y = self.all_params['domain_y']

        #TODO: everything breaks is 'loss' is false, because v_refindex become real which breaks passing to fortran
        self.loss = float(kwargs.get('loss', True))       # TODO: Boolean. Needs better name


        if len(largs)>3: self.all_params['inc_a_x'] = largs[3]
        if len(largs)>4: self.all_params['inc_a_y'] = largs[4]




        # needed before building
        self.symmetry_flag = int(kwargs.get('symmetry_flag', 0))

        # Minor stuff
        self.shift_em_x = 0  # user requested offsets to coord-system
        self.shift_em_y = 0
        self.shift_ac_x = 0  # user requested offsets to coord-system
        self.shift_ac_y = 0


    def __init__(self, *largs, **kwargs):
        """Initialise the Structure object."""

        numbat.assert_numbat_object_created()

        self._clean_and_handle_args_and_parameters(*largs, **kwargs)

        # locations of gmsh input and output files
        self.msh_location_in = numbat.NumBATApp().path_mesh_templates()
        self.msh_location_out = self.msh_location_in / 'build'

        if not self.msh_location_out.is_dir():
            self.msh_location_out.mkdir()


        mat_vac = materials.make_material('Vacuum')
        self.d_materials = {'bkg': kwargs['material_bkg']}  # this one must come first

        # fill up materials with vacuums
        for letter in range(ord('a'), ord('r')+1):
            self.d_materials[chr(letter)] = mat_vac

        # overwrite vacuums for those materials provided by user
        for k,v in kwargs.items():
            if k.startswith('material_') and k !='material_bkg':
                self.d_materials[k[9:]] = v

        n_mats_em = self._build_waveguide_geometry()

        # Build the whole mesh (A major step involving Fortran)
        self._build_mesh()

        self.optical_props = OpticalProps(self.d_materials, n_mats_em, self.loss)

        self.elastic_props = ElasticProps(self, self.symmetry_flag)


    def get_material(self, k):
        return self.d_materials[k]

    def get_optical_materials(self):
        return list(self.d_materials.values())[: self.optical_props.n_mats_em]

    def set_xyshift_em(self, x, y):
        # Sets shift in grid from user perspective in nm
        self.shift_em_x = x * SI_nm
        self.shift_em_y = y * SI_nm

    def set_xyshift_ac(self, x, y):
        self.shift_ac_x = x * SI_nm
        self.shift_ac_y = y * SI_nm

    def using_linear_elements(self):
        return self.inc_shape in self.linear_element_shapes

    def using_curvilinear_elements(self):
        return self.inc_shape in self.curvilinear_element_shapes


    def _build_waveguide_geometry(self):
        """ Take the parameters specified in python and make a User waveguide class.
        """

        # Find and instantiate the class that defines this waveguide model
        wg_cls = mshtemplates.get_waveguide_template_class(self.inc_shape)
        wg_geom = wg_cls(self.all_params, self.d_materials)
        wg_geom.init_geometry()   # run the user code to set up the geometry

        # This is number of distinct materials == element types materials declared by the template
        n_mats_em = wg_geom._num_materials

        assert n_mats_em > 0, 'No active materials defined in the waveguide geometry.'

        wg_geom.check_parameters(self.all_params)
        wg_geom.validate_dimensions()

        self.wg_geom = wg_geom

        self.linear_element_shapes = []
        self.curvilinear_element_shapes = []

        if wg_geom.is_curvilinear:
            self.curvilinear_element_shapes.append(wg_geom._shape_name)
        else:
            self.linear_element_shapes.append(wg_geom._shape_name)

        return n_mats_em

    def _build_mesh(self):
        """Instantiates generic template gmsh file to aspecific gmsh then runs conv_gmsh
        to generate the NumBAT .mail file."""

        print(' Calling Gmsh: ', end='')

        msh_tmpl_fname = self.wg_geom.gmsh_template_filename()  # template we're reading from
        msh_inst_fname = self.wg_geom.get_instance_filename()      # unique name for this instance

        # create string in .geo format from the template file with all parameters adjusted for our design
        geo_str = self.wg_geom.make_geometry(self.msh_location_in)

        path_geo_out = self.msh_location_out / msh_inst_fname

        with open(str(path_geo_out) + '.geo', 'w') as fout:
            fout.write(geo_str)

        # Convert our Gmsh .geo file into Gmsh .msh
        gmsh_exe =  numbat.NumBATApp().path_gmsh()
        args =f' -2 -order 2 -v 0 -o {msh_inst_fname}.msh {msh_inst_fname}.geo'
        cmd = [str(gmsh_exe)]
        cmd.extend(args.split())
        run_subprocess(cmd, 'Gmsh', cwd=self.msh_location_out)

        # And now to NumBAT .mail format
        assertions_on = False
        err_no, err_msg = nb_fortran.conv_gmsh(str(path_geo_out), assertions_on)
        if err_no:
            s = f'Terminating after Fortran error in processing .geo file "{path_geo_out}.geo".'
            if len(err_msg):
                s += f'\nMessage was:\n {err_msg}'
            s += f"""

            Is the mesh template file "backend/fortran/msh/{msh_tmpl_fname}_msh_template.geo" designed correctly?
            To help diagnose the problem, try viewing the generated mesh file in gmsh by running:
                gmsh {path_geo_out}.geo"""
            reporting.report_and_exit(s)

        self.mesh_mail_fname = str(path_geo_out) + '.mail'
        # TODO: curently used onyl to generate filenames for plot_mesh. Needed? Fix the filenames.
        self.msh_name = msh_inst_fname

    def plot_mail_mesh(self, outpref):
        """Visualise the mesh in .mail format."""

        pltgmsh.plot_mail_mesh(self.mesh_mail_fname, outpref)

    def plot_mesh(self, outpref, combo_plot=True):
        """Visualise mesh with gmsh and save to a file.

        If combo_plot==False, wire frame and mesh node plots are stored in separate files."""

        return pltgmsh.plot_mesh(self.msh_location_in, self.msh_location_out,
                          self.msh_name, outpref, combo_plot)

    def check_mesh(self):
        """Visualise geometry and mesh with gmsh."""

        pltgmsh.check_mesh(self.msh_location_out, self.msh_name)


    def calc_EM_modes(self, num_modes, wl_nm, n_eff, Stokes=False, debug=False,
                      **args):
        """ Run a simulation to find the Structure's EM modes.

            Args:
                num_modes  (int): Number of EM modes to solve for.

                wl_nm  (float): Wavelength of EM wave in vacuum in nanometres.

                n_eff  (float): Guesstimated effective index of
                    fundamental mode, will be origin of FEM search.

            Returns:
                ``EMSimResult`` object
        """
        return em_mode_calculation(self, num_modes, wl_nm, n_eff, Stokes, debug, **args)

    def calc_AC_modes(self, num_modes, q_AC,
                      shift_Hz=None, EM_sim=None, bcs=None, debug=False, **args):
        """ Run a simulation to find the Structure's acoustic modes.

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
                ``ACSimResult`` object
        """

        return ac_mode_calculation(self, num_modes, q_AC, shift_Hz, EM_sim, bcs, debug, **args)


    # def _make_refindex_plotter(self, as_epsilon, n_points):


    #     v_neffeps = self.optical_props.v_refindexn  # mapping from material index to refractive index

    #     if as_epsilon:
    #         v_neffeps = v_neffeps**2
    #         nm_eng = 'Dielectric constant'
    #         nm_math=r'$\epsilon(\vec x)$'
    #         fname_suffix='dielectric_constant'
    #     else:
    #         nm_eng = 'Refractive index'
    #         nm_math=r'$n(\vec x)$'
    #         fname_suffix='refractive_index'

    #     fsfp = femmesh.FEMScalarFieldPlotter(self, n_points)

    #     unit=''

    #     fsfp.setup_scalar_properties(nm_eng, unit, nm_math, fname_suffix)
    #     fsfp.fill_quantity_by_material_index(v_neffeps)

    #     return fsfp

    def plot_refractive_index_profile(self, pref):
        """Draw 2D plot of refractive index profile."""
        pl_ref = self.get_structure_plotter_refractive_index()
        pl_ref.make_plot_2D(pref)

    def get_structure_plotter_refractive_index(self, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the refractive index profile."""
        return pltprof._make_dielectric_plotter(self, self.optical_props, False, n_points)


    def get_structure_plotter_epsilon(self, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the dielectric constant profile."""
        return pltprof._make_dielectric_plotter(self, self.optical_props, True, n_points)

    def get_structure_plotter_stiffness(self, c_I, c_J, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the elastic stiffness."""
        return pltprof.get_structure_plotter_stiffness(self, c_I, c_J, n_points)


    def get_structure_plotter_acoustic_velocity(self, vel_index=0, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the elastic acoustic phase speed.

        Args:
            vel_index (0,1,2): Index of the elastic mode phase speed to plot.

        Currently only works for isotropic materials."""

        return pltprof.get_structure_plotter_acoustic_velocity(self, self.d_materials,
                                                               vel_index, n_points)


    def plot_refractive_index_profile_rough(self, prefix, n_points = 200, as_epsilon=False):
        """ Draws refractive index profile by primitive sampling, not proper triangular mesh sampling"""

        pltprof.plot_refractive_index_profile_rough(
            self.mesh_mail_fname,
            self.d_materials,
            prefix, n_points , as_epsilon)

