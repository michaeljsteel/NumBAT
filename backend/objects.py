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

import subprocess
import copy
import traceback
#import itertools


import tempfile
from pathlib import Path
import json
import importlib
import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker

import numpy as np
import scipy.interpolate


import numbat
import reporting
from nbtypes import SI_nm
from numbattools import run_subprocess, f2f_with_subs

import materials
from mode_calcs import em_mode_calculation, ac_mode_calculation
import nbgmsh
import plotting
import plottools
import femmesh
from materialprops import OpticalProps, ElasticProps

from fortran import nb_fortran

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

        if 'plt_mesh' in kwargs:
            print("\n Warning: Option 'plt_mesh' is deprecated. Call method .plot_mesh() on your Struct object.\n\n")

        if 'check_mesh' in kwargs:
            print("\n Warning: Option 'check_mesh' is deprecated. Call method .check_mesh() on your Struct object.\n\n")

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

        mat_vac = materials.make_material('Vacuum')
        self.d_materials = {'bkg': kwargs['material_bkg']}  # this one must come first

        # fill up materials with vacuums
        for letter in range(ord('a'), ord('r')+1):
            self.d_materials[chr(letter)] = mat_vac

        # overwrite vacuums for those materials provided by user
        for k,v in kwargs.items():
            if k.startswith('material_') and k !='material_bkg':
                self.d_materials[k[9:]] = v


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

        n_mats_em = self._build_waveguide_geometry()

        # Build the whole mesh (A major step involving Fortran)
        self._build_mesh()


        self.optical_props = OpticalProps(list(self.d_materials.values()), n_mats_em, self.loss)

        # construct list of materials with nonzero density, ie with acoustic properties likely defined
        # Any material not given v_acoustic_mats assumed to be vacuum.
        #v_acoustic_mats = [m for m in self.d_materials.values() if m.has_elastic_properties()]

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
        """ Take the parameters specified in python and make a Gmsh FEM mesh.
            Creates a .geo and .msh file from the .geo template,
            then uses Fortran conv_gmsh routine
            to convert .msh into .mail, which is used in NumBAT FEM routine.
        """

        # full path to backend directory that this code file is in
        this_directory = Path(__file__).resolve().parent

        # locations of gmsh input and output files
        self.msh_location_in = this_directory / 'msh'
        self.msh_location_out = self.msh_location_in / 'build'


        if not self.msh_location_out.is_dir():
            self.msh_location_out.mkdir()

        self.wg_geom = None

        for wg in g_waveguide_templates:
            if self.inc_shape in wg['inc_shape']:  # is the desired shape supported by this template class?

                # Instantiate the class that defines this waveguide model
                wg_geom = wg['wg_template_cls'](self.all_params, self.d_materials)

                wg_geom.init_geometry()   # run the user code to set up the geometry

                n_mats_em = wg_geom._num_materials  # This is number of distinct materials == element types materials declared by the template

                assert n_mats_em > 0, 'No active materials defined in the waveguide geometry.'

                break

        else:  # didn't find required wg
            raise NotImplementedError(f"\n Selected inc_shape = '{self.inc_shape}' "
                                      'is not currently implemented. \nPlease make a mesh with gmsh and '
                                      'consider contributing this to NumBAT via github.')

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

        msh_template = self.wg_geom.gmsh_template_filename()  # template we're reading from
        msh_fname = self.wg_geom.get_instance_filename()      # unique name for this instance

        # create string in .geo format from the template file with all parameters adjusted for our design
        geo_str = self.wg_geom.make_geometry(numbat.NumBATApp().path_mesh_templates())
        fname = Path(self.msh_location_out) / msh_fname
        with open(str(fname) + '.geo', 'w') as fout:
            fout.write(geo_str)

        # Convert our Gmsh .geo file into Gmsh .msh
        gmsh_exe =  numbat.NumBATApp().path_gmsh()
        args =f' -2 -order 2 -v 0 -o {msh_fname}.msh {msh_fname}.geo'
        cmd = [gmsh_exe]
        cmd.extend(args.split())
        comp_stat = subprocess.run(cmd, cwd=self.msh_location_out)
        if comp_stat.returncode:
            tcmd = ' '.join(cmd)
            reporting.report_and_exit(f'Gmsh call failed executing: "{tcmd}".')

        # And now to NumBAT .mail format
        assertions_on = False
        err_no, err_msg = nb_fortran.conv_gmsh(str(fname), assertions_on)
        if err_no:
            s = f'Terminating after Fortran error in processing .geo file "{fname}.geo".'
            if len(err_msg):
                s += f'\nMessage was:\n {err_msg}'
            s += f"""

            Is the mesh template file "backend/fortran/msh/{msh_template}_msh_template.geo" designed correctly?
            To help diagnose the problem, try viewing the generated mesh file in gmsh by running:
                gmsh {fname}.geo"""
            reporting.report_and_exit(s)

        self.mesh_mail_fname = str(fname) + '.mail'
        # TODO: curently used onyl to generate filenames for plot_mesh. Needed? Fix the filenames.
        self.msh_name = msh_fname

    def plot_mail_mesh(self, outpref):
        """Visualise the mesh in .mail format."""
        path = numbat.NumBATApp().outpath()
        mail_data = nbgmsh.MailData(self.mesh_mail_fname)
        mail_data.plot_mesh(path)

    def plot_mesh(self, outpref):
        """Visualise mesh with gmsh and save to a file."""

        # Manipulate scripts in backend/fortran/build
        # Writes final png file to user directory

        nbapp = numbat.NumBATApp()
        gmsh_exe = nbapp.path_gmsh()

        outprefix = Path(numbat.NumBATApp().outdir(), outpref)
        tdir = tempfile.TemporaryDirectory()
        tmpoutpref = str(Path(tdir.name, outpref))

        # Make the wire frame image
        fn_in = Path(self.msh_location_in) / 'geo2png.scr'
        fn_out = Path(self.msh_location_out) / (self.msh_name + '_geo2png.scr')
        f2f_with_subs(fn_in, fn_out, {'tmp': str(tmpoutpref) + '-entities'})
        cmd = [gmsh_exe, self.msh_name + '.geo', fn_out.name]
        run_subprocess(cmd, 'Gmsh', cwd=self.msh_location_out)

        # Make the mesh image
        fn_in = Path(self.msh_location_in) / 'msh2png.scr'
        fn_out = Path(self.msh_location_out) / (self.msh_name + '_msh2png.scr')
        f2f_with_subs(fn_in, fn_out, {'tmp': str(tmpoutpref) + '-mesh_nodes'})

        cmd = [gmsh_exe, self.msh_name + '.msh', fn_out.name]
        run_subprocess(cmd, 'Gmsh', cwd=self.msh_location_out)

        # Join the two images
        plottools.join_figs([tmpoutpref+'-entities.png',
                          tmpoutpref+'-mesh_nodes.png',],
                          str(outprefix)+'-mesh.png',
                          clip=(60,50,60,50))



    def check_mesh(self):
        """Visualise geometry and mesh with gmsh."""

        nbapp = numbat.NumBATApp()
        gmsh_exe = str(nbapp.path_gmsh())

        gmsh_cmd = [gmsh_exe, f'{self.msh_location_out}/{self.msh_name}.geo']
        run_subprocess(gmsh_cmd, 'Gmsh', cwd=self.msh_location_out)

        gmsh_cmd = [gmsh_exe, f'{self.msh_location_out}/{self.msh_name}.msh']
        run_subprocess(gmsh_cmd, 'Gmsh', cwd=self.msh_location_out)


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


    def _make_refindex_plotter(self, as_epsilon, n_points):


        v_neffeps = self.optical_props.v_refindexn  # mapping from material index to refractive index

        if as_epsilon:
            v_neffeps = v_neffeps**2
            nm_eng = 'Dielectric constant'
            nm_math=r'$\epsilon(\vec x)$'
            fname_suffix='dielectric_constant'
        else:
            nm_eng = 'Refractive index'
            nm_math=r'$n(\vec x)$'
            fname_suffix='refractive_index'

        fsfp = femmesh.FEMScalarFieldPlotter(self, n_points)

        unit=''

        fsfp.setup_scalar_properties(nm_eng, unit, nm_math, fname_suffix)
        fsfp.fill_quantity_by_material_index(v_neffeps)

        return fsfp

    def plot_refractive_index_profile(self, pref):
        """Draw 2D plot of refractive index profile."""
        pl_ref = self.get_structure_plotter_refractive_index()
        pl_ref.make_plot_2D(pref)

    def get_structure_plotter_refractive_index(self, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the refractive index profile."""
        return self._make_refindex_plotter(False, n_points)

    def get_structure_plotter_epsilon(self, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the dielectric constant profile."""
        return self._make_refindex_plotter(True, n_points)

    def get_structure_plotter_stiffness(self, c_I, c_J, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the elastic stiffness."""
        if c_I not in range(1,7) or c_J not in range(1,7):
            reporting.report_and_exit('Stiffness tensor indices c_I, c_J must be in the range 1..6.')

        v_stiff = np.zeros(5) # fill me

        fsfp = femmesh.FEMScalarFieldPlotter(self, n_points)
        qname = 'Stiffness $c_{'+f'{c_I},{c_J}' +'}$'
        suffname = f'stiffness_c_{c_I}{c_J}'
        fsfp.set_quantity_name(qname, suffname)
        fsfp.fill_scalar_by_material_index(v_stiff)

    def get_structure_plotter_acoustic_velocity(self, vel_index=0, n_points=500):
        """Make plotter for arbitrary 1D and 2D slices of the elastic acoustic phase speed.

        Args:
            vel_index (0,1,2): Index of the elastic mode phase speed to plot.

        Currently only works for isotropic materials."""

        v_mats = list(self.d_materials.values())
        v_acvel = np.zeros([len(v_mats),3])

        for i in range(len(v_mats)): # get all 3 phase velocities for each material
            if v_mats[i].has_elastic_properties():
                v_acvel[i,:] = v_mats[i].Vac_phase()

        fsfp = femmesh.FEMScalarFieldPlotter(self, n_points)
        fsfp.setup_vector_properties(3, 'Elastic velocity', '[km/s]', r'$v_i$',
                                     [r'$v_0$', r'$v_1$', r'$v_2$'],
                                     'elastic_velocity', ['v0', 'v1', 'v2'])

        fsfp.fill_quantity_by_material_index(v_acvel)
        return fsfp


    def plot_refractive_index_profile_rough(self, prefix, n_points = 200, as_epsilon=False):
        """ Draws refractive index profile by primitive sampling, not proper triangular mesh sampling"""

        print('\n\nPlotting ref index')

        mail = nbgmsh.MailData(self.mesh_mail_fname)
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
            v_refindex[i] = v_elt_refindex[elt-1]  # the type of element is labelled by gmsh from 1.

        # Now we have an irregular x,y,n array to interpolate onto.


        # Construct a regular rect array with n_pts_x * n_pts_y ~ n_points**2
        # and with approximately square pixels
        x_min = np.min(v_x)
        x_max = np.max(v_x)
        y_min = np.min(v_y)
        y_max = np.max(v_y)

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

        #v_regindex = np.where(v_regindex==0, 1, v_regindex)
        v_regindex = np.nan_to_num(v_regindex, nan=1.0)

        if as_epsilon:
            v_regindex = v_regindex**2
            fig.suptitle('Dielectric constant')
        else:
            fig.suptitle('Refractive index')

        cmap='cool'
        cf=ax.imshow(v_regindex, cmap=cmap, vmin=1.0, vmax=np.nanmax(v_regindex), origin='lower',
                     extent = [x_min, x_max, y_min, y_max])
        #cf=ax.contourf(m_regx, m_regy, v_regindex, cmap=cmap, vmin=1.0, vmax=np.nanmax(v_regindex))
        ax.set_xlabel(r'$x$ [μm]')
        ax.set_ylabel(r'$y$ [μm]')
        cb = fig.colorbar(cf)
        cf.set_clim(1,np.nanmax(v_regindex))
        cb.outline.set_linewidth(.5)
        cb.outline.set_color('gray')


        plotting.save_and_close_figure(fig, prefix+'refn.png')




def print_waveguide_help(inc_shape):
    for wg in g_waveguide_templates:
        if inc_shape in wg['inc_shape']:  # is the desired shape supported by this template class?
            #found = True

            # Instantiate the class that defines this waveguide model
            wg_geom = wg['wg_template_cls'](None, None)
            wg_geom.init_geometry()
            print(wg_geom.get_parameter_help_summary())
