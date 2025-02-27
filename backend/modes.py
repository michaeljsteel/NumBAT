
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


from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

import numbat
from nbtypes import FieldType, FieldCode, FieldTag, SI_um, SI_vacuum_impedance_Z0
import reporting

from numbattools import int2D_trapz
import femmesh
import plotmodes
import plotmoderaw
import plottools



class ModeInterpolator:
    '''Helper class for plotting modes.
       Factors common info from Simulation that each mode can draw on, but we only need to do once for each Sim.
       '''

    def __init__(self, simresult):
        self.sim_result = simresult

        self.setup_for_2D_n_pts = 0
        self.setup_for_1D = ()

        # Points for the plot grid in absolute coords and units of meters
        self.xy_raw = {}

        # Points for the plot grid in shifted coords and units of microns
        self.xy_out = {}

        # x-axis for any 1D cut plots
        self.v_x_axis = None


        self.interper_f_2D = None
        self.interper_f_1D = None

        self.zero_arrays()

    def zero_arrays(self):
        self.interper_f_2D = None
        self.interper_f_1D = None

        self.xy_raw = {}
        self.xy_out = {}

    def cleanup(self):
        if self.interper_f_2D is not None: del self.interper_f_2D
        if self.interper_f_1D is not None: del self.interper_f_1D

        self.zero_arrays()


    def interpolate_mode_i(self, md, field_type, dims=2):
        # construct the meshed field from fortran solution

        simres = self.sim_result

        (v_Fx6p, v_Fy6p, v_Fz6p, v_Fa6p) = simres.get_mode_fields_on_fem_mesh(
            md, field_type)

        interper_f = self.interper_f_2D if dims == 2 else self.interper_f_1D

        # Always need these ones.
        m_ReFx = interper_f(v_Fx6p.real)
        m_ReFy = interper_f(v_Fy6p.real)
        m_ImFz = interper_f(v_Fz6p.imag)
        m_AbsF = interper_f(v_Fa6p)

        # often not needed for plotting, but are used for measuring fractions. (Could fix taht?)
        m_ImFx = interper_f(v_Fx6p.imag)
        m_ImFy = interper_f(v_Fy6p.imag)
        m_ReFz = interper_f(v_Fz6p.real)

        d_fields = {'Fxr': m_ReFx, 'Fxi': m_ImFx,
                    'Fyr': m_ReFy, 'Fyi': m_ImFy,
                    'Fzr': m_ReFz, 'Fzi': m_ImFz,
                    'Fa': m_AbsF}

        if field_type == FieldType.EM_H:  # scale H fields by Z0 to get common units and amplitude with E
            for m_F in d_fields.values():   # Do this when they are first made
                m_F *= SI_vacuum_impedance_Z0

        return d_fields


    def _choose_plot_points(self, n_pts):
        '''Picks actual data points for the plot grid based on requested resolution.'''
        self.setup_for_2D_n_pts = n_pts

        fm = self.sim_result.fem_mesh

        x_min, x_max, y_min, y_max = fm.get_xy_limits()

        area = abs((x_max-x_min)*(y_max-y_min))
        n_pts_x = int(n_pts*abs(x_max-x_min)/np.sqrt(area))
        n_pts_y = int(n_pts*abs(y_max-y_min)/np.sqrt(area))

        # Now use the coords user would like to think in
        shiftx, shifty = self.sim_result.get_xyshift()

        # These are the actual x and y domains of the final plots
        v_x = np.linspace(x_min, x_max, n_pts_x)
        v_y = np.linspace(y_min, y_max, n_pts_y)
        m_X, m_Y = np.meshgrid(v_x, v_y)

        # raw means in meters and with no position shifts
        self.xy_raw = {'v_x': v_x, 'v_y': v_y, 'm_x': m_X, 'm_y': m_Y}

        v_x_out = (v_x + shiftx) / SI_um
        v_y_out = (v_y + shifty) / SI_um
        m_X_out, m_Y_out = np.meshgrid(v_x_out, v_y_out)

        self.xy_out = {'v_x': v_x_out, 'v_y': v_y_out,
                       'm_x': m_X_out, 'm_y': m_Y_out}

        print('Structure has raw domain (x,y) = '
              + f' [{v_x[0]/SI_um:.5f}, {v_x[-1]/SI_um:.5f}] x [{v_y[0] /SI_um:.5f}, {v_y[-1]/SI_um:.5f}] (μm),'
              + "\n             mapped to (x',y') = "
              + f' [{v_x_out[0]:.5f}, {v_x_out[-1]:.5f}] x [{v_y_out[0]:.5f}, {v_y_out[-1]:.5f}] (μm)')
        return shiftx, shifty



    def define_plot_grid_2D(self, n_pts=501):
        '''Define interpolation plotting grids for a nominal n_pts**2 points distributed evenly amongst x and y.'''

        if self.setup_for_2D_n_pts == n_pts:
            return  # only need to repeat if the grid density changes

        fm = self.sim_result.fem_mesh

        shiftx, shifty = self._choose_plot_points(n_pts)

        v_x_flat = self.xy_raw['m_x'].flatten('F') - shiftx
        v_y_flat = self.xy_raw['m_y'].flatten('F') - shifty

        nx, ny = len(self.xy_out['v_x']), len(self.xy_out['v_y'])

        self.interper_f_2D = fm.make_interpolator_for_grid(v_x_flat, v_y_flat, nx, ny)

    def define_plot_grid_1D(self, s_cut, val1, val2, n_pts):

        if self.setup_for_1D == (s_cut, val1, val2, n_pts):
            return  # no update to interpolator needed

        self.setup_for_1D=(s_cut, val1, val2, n_pts)

        fm = self.sim_result.fem_mesh
        x_min, x_max, y_min, y_max = fm.get_xy_limits()


        match s_cut.lower():
            case 'x':
                v_y_flat = np.linspace(y_min, y_max, n_pts)
                v_x_flat = val1 * SI_um + np.zeros(n_pts)
            case 'y':
                v_x_flat = np.linspace(x_min, x_max, n_pts)
                v_y_flat = val1 * SI_um + np.zeros(n_pts)
            case 'line':
                v_x_flat = np.linspace(val1[0]* SI_um, val2[0]* SI_um, n_pts)
                v_y_flat = np.linspace(val1[1]* SI_um, val2[1]* SI_um, n_pts)

        nx = len(v_x_flat)
        ny = 1
        self.interper_f_1D = fm.make_interpolator_for_grid(
            v_x_flat, v_y_flat, nx, ny)

        shiftx, shifty = self.sim_result.get_xyshift()
        v_x_flat -= shiftx
        v_y_flat -= shifty

        v_z_1D = np.sqrt((v_x_flat - v_x_flat[0]) ** 2 + (v_y_flat - v_y_flat[0]) ** 2)


        match s_cut.lower():
            case 'x':
                v_x_ax = v_y_flat
            case 'y':
                v_x_ax = v_x_flat
            case 'line':
                v_x_ax = v_z_1D
        #    return v_x_ax
        self.v_x_axis = v_x_ax / SI_um


class Mode:
    '''This is a base class for both EM and AC modes.'''

    def __init__(self, simres, m):
        self.sim_result = simres
        self.mode_num = m

        self.d_fields_2D = {}

        self.interpolated_2D = {FieldType.EM_E: False,
                             FieldType.EM_H: False, FieldType.AC: False}

        self.poln_fracs = []  # fx, fy, ft, fz
        self.r0 = None  # centre of mass

        self.r0_offset = (0.0, 0.0)
        self.w2 = None  # second moment width
        self._width_r0_ref = None

        self.extra_data = {}
        self.analysed = False

        self.mode_interpolator = self.sim_result.get_mode_interpolator()

        self.clear_mode_plot_data()

    def clear_mode_plot_data(self):
        for k in self.d_fields_2D.keys():
            self.d_fields_2D[k] = None
        self.interpolated_2D = {FieldType.EM_E: False,
                             FieldType.EM_H: False, FieldType.AC: False}

    def plot_mode(self, comps=(), field_type=FieldType.EM_E, ax=None, n_pts=501, decorator=None, prefix='',
                  plot_params=plotmodes.PlotParams2D()):

        if prefix: plot_params['prefix'] = prefix
        if decorator is not None: plot_params['decorator'] = decorator

        fc = FieldCode(field_type, self.is_AC())
        ft = fc.as_field_type()
        plot_params['field_type'] = ft  # awkward way to get the correct filename

        self._interpolate_mode_2D(n_pts, ft)
        self._plot_me_2D(n_pts, comps, fc, plot_params, ax)

        #self.clear_mode_plot_data()

    def plot_mode_H(self, comps, ax=None, n_pts=501, decorator=None):
        """Plot magnetic field for EM modes."""
        self.plot_mode(comps, field_type=FieldType.EM_H, ax=ax, n_pts=n_pts, decorator=decorator)



    def plot_mode_1D(self, s_cut, val1, val2=None, comps=(), field_type=FieldType.EM_E, n_pts=501, prefix='',
        plot_params=plotmodes.PlotParams2D()):

        """Make a 1D plot of the field components of this mode along a line.

        s_cut is a string determining the cut direction, with allowed values: 'x', 'y', or 'line'.

        For s_cut='x', the function is plotted along y at fixed x=val1, where val1 is a float.

        For s_cut='y', the function is plotted along x at fixed x=val1, where val1 is a float.

        For s_cut='line', the function is plotted along a straight line from val1 to val2, where the latter are float tuples (x0,y0) and (x1,y1).
        """

        if prefix: plot_params['prefix'] = prefix

        femmesh.validate_1D_cut_format(s_cut, val1, val2)
        fc = FieldCode(field_type, self.is_AC())
        ft = fc.as_field_type()

        plot_params['field_type'] = ft  # awkward way to get the correct filename

        self._interpolate_mode_1D(s_cut, val1, val2, n_pts, ft)

        if s_cut.lower() == "x":
            xlab = '$y$ (μm)'        # apparently reversed labels is correct!
        elif s_cut.lower() == "y":
            xlab = '$x$ (μm)'
        else:
            xlab = '$z$ (μm)'

        self._plot_mode_1D_gencut(comps, fc, ft, xlab, s_cut, plot_params)


    def _plot_mode_1D_gencut(self, comps, fc, ft, xlab, cut, plot_params):

        logx = plot_params.get('logx', False)
        logy = plot_params.get('logy', False)

        tags=[]
        if not len(comps):
            comps = ['x', 'y', 'z', 'a']

        tags = [FieldTag.make_from_field_and_component(ft, comp) for comp in comps]

        # add minor tag components if required
        all_tags = tags.copy()
        if not plot_params['suppress_imimre']:
            for tag in tags:
                tag.set_to_minor()
                all_tags.append(tag)

        fig, ax = plt.subplots()
        v_genx = self.mode_interpolator.v_x_axis


        for tag in all_tags:
            lc = tag.linecolor()
            ls = tag.linestyle(comps)
            co = tag.component_as_F()

            v_y = self.d_fields_1D[co]
            if logy: v_y = np.abs(v_y)
            ax.plot(v_genx, v_y, label=tag.get_tex_plot_label(), color=lc, linestyle=ls)
            if logx: ax.set_xscale('log')
            if logy: ax.set_yscale('log')

        ax.set_xlabel(xlab)
        ax.set_ylabel('fields')
        ax.legend()

        fn = plotmodes.modeplot_filename_1D(fc, plot_params, self.mode_num, cut, label='')

        plottools.save_and_close_figure(fig, fn)


    def _interpolate_mode_2D(self, n_pts, field_type):
        """Extracts fields from FEM grid to desired rectangular grid for either plotting or analysis."""

        mh = self.mode_interpolator

        mh.define_plot_grid_2D(n_pts=n_pts)

        if not self.interpolated_2D[field_type]:
            self.d_fields_2D = mh.interpolate_mode_i(self.mode_num, field_type, dims=2)
            self.interpolated_2D[field_type] = True

    def _interpolate_mode_1D(self, s_cut, val1, val2, n_pts,
                             field_type):
        mh = self.mode_interpolator

        mh.define_plot_grid_1D(s_cut, val1, val2, n_pts)

        self.d_fields_1D = mh.interpolate_mode_i(self.mode_num, field_type, dims=1)


    def _plot_me_2D(self, n_pts, comps, field_code, plot_params, ax=None):

        # TODO: weirdly, we only ax != None when there is one component to plot
        if ax is not None and len(comps) != 1:
            reporting.report_and_exit('When providing an axis to plot on, must specify exactly one modal component.')

        mh = self.mode_interpolator

        #decorator = mh.plot_params['decorator']
        decorator = plot_params['decorator']
        # All components in one plot

        decorator.set_for_multi()
        # TODO this is a kludgy way of doing this. send it through separately
        #plot_params['EM_AC'] = field_code.as_field_type()

        # can't do multiplots on a provided axis (would need a provided figure)
        if ax is None:
            plotmodes.plot_all_components(field_code, mh.xy_out, self.d_fields_2D,
                                          plot_params, self.sim_result, self.mode_num)

        # Individual component plots

        if len(comps):
            decorator.set_for_single()
            # options are ['x', 'y', 'z', 'a', 't']
            for comp in comps:  # change so this takes field type and just the x,y,z...
                cc = FieldTag.make_from_field_and_component(field_code.as_field_type(), comp)
                plotmodes.plot_one_component(field_code,
                    mh.xy_out, self.d_fields_2D, plot_params, self.mode_num, cc, ax)





    def add_mode_data(self, d):
        '''Adds a dictionary of user-defined information about a mode.

           :param dict d: Dict of (str, data) tuples of user-defined information about a mode.
      '''
        self.extra_data.update(d)

    def get_mode_data(self):
        '''Return dictionary of user-defined information about the mode.

           :return: Dictionary of user-defined information about the mode.
           :rtype: dict(str, obj)
           '''
        return self.extra_data

    def field_fracs(self):
        '''Returns tuple (*fx*, *fy*, *fz*, *ft*) of "fraction" of mode contained in *x*, *y*, *z* or *t* (sum of transverse *x+y*) components.

           Note that *fraction* is defined through a simple overlap integral. It does not necessarily represent the fraction of energy density in the component.

           :return: Tuple of mode fractions
           :rtype: tuple(float, float, float, float)
        '''
        if not self.analysed:
            reporting.report_and_exit('mode has not being analysed')

        return self.poln_fracs

    def __str__(self):
        '''String representation of the mode.'''
        s = 'Abstract mode class'
        return s

    def is_poln_ex(self):
        '''Returns true if mode is predominantly x-polarised (ie if *fx*>0.7).

           :rtype: bool
           '''
        polthresh = .7
        return self.poln_fracs[0] > polthresh

    def is_poln_ey(self):
        '''Returns true if mode is predominantly y-polarised (ie if *fy*>0.7).

           :rtype: bool
           '''
        polthresh = .7
        return self.poln_fracs[1] > polthresh

    def is_poln_indeterminate(self):
        '''Returns true if transverse polarisation is neither predominantly *x* or *y* oriented.

           :rtype: bool
           '''
        return not (self.is_poln_ex() or self.is_poln_ey())

    def is_EM(self):
        '''Returns true if the mode is an electromagnetic mode.

           :rtype: bool
           '''
        return self.sim_result.is_EM()

    def is_AC(self):
        '''Returns true if the mode is an acoustic mode.

           :rtype: bool
           '''
        return self.sim_result.is_AC()

    def center_of_mass(self):
        '''Returns the centre of mass of the mode relative to the specified origin.

           :rtype: float
        '''
        return self.r0-self.r0_offset

    def second_moment_widths(self):
        r'''Returns the second moment widths :math:`(w_x, w_y, \sqrt{w_x^2+w_y^2})` of the mode relative to the specified origin.

           :rtype: (float, float, float)
           '''
        return self.w2

    def center_of_mass_x(self):
        '''Returns the $x$ component moment of the centre of mass  of the mode.

           :rtype: float
           '''
        return self.r0[0]-self.r0_offset[0]

    def center_of_mass_y(self):
        '''Returns the $y$ component moment of the centre of mass  of the mode.

           :rtype: float
           '''
        return self.r0[1]-self.r0_offset[1]

    def wx(self):
        '''Returns the $x$ component moment of the second moment width.

           :rtype: float
           '''
        return self.w2[0]

    def wy(self):
        '''Returns the $y$ component moment of the second moment width.

           :rtype: float
           '''
        return self.w2[1]

    def w0(self):
        r'''Returns the combined second moment width :math:`\sqrt{w_x^2+w_y^2}`.

           :rtype: float
           '''
        return self.w2[2]

    def set_width_r0_reference(self, x0, y0):
        '''Set reference point for calculation of second moment width.

        Positions are measured in microns.'''

        #self._width_r0_ref=(x0/SI_um, y0/SI_um)
        self._width_r0_ref=(x0, y0)

    def set_r0_offset(self, x0, y0):
        '''Sets the transverse position in the grid that is to be regarded as the origin for calculations of center-of-mass.

           This can be useful in aligning the FEM coordinate grid with a physically sensible place in the waveguide.

           :param float x0: *x* position of nominal origin.
           :param float y0: *y* position of nominal origin.
      '''
        self.r0_offset = (x0, y0)

    def analyse_mode(self, n_pts=501, ft=FieldType.EM_E):
        '''Perform a series of measurements on the mode *f* to determine polarisation fractions, second moment widths etc.

           :param array v_x: Vector of x points.
           :param array v_y: Vector of y points.
           :param array m_Refx: Matrix of real part of fx.
           :param array m_Refy: Matrix of real part of fy.
           :param array m_Refz: Matrix of real part of fz.
           :param array m_Imfx: Matrix of imaginary part of fx.
           :param array m_Imfy: Matrix of imaginary part of fy.
           :param array m_Imfz: Matrix of imaginary part of fz.
           '''

        self._interpolate_mode_2D(n_pts, ft)

        self.analysed = True

        mh = self.mode_interpolator

        # Tranposed indexing to get image style ordering
        # These are 'output' domains so in microns
        m_x = mh.xy_out['m_x'].T
        m_y = mh.xy_out['m_y'].T

        dx = m_x[1, 0]-m_x[0, 0]
        dy = m_y[1, 1]-m_y[1, 0]

        #
        # Polarisation fractions
        #
        mFs = self.d_fields_2D  # the incoming fields, not necessarily normalised in any way

        # unit = [|F|^2], mag. \approx 1
        m_Fx2 = mFs['Fxr']**2 + mFs['Fxi']**2
        m_Fy2 = mFs['Fyr']**2 + mFs['Fyi']**2
        m_Fz2 = mFs['Fzr']**2 + mFs['Fzi']**2
        m_Fall2 = m_Fx2 + m_Fy2 + m_Fz2          # unit = [|F|^2]

        # unit = [|F|^2] um^2, mag. \approx 1
        s_fx = int2D_trapz(m_Fx2, dx=dx, dy=dy)
        s_fy = int2D_trapz(m_Fy2, dx=dx, dy=dy)
        s_fz = int2D_trapz(m_Fz2, dx=dx, dy=dy)
        s_f = s_fx+s_fy+s_fz

        f_x = s_fx/s_f                           # dimensionless
        f_y = s_fy/s_f
        f_z = s_fz/s_f
        f_t = f_x+f_y

        self.poln_fracs = [f_x, f_y, f_t, f_z]


        #
        # Positions and widths        #

        # Flipping upside down y to get sensible values for r0 position.
        m_yud = np.flipud(m_y)

        m_xmod = m_x * m_Fall2  # could do this by broadcasting without meshgrid?
        m_ymod = m_yud * m_Fall2

        x0 = int2D_trapz(m_xmod, dx, dy)/s_f          # unit = um
        y0 = int2D_trapz(m_ymod, dx, dy)/s_f

        # This just makes the outdata look cleaner on plots, by avoiding a -0.000
        if abs(x0) < 1e-6:
            x0 = 0.0
        if abs(y0) < 1e-6:
            y0 = 0.0

        # unit = um^2 [|F|^2]
        if self._width_r0_ref is None: # allow user setting of the second width moment reference point
            w_x0 = x0
            w_y0 = y0
        else:
            w_x0 = self._width_r0_ref[0] #
            w_y0 = self._width_r0_ref[1]

        m_x2mod = np.power((m_x - w_x0), 2) * m_Fall2
        m_y2mod = np.power((m_yud - w_y0), 2) * m_Fall2
        w2x = sqrt(int2D_trapz(m_x2mod, dx, dy)/s_f)
        w2y = sqrt(int2D_trapz(m_y2mod, dx, dy)/s_f)

        w2 = sqrt(w2x*w2x+w2y*w2y)
        self.r0 = np.array([x0, y0])
        self.w2 = np.array([w2x, w2y, w2])



    def write_mode(self, prefix='', n_points=501,
                   field_type=FieldType.EM_E):

        ft = FieldType.AC if self.is_AC() else FieldType(field_type)

        self._interpolate_mode_2D(n_points, ft)

        md_interp = self.mode_interpolator

        # TODO sort out wanting H fields on a singl direct call to this fnc

        v_x = md_interp.xy_out['v_x']
        v_y = md_interp.xy_out['v_y']
        s_xy= f'vx: [{v_x[0]:.6f},{v_x[-1]:.6f}: {len(v_x)}], ' + f'vy: [{v_y[0]:.6f},{v_y[-1]:.6f}: {len(v_y)}]'

        ccs = ('x', 'y', 'z')
        for cc in ccs:
            ftag = FieldTag.make_from_field_and_component(ft, cc)
            ftlab = ftag.field_type_label()
            pref=str(numbat.NumBATApp().outdir_fields_path(prefix))

            longpref=f'{pref}/{ftlab}_mode_{self.mode_num:02d}_{ftag.field_component()}'

            _write_one_component_to_file(self.mode_num,
                                         longpref, ftag, s_xy, self.d_fields_2D)

    def write_mode_1D(self, s_cut, val1, val2=None, n_points=501,
                      prefix='', field_type=FieldType.EM_E):

        self._interpolate_mode_1D(s_cut, val1, val2, n_points, field_type)

        pref=str(numbat.NumBATApp().outdir_fields_path(prefix))
        ftag = FieldTag.make_from_field(field_type)
        ftlab = ftag.field_type_label()

        longpref=f'{pref}/{ftlab}_mode_{self.mode_num:02d}_{s_cut}cut'
        fname = longpref+'.txt'

        v_genx = self.mode_interpolator.v_x_axis
        nr = len(v_genx)
        nc = 8 # v_genx, vFxr, Fxi, Fyr, Fyi, Fzr, Fzi, Fa

        flds = np.zeros([nr,nc], dtype=np.float64)
        header='x Fxr Fxi Fyr Fyi Fzr Fzi Fa'

        flds[:,0] = v_genx
        flds[:,1] = self.d_fields_1D['Fxr']
        flds[:,2] = self.d_fields_1D['Fxi']
        flds[:,3] = self.d_fields_1D['Fyr']
        flds[:,4] = self.d_fields_1D['Fyi']
        flds[:,5] = self.d_fields_1D['Fzr']
        flds[:,6] = self.d_fields_1D['Fzi']
        flds[:,7] = self.d_fields_1D['Fa']

        np.savetxt(fname, flds, header=header, fmt='%16.8f')


    def plot_mode_raw_fem(self, comps):
        '''Plot the requested field components on the sim mesh with no interpolation.'''

        simres = self.sim_result
        fem_evecs = simres.fem_evecs_for_ft(self.field_code.as_field_type())
        plotmoderaw.do_raw_fem_mode_plot(comps, self.mode_interpolator,
                                         self.sim_result.fem_mesh, fem_evecs, self.mode_num)

    def plot_strain(self):
        if not self.sim_result.is_AC():
            print("Doing strain in an EM sim.!")
        print('doing strain')
        self.mode_interpolator.plot_strain_mode_i(self.mode_num)


class ModeEM(Mode):
    '''Class representing a single electromagnetic (EM) mode.'''

    def __init__(self, sim, m):
        super().__init__(sim, m)
        self.field_code = FieldCode(FieldType.EM_E)

    def __str__(self):
        s = f'EM mode # {self.mode_num}'
        return s

    def __repr__(self):
        s = f'ModeEM(simres, mode_num={self.mode_num})'
        return s



class ModeAC(Mode):
    '''Class representing a single acoustic (AC) mode.'''

    def __init__(self, sim, m):
        super().__init__(sim, m)

        self.field_code = FieldCode(FieldType.AC)

        self.gain = {}  # { (EM_p_i, EM_s_j): gain}
        self.gain_PE = {}
        self.gain_MB = {}

    def __str__(self):
        s = 'AC mode # {self.mode_num}'
        return s

    def __repr__(self):
        s = f'ModeAC(simres, mode_num={self.mode_num})'
        return s


def _write_one_component_to_file(md_index, longpref, ftag, s_xy, d_fields):

    fc = ftag.component_as_F() # eg 'Fabs'

    if ftag.is_abs():  # a real valued quantity
        fname = longpref+'.txt'
        fld = d_fields[fc]
        header=f'Mode {md_index:02d}, {fc}, '  + s_xy
        np.savetxt(fname, fld,header=header)

    else:
        fname = longpref+'_re.txt'
        tagmaj = f'F{ftag._xyz}r'
        tagmin = f'F{ftag._xyz}i'

        fld = d_fields[tagmaj]  # eg 'Fxr'
        header=f'Mode {md_index:02d}, {fc}_re, ' + s_xy
        np.savetxt(fname, fld,header=header)

        fname = longpref+'_im.txt'
        fld = d_fields[tagmin]  # eg 'Fxr'
        header=f'Mode {md_index:02d}, {fc}_im, ' + s_xy
        np.savetxt(fname, fld,header=header)
