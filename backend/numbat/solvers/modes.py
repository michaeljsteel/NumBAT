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
"""Mode-domain state, analysis, and public compatibility façades.

Field sampling, rendering, and field-file serialization are delegated to
dedicated solver, plotting, and I/O modules. The compatibility helper imports
in this module preserve established public `Mode` and `ModeInterpolator`
entry points.
"""
from math import sqrt

import numpy as np

from numbat.field_selection import (
    complex_storage_keys,
    resolve_field_type,
    storage_key,
)
from numbat.io import mode_fields, strain_fields
from numbat.mesh import fem as femmesh
from numbat.nbtypes import (
    ComplexPart,
    FieldComponent,
    FieldSelection,
    FieldType,
)
from numbat.plotting import mode_rendering, strain_rendering
from numbat.plotting import modes as plotmodes
from numbat.plotting import raw_modes as plotmoderaw
from numbat.solvers import strain_sampling
from numbat.solvers.mode_sampling import ModeInterpolator  # noqa: F401
from numbat.strain_selection import (
    complex_storage_keys as strain_complex_storage_keys,
)
from numbat.strain_selection import (
    component_index,
    parse_strain_components,
)
from numbat.utils.numeric import int2D_trapz

from .. import reporting


def _field_storage_key(field_type, component, complex_part=None):
    return storage_key(FieldSelection(field_type, component, complex_part))


class Mode:
    """Base domain object for electromagnetic and acoustic modes.

    It owns mode state and analysis while delegating sampled field production,
    rendering, and text export to dedicated modules.
    """

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

        self.clear_interpolation_data()

    def clear_interpolation_data(self):
        for k in self.d_fields_2D:
            self.d_fields_2D[k] = None
        self.interpolated_2D = {FieldType.EM_E: False,
                             FieldType.EM_H: False, FieldType.AC: False}

    def plot_mode(self, comps=(), field_type=FieldType.EM_E, ax=None, n_pts=501, decorator=None, prefix='', plot_params=None, trim=None):
        """Plot selected components using a compact string or `FieldComponent` values.

        With no components, write the combined field overview. Use a compact
        string such as ``comps="exyat"`` or an iterable such as
        ``comps=(FieldComponent.EVERYTHING, FieldComponent.X)`` to request
        the overview and individual component plots. Without `EVERYTHING`,
        only individual component plots are written. Supported individual
        components are `X`/`"x"`, `Y`/`"y"`, `Z`/`"z"`, `TRANSVERSE`/`"t"`,
        and `MAGNITUDE`/`"a"`.

        `trim=(left, bottom, right, top)` trims the displayed plot bounds by
        fractional amounts and takes precedence over legacy limit keys in
        `plot_params`.
        """

        include_overview, comps = mode_rendering.split_plot_components(comps)

        if plot_params is None:
            plot_params = plotmodes.PlotParams2D()

        if prefix: plot_params['prefix'] = prefix
        if decorator is not None: plot_params['decorator'] = decorator
        if trim is not None: plot_params["trim"] = trim

        ft = resolve_field_type(field_type, simulation_is_acoustic=self.is_AC())
        plot_params['field_type'] = ft  # awkward way to get the correct filename

        self._interpolate_mode_2D(n_pts, ft)
        mode_rendering.render_mode_2d(
            field_type=ft,
            coordinates=self.mode_interpolator.xy_out,
            fields=self.d_fields_2D,
            include_overview=include_overview,
            components=comps,
            simulation_result=self.sim_result,
            mode_index=self.mode_num,
            plot_params=plot_params,
            axis=ax,
        )

        #self.clear_mode_plot_data()

    def plot_mode_H(self, comps=(), ax=None, n_pts=501, decorator=None):
        """Plot magnetic field for EM modes."""
        self.plot_mode(comps, field_type=FieldType.EM_H, ax=ax, n_pts=n_pts, decorator=decorator)



    def plot_mode_1D(self, s_cut, val1, val2=None, comps=(), field_type=FieldType.EM_E, n_pts=501, prefix='', plot_params=None):

        """Make a 1D plot of the field components of this mode along a line.

        s_cut is a string determining the cut direction, with allowed values: 'x', 'y', or 'line'.

        Field components may be a compact string such as ``comps="xyat"`` or
        an iterable of `FieldComponent` values. Supported values are
        `X`/`"x"`, `Y`/`"y"`, `Z`/`"z"`, `TRANSVERSE`/`"t"`, and
        `MAGNITUDE`/`"a"`. `EVERYTHING`/`"e"` is available only for 2D mode
        plotting and is rejected here.

        For s_cut='x', the function is plotted along y at fixed x=val1, where val1 is a float.

        For s_cut='y', the function is plotted along x at fixed x=val1, where val1 is a float.

        For s_cut='line', the function is plotted along a straight line from val1 to val2, where the latter are float tuples (x0,y0) and (x1,y1).
        """

        if plot_params is None:
            plot_params = plotmodes.PlotParams2D()

        comps = mode_rendering.components_for_1d_plot(comps)

        if prefix: plot_params['prefix'] = prefix

        femmesh.validate_1D_cut_format(s_cut, val1, val2)
        ft = resolve_field_type(field_type, simulation_is_acoustic=self.is_AC())

        plot_params['field_type'] = ft  # awkward way to get the correct filename

        self._interpolate_mode_1D(s_cut, val1, val2, n_pts, ft)

        if s_cut.lower() == "x":
            xlab = '$y$ (μm)'        # apparently reversed labels is correct!
        elif s_cut.lower() == "y":
            xlab = '$x$ (μm)'
        else:
            xlab = '$z$ (μm)'

        mode_rendering.render_mode_1d(
            field_type=ft,
            fields=self.d_fields_1D,
            x_axis=self.mode_interpolator.v_x_axis,
            components=comps,
            x_label=xlab,
            cut=s_cut,
            mode_index=self.mode_num,
            plot_params=plot_params,
        )


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
        '''Returns tuple (*fx*, *fy*, *ft*, *fz*) of "fraction" of mode contained in *x*, *y*, *t* (sum of transverse *x+y*) or *z* components.

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

    def set_width_r0_reference(self, wid_r0_ref=None):
        '''Set reference point for calculation of second moment width.

        Positions are measured in microns.'''

        #self._width_r0_ref=(x0/SI_um, y0/SI_um)
        if wid_r0_ref is not None:
            self._width_r0_ref=wid_r0_ref

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

        ft = resolve_field_type(ft, simulation_is_acoustic=self.is_AC())
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
        m_Fx2 = (
            mFs[_field_storage_key(ft, FieldComponent.X, ComplexPart.PRIMARY)]**2
            + mFs[_field_storage_key(ft, FieldComponent.X, ComplexPart.SECONDARY)]**2
        )
        m_Fy2 = (
            mFs[_field_storage_key(ft, FieldComponent.Y, ComplexPart.PRIMARY)]**2
            + mFs[_field_storage_key(ft, FieldComponent.Y, ComplexPart.SECONDARY)]**2
        )
        m_Fz2 = (
            mFs[_field_storage_key(ft, FieldComponent.Z, ComplexPart.PRIMARY)]**2
            + mFs[_field_storage_key(ft, FieldComponent.Z, ComplexPart.SECONDARY)]**2
        )
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

        # This just makes the outdata look cleaner on plots, by avoiding a -0.000.
        # But is changing actual values of r0 !
       # if abs(x0) < 1e-6:
       #     x0 = 0.0
       # if abs(y0) < 1e-6:
       #     y0 = 0.0

        # unit = um^2 [|F|^2]
        if self._width_r0_ref is None: # allow user setting of the second width moment reference point
            w_x0 = x0
            w_y0 = y0
        else:
            w_x0 = self._width_r0_ref[0] #
            w_y0 = self._width_r0_ref[1]

        m_x2mod = np.power((m_x - w_x0), 2) * m_Fall2
        m_y2mod = np.power((m_yud - w_y0), 2) * m_Fall2


        #m_y2modb = np.power((m_yud -y0), 2) * m_Fall2
        #print('cols', m_yud[:,0])
        #print('rows', m_yud[0,:])
        #print('colsb', m_yud[:,0]-w_y0)
        #print('\nrowsa',
        #      np.transpose(np.array([
        #          m_yud[0,:],
        #      np.sum(m_y2mod, axis=0)/s_f*dx*dy,
        #      np.sum(m_y2modb, axis=0)/s_f*dx*dy
        #      ])
        #      ))
        #sys.exit(1)

        w2x = sqrt(int2D_trapz(m_x2mod, dx, dy)/s_f)
        w2y = sqrt(int2D_trapz(m_y2mod, dx, dy)/s_f)

        w2 = sqrt(w2x*w2x+w2y*w2y)
        self.r0 = np.array([x0, y0])
        self.w2 = np.array([w2x, w2y, w2])

        #print('rowsb', self.mode_num, x0, y0, w_x0, w_y0, w2y)# m_yud[0,:]-w_y0)



    def write_mode(self, prefix='', n_points=501,
                   field_type=FieldType.EM_E):

        ft = resolve_field_type(field_type, simulation_is_acoustic=self.is_AC())

        self._interpolate_mode_2D(n_points, ft)
        mode_fields.write_mode_fields_2d(
            prefix=prefix,
            field_type=ft,
            mode_index=self.mode_num,
            coordinates=self.mode_interpolator.xy_out,
            fields=self.d_fields_2D,
        )

    def write_mode_1D(self, s_cut, val1, val2=None, n_points=501,
                      prefix='', field_type=FieldType.EM_E):

        field_type = resolve_field_type(
            field_type, simulation_is_acoustic=self.is_AC()
        )
        self._interpolate_mode_1D(s_cut, val1, val2, n_points, field_type)
        mode_fields.write_mode_fields_1d(
            prefix=prefix,
            field_type=field_type,
            mode_index=self.mode_num,
            cut=s_cut,
            x_axis=self.mode_interpolator.v_x_axis,
            fields=self.d_fields_1D,
        )

    def get_mode_1D(self, s_cut, val1, val2=None, n_points=501, field_type=FieldType.EM_E):

        field_type = resolve_field_type(
            field_type, simulation_is_acoustic=self.is_AC()
        )
        self._interpolate_mode_1D(s_cut, val1, val2, n_points, field_type)

        v_genx = self.mode_interpolator.v_x_axis
        nr = len(v_genx)
        nc = 4 # vFx, Fy, Fz, Fa

        flds = np.zeros([nr,nc], dtype=np.complex128)
        #header='x Fxr Fxi Fyr Fyi Fzr Fzi Fa'

        for index, component in enumerate(
            (FieldComponent.X, FieldComponent.Y, FieldComponent.Z)
        ):
            real_key, imaginary_key = complex_storage_keys(field_type, component)
            flds[:,index] = self.d_fields_1D[real_key] + 1j * self.d_fields_1D[imaginary_key]
        flds[:,3] = self.d_fields_1D[_field_storage_key(field_type, FieldComponent.MAGNITUDE)]

        return v_genx, flds


    def plot_mode_raw_fem(self, comps):
        '''Plot the requested field components on the sim mesh with no interpolation.'''

        simres = self.sim_result
        fem_evecs = simres.fem_evecs_for_ft(self.field_type)
        plotmoderaw.do_raw_fem_mode_plot(
            comps, self.sim_result.fem_mesh, fem_evecs, self.mode_num
        )

    def plot_strain(self):
        if not self.sim_result.is_AC():
            print("Doing strain in an EM sim.!")
        print('doing strain')
        self.mode_interpolator.plot_strain_mode_i(self.mode_num)


class ModeEM(Mode):
    '''Class representing a single electromagnetic (EM) mode.'''

    def __init__(self, sim, m):
        super().__init__(sim, m)
        self.field_type = FieldType.EM_E

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

        self.field_type = FieldType.AC

        self.gain = {}  # { (EM_p_i, EM_s_j): gain}
        self.gain_PE = {}
        self.gain_MB = {}

    def __str__(self):
        s = 'AC mode # {self.mode_num}'
        return s

    def __repr__(self):
        s = f'ModeAC(simres, mode_num={self.mode_num})'
        return s

    def plot_strain(
        self, comps=(), ax=None, n_pts=501, prefix="", plot_params=None, trim=None
    ):
        """Plot native-FEM acoustic strain on the standard 2D sampling grid.

        `trim=(left, bottom, right, top)` trims the displayed plot bounds by
        fractional amounts.
        """
        if plot_params is None:
            plot_params = plotmodes.PlotParams2D()
        if prefix:
            plot_params["prefix"] = prefix
        if trim is not None:
            plot_params["trim"] = trim

        components = parse_strain_components(comps)
        coordinates, fields = self._sample_strain_2d(n_pts)
        strain_rendering.render_strain_2d(
            coordinates=coordinates,
            fields=fields,
            components=components,
            simulation_result=self.sim_result,
            mode_index=self.mode_num,
            plot_params=plot_params,
            axis=ax,
        )

    def write_strain(self, prefix="", n_pts=501):
        """Write sampled native-FEM acoustic strain as twelve 2D text files."""
        coordinates, fields = self._sample_strain_2d(n_pts)
        strain_fields.write_strain_fields_2d(
            prefix=prefix,
            mode_index=self.mode_num,
            coordinates=coordinates,
            fields=fields,
        )

    def plot_strain_1D(
        self, s_cut, val1, val2=None, comps=(), n_pts=501, prefix="", plot_params=None
    ):
        """Plot display-primary strain components along an acoustic mode cut."""
        if plot_params is None:
            plot_params = plotmodes.PlotParams1D()
        if prefix:
            plot_params["prefix"] = prefix

        components = parse_strain_components(comps)
        x_axis, fields = self._sample_strain_1d(s_cut, val1, val2, n_pts)
        x_label = {
            "x": "$y$ (μm)",
            "y": "$x$ (μm)",
            "line": "$z$ (μm)",
        }[s_cut.lower()]
        strain_rendering.render_strain_1d(
            x_axis=x_axis,
            fields=fields,
            components=components,
            x_label=x_label,
            cut=s_cut,
            mode_index=self.mode_num,
            plot_params=plot_params,
        )

    def write_strain_1D(self, s_cut, val1, val2=None, n_pts=501, prefix=""):
        """Write all physical strain parts along an acoustic mode cut."""
        x_axis, fields = self._sample_strain_1d(s_cut, val1, val2, n_pts)
        strain_fields.write_strain_fields_1d(
            prefix=prefix,
            mode_index=self.mode_num,
            cut=s_cut,
            x_axis=x_axis,
            fields=fields,
        )

    def _sample_strain_2d(self, n_pts):
        """Sample all physical strain parts onto the standard 2D plot grid."""
        mode_interpolator = self.mode_interpolator
        mode_interpolator.define_plot_grid_2D(n_pts=n_pts)
        shift_x, shift_y = self.sim_result.get_xyshift()
        sample_points = np.column_stack((
            mode_interpolator.xy_raw["m_x"].flatten("F") - shift_x,
            mode_interpolator.xy_raw["m_y"].flatten("F") - shift_y,
        ))
        native_samples = strain_sampling.calculate_mode_strain(
            self.sim_result, self.mode_num
        )
        sampled_values = strain_sampling.make_native_strain_interpolator(
            native_samples
        )(sample_points)

        nx = len(mode_interpolator.xy_out["v_x"])
        ny = len(mode_interpolator.xy_out["v_y"])
        values = sampled_values.reshape(nx, ny, len(strain_sampling.StrainComponent))
        return mode_interpolator.xy_out, self._strain_fields_from_values(values)

    def _sample_strain_1d(self, s_cut, val1, val2, n_pts):
        """Sample all physical strain parts along a native one-sided cut."""
        x_axis, values = strain_sampling.sample_mode_strain_cut(
            self.sim_result,
            self.mode_num,
            s_cut,
            val1,
            val2,
            n_points=n_pts,
        )
        return x_axis, self._strain_fields_from_values(values)

    @staticmethod
    def _strain_fields_from_values(values):
        fields = {}
        for component in strain_sampling.StrainComponent:
            real_key, imaginary_key = strain_complex_storage_keys(component)
            component_values = values[..., component_index(component)]
            fields[real_key] = component_values.real
            fields[imaginary_key] = component_values.imag
        return fields
