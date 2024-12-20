
from math import sqrt
import numpy as np

import numbat
from nbtypes import FieldType, component_t, SI_um, SI_vacuum_impedance_Z0
import reporting

from numbattools import int2d_trapz, np_min_max
import plotmodes
import plotmoderaw

# # Checks of mesh and triangles satisfy conditions for triangulation
# # Quadratic algorithm. Use on the smallest grid possible
# def check_triangulation(vx, vy, triangs):
#     # are points unique
#     print('\n\nChecking triangulation goodness')
#     npts = len(vx)
#     dsepmin = 1e6
#     dsi = 0
#     dsj = 0
#     for i in range(npts):
#         for j in range(i+1, npts):
#             dsep = sqrt((vx[i]-vx[j])**2 + (vy[i]-vy[j])**2)
#             if dsep < dsepmin:
#                 dsepmin = dsep
#                 dsi = i
#                 dsj = j

#     print('  Closest space of triangle points was', dsepmin)
#     if dsepmin < 1e-11:
#         msg = f'Point collision at {dsi}, {
#             dsj}: ({vx[dsi]},{vy[dsi]}) =  ({vx[dsj]},{vy[dsj]}).'
#         msg += '\nIt seems the mesh grid reordering has failed.'
#         reporting.report_and_exit(msg)

#     # is list of triangles unique
#     s_vtri = set()
#     clean = True
#     for tri in triangs:
#         stri = str(tri)
#         if stri in s_vtri:
#             print("        Double triangle at", stri)
#             clean = False
#         else:
#             s_vtri.add(stri)
#     if clean:
#         print("  No doubled triangles found")
#     else:
#         print("  Found doubled triangles")


class ModePlotHelper:
    '''Helper class for plotting modes.
       Factors common info from Simulation that each mode can draw on, but we only need to do once for each Sim.
       '''

    def __init__(self, simresult):  # , field_type):
        self.sim_result = simresult
        self.setup_for_npoints = 0

        self.plot_params = {}

        self.xy_raw = {}
        self.xy_out = {}

        self.interper_f = None

        self._init_plot_params()

        self.zero_arrays()

    def zero_arrays(self):
        self.interper_f = None

        self.xy_raw = {}
        self.xy_out = {}

    def cleanup(self):
        del self.interper_f

        self.xy_raw = {}
        self.xy_out = {}

        self.zero_arrays()

    def _init_plot_params(self):
        self.plot_params = {'xlim_min': 0, 'xlim_max': 0, 'ylim_min': 0, 'ylim_max': 0,
                            'aspect': 1.0,
                            'ticks': True, 'num_ticks': None,
                            'colorbar': True, 'contours': False, 'contour_lst': None,
                            'EM_AC': FieldType.EM_E,
                            'hide_vector_field': False,
                            'prefix': 'tmp', 'suffix': '',
                            'decorator': plotmodes.Decorator(),
                            'suppress_imimre': True,
                            'quiver_points': 30
                            }

    def update_plot_params(self, d_params):
        self.plot_params.update(d_params)

    def interpolate_mode_i(self, ival, field_type):
        # construct the meshed field from fortran solution

        simres = self.sim_result

        fm = simres.fem_mesh

        # extract the field data at every node of every elt for the desired mode and field

        v_Fx6p = np.zeros(6*fm.n_msh_el, dtype=np.complex128)
        v_Fy6p = np.zeros(6*fm.n_msh_el, dtype=np.complex128)
        v_Fz6p = np.zeros(6*fm.n_msh_el, dtype=np.complex128)

        fem_evecs = simres.fem_evecs_H if field_type == FieldType.EM_H else simres.fem_evecs

        i = 0
        for i_el in range(fm.n_msh_el):
            for i_node in range(6):  # TODO: make this one xyz array so we can broadcast
                v_Fx6p[i] = fem_evecs[0, i_node, ival, i_el]
                v_Fy6p[i] = fem_evecs[1, i_node, ival, i_el]
                v_Fz6p[i] = fem_evecs[2, i_node, ival, i_el]

                i += 1

        v_F6p = np.sqrt(np.abs(v_Fx6p)**2 +
                        np.abs(v_Fy6p)**2 + np.abs(v_Fz6p)**2)

        # Always need these ones.

        m_ReFx = self.interper_f(v_Fx6p.real)
        m_ReFy = self.interper_f(v_Fy6p.real)
        m_ImFz = self.interper_f(v_Fz6p.imag)
        m_AbsF = self.interper_f(v_F6p)

        # often not needed for plotting, but are used for measuring fractions. (Could fix taht?)
        m_ImFx = self.interper_f(v_Fx6p.imag)
        m_ImFy = self.interper_f(v_Fy6p.imag)
        m_ReFz = self.interper_f(v_Fz6p.real)

        d_fields = {'Fxr': m_ReFx, 'Fxi': m_ImFx, 'Fyr': m_ReFy, 'Fyi': m_ImFy,
                    'Fzr': m_ReFz, 'Fzi': m_ImFz, 'Fabs': m_AbsF}

        return d_fields

    def _choose_plot_points(self, n_pts):
        '''Picks actual data points for the maplot grid based on requested resolution.'''
        self.setup_for_npoints = n_pts

        fm = self.sim_result.fem_mesh

        x_min, x_max = np_min_max(fm.v_nd_xy[0, :])
        y_min, y_max = np_min_max(fm.v_nd_xy[1, :])

        area = abs((x_max-x_min)*(y_max-y_min))
        n_pts_x = int(n_pts*abs(x_max-x_min)/np.sqrt(area))
        n_pts_y = int(n_pts*abs(y_max-y_min)/np.sqrt(area))

        # Now use the coords user would like to think in
        shiftx, shifty = self.sim_result.get_xyshift()
        # self.shiftx, self.shifty = shiftx, shifty  # TODO: get rid of these.

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



    def setup_plot_grid(self, n_pts=501):
        '''Define interpolation plotting grids for a nominal n_pts**2 points distributed evenly amongst x and y.'''

        if self.setup_for_npoints == n_pts:
            return  # only need to repeat if the grid density changes

        fm = self.sim_result.fem_mesh

        shiftx, shifty = self._choose_plot_points(n_pts)

        v_x_flat = self.xy_raw['m_x'].flatten('F') - shiftx
        v_y_flat = self.xy_raw['m_y'].flatten('F') - shifty

        nx, ny = len(self.xy_out['v_x']), len(self.xy_out['v_y'])

        self.interper_f = fm.make_interpolator_for_grid(v_x_flat, v_y_flat, nx, ny)

class Mode:
    '''This is a base class for both EM and AC modes.'''

    def __init__(self, simres, m):
        self.mode_num = m
        self.sim_result = simres
        self.field_type = None

        self.fracs = []  # fx, fy, ft, fz
        self.r0 = None  # centre of mass
        self.w2 = None  # second moment width
        self.r0_offset = (0.0, 0.0)
        self._width_r0_ref = None
        self.extra_data = {}
        self.analysed = False
        self.interpolated = {FieldType.EM_E: False,
                             FieldType.EM_H: False, FieldType.AC: False}
        self.d_fields = {}
        self.clear_mode_plot_data()

    def get_mode_helper(self):
        return self.sim_result.get_mode_helper()

    def prepare_mode(self, n_pts, field_type):

        mh = self.get_mode_helper()

        mh.setup_plot_grid(n_pts=n_pts)

        if self.is_AC():
            self.field_type = FieldType.AC
        else:
            self.field_type = field_type

        if not self.interpolated[field_type]:
            self.interpolate_mode(mh)
            self.interpolated[field_type] = True

    def plot_mode(self, comps, field_type=FieldType.EM_E, ax=None,
                  n_pts=501, decorator=None):  # TODO get this random parameters hooked better into mode_helper.plot_params

        self.prepare_mode(n_pts, field_type)

        mh = self.get_mode_helper()

        # FIX ME
        if decorator is not None:
            mh.plot_params['decorator'] = decorator
        elif mh.plot_params['decorator'] is None:
            # don't want to do this.
            mh.plot_params['decorator'] = plotmodes.Decorator()

        self._plot_me(mh, comps, field_type, ax)

        #self.clear_mode_plot_data()

    def plot_mode_H(self, comps):  # plot magnetic field for EM modes
        self.plot_mode(comps, field_type=FieldType.EM_H)

    def plot_mode_raw_fem(self, comps):
        '''Plot the requested field components on the sim mesh with no interpolation.'''

        simres = self.sim_result
        ft = self.field_type
        mh = self.get_mode_helper()

        fem_evecs = simres.fem_evecs_H if ft == FieldType.EM_H else simres.fem_evecs
        plotmoderaw.do_raw_fem_mode_plot(comps, mh,
                                         self.sim_result.fem_mesh, fem_evecs, self.mode_num)

    def _write_one_component_to_file(self, longpref, comp, s_xy, d_fields):
        if comp.is_abs():  # a real valued quantity
            fname = longpref+'.txt'
            fld = d_fields[comp._f_code]  # eg 'Fabs'
            header=f'# {comp._f_code}, '  + s_xy
            np.savetxt(fname, fld,header=header)

        else:
            fname = longpref+'_re.txt'
            fld = d_fields['F'+comp._xyz+'r']  # eg 'Fxr'
            print('the field', fld, d_fields['Fxr'])
            header=f'# {comp._f_code}_re, '  + s_xy
            np.savetxt(fname, fld,header=header)

            fname = longpref+'_im.txt'
            fld = d_fields['F'+comp._xyz+'i']  # eg 'Fxr'
            header=f'# {comp._f_code}__im, '  + s_xy
            np.savetxt(fname, fld,header=header)



    def write_to_file(self):
        #comps = ('exr', 'exi','eyr', 'eyi','ezr', 'ezi',
        #         'hxr', 'ehi','ehr', 'ehi','ehr', 'ehi')
        ccs = ('Fx', 'Fy', 'Fz')
        ft=self.field_type
        mh = self.get_mode_helper()
        #print('mhprops', self.analysed, self.d_fields)
        v_x = mh.xy_out['v_x']
        v_y = mh.xy_out['v_y']
        s_xy= f'v_x: {v_x[0]:.8f}, {v_x[-1]:.8f}, {len(v_x)}, ' + f'v_y: {v_y[0]:.8f}, {v_y[-1]:.8f}, {len(v_y)}'

        for cc in ccs:
            comp = component_t.make_comp_noreim(ft, cc)
            pref=numbat.NumBATApp().outpath_fields()
            longpref=f'{pref}/{comp.emac()}_mode_{self.mode_num:02d}_{comp._user_code}'
            #print('writing data file', ft,longpref, self.d_fields.keys())
            self._write_one_component_to_file(longpref, comp, s_xy, self.d_fields)

    def plot_strain(self):
        if not self.sim_result.is_AC():
            print("Doing strain in an EM sim.!")
        print('doing strain')
        mh = self.get_mode_helper()
        mh.plot_strain_mode_i(self.mode_num)

    def clear_mode_plot_data(self):
        for k in self.d_fields.keys():
            self.d_fields[k] = None
        self.interpolated = {FieldType.EM_E: False,
                             FieldType.EM_H: False, FieldType.AC: False}

    def interpolate_mode(self, mode_helper):
        mh = mode_helper
        self.d_fields = mh.interpolate_mode_i(self.mode_num, self.field_type)

        if self.field_type == FieldType.EM_H:  # scale H fields by Z0 to get common units and amplitude with E
            for m_F in self.d_fields.values():   # Do this when they are first made
                m_F *= SI_vacuum_impedance_Z0

    def _plot_me(self, mode_helper, comps, field_type, ax=None):

        # TODO: weirdly, we only ax != None when there is one component to plot
        if ax is not None and len(comps) != 1:
            print(
                '\nError: when providing an axis to plot on, must specify exactly one modal component.')
            return

        mh = mode_helper
        decorator = mh.plot_params['decorator']

        decorator.set_for_multi()
        # TODO this is a kludgy way of doing this. send it through separately
        mh.plot_params['EM_AC'] = field_type

        # can't do multiplots on a provided axis (would need a provided figure)
        if ax is None:
            plotmodes.plot_all_components(mh.xy_out, self.d_fields,
                                          mh.plot_params, self.sim_result, self.mode_num)

        if len(comps):
            decorator.set_for_single()
            # options are ['x', 'y', 'z', 'abs', 't']
            for comp in comps:  # change so this takes field type and just the x,y,z...
                cc = component_t.make_comp_from_component(field_type, comp)
                plotmodes.plot_one_component(
                    mh.xy_out, self.d_fields, mh.plot_params, self.mode_num, cc, ax)

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

        return self.fracs

    def __str__(self):
        '''String representation of the mode.'''
        s = 'Abstract mode class'
        return s

    def is_poln_ex(self):
        '''Returns true if mode is predominantly x-polarised (ie if *fx*>0.7).

           :rtype: bool
           '''
        polthresh = .7
        return self.fracs[0] > polthresh

    def is_poln_ey(self):
        '''Returns true if mode is predominantly y-polarised (ie if *fy*>0.7).

           :rtype: bool
           '''
        polthresh = .7
        return self.fracs[1] > polthresh

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

    def analyse_mode(self, n_pts=501, EM_field=FieldType.EM_E):
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

        self.prepare_mode(n_pts, EM_field)

        self.analysed = True

        mh = self.get_mode_helper()

        # Tranposed indexing to get image style ordering
        # These are 'output' domains so in microns
        m_x = mh.xy_out['m_x'].T
        m_y = mh.xy_out['m_y'].T

        dx = m_x[1, 0]-m_x[0, 0]
        dy = m_y[1, 1]-m_y[1, 0]

        mFs = self.d_fields  # the incoming fields, not necessarily normalised in any way

        # unit = [|F|^2], mag. \approx 1
        m_Fx2 = mFs['Fxr']**2 + mFs['Fxi']**2
        m_Fy2 = mFs['Fyr']**2 + mFs['Fyi']**2
        m_Fz2 = mFs['Fzr']**2 + mFs['Fzi']**2
        m_Fall2 = m_Fx2 + m_Fy2 + m_Fz2          # unit = [|F|^2]

        # unit = [|F|^2] um^2, mag. \approx 1
        s_fx = int2d_trapz(m_Fx2, dx=dx, dy=dy)
        s_fy = int2d_trapz(m_Fy2, dx=dx, dy=dy)
        s_fz = int2d_trapz(m_Fz2, dx=dx, dy=dy)
        s_f = s_fx+s_fy+s_fz

        f_x = s_fx/s_f                           # dimensionless
        f_y = s_fy/s_f
        f_z = s_fz/s_f
        f_t = f_x+f_y

        self.fracs = [f_x, f_y, f_t, f_z]

        # Flipping upside down y to get sensible values for r0 position.
        m_yud = np.flipud(m_y)

        m_xmod = m_x * m_Fall2  # could do this by broadcasting without meshgrid?
        m_ymod = m_yud * m_Fall2

        x0 = int2d_trapz(m_xmod, dx, dy)/s_f          # unit = um
        y0 = int2d_trapz(m_ymod, dx, dy)/s_f

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
        w2x = sqrt(int2d_trapz(m_x2mod, dx, dy)/s_f)
        w2y = sqrt(int2d_trapz(m_y2mod, dx, dy)/s_f)

        w2 = sqrt(w2x*w2x+w2y*w2y)
        self.r0 = np.array([x0, y0])
        self.w2 = np.array([w2x, w2y, w2])


class ModeEM(Mode):
    '''Class representing a single electromagnetic (EM) mode.'''

    def __init__(self, sim, m):
        super().__init__(sim, m)

    def __str__(self):
        s = f'EM mode # {self.mode_num}'
        return s


class ModeAC(Mode):
    '''Class representing a single acoustic (AC) mode.'''

    def __init__(self, sim, m):
        super().__init__(sim, m)

        self.gain = {}  # { (EM_p_i, EM_s_j): gain}
        self.gain_PE = {}
        self.gain_MB = {}

    def __str__(self):
        s = 'AC mode # {self.mode_num}'
        return s
