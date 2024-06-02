
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri

import numbat
from nbtypes import FieldType, component_t, SI_um, vacuum_impedance_Z0
from numbattools import int2d, np_min_max, save_and_close_figure

import plotmodes
import reporting


# Checks of mesh and triangles satisfy conditions for triangulation
# Quadratic algorithm. Use on the smallest grid possible
def check_triangulation(vx, vy, triangs):
    # are points unique
    print('\n\nChecking triangulation goodness')
    npts = len(vx)
    dsepmin = 1e6
    dsi=0
    dsj=0
    for i in range(npts):
        for j in range(i+1, npts):
            dsep = sqrt( (vx[i]-vx[j])**2 +(vy[i]-vy[j])**2)
            if dsep < dsepmin:
                dsepmin = dsep
                dsi=i
                dsj=j


    print('  Closest space of triangle points was', dsepmin)
    if dsepmin < 1e-11:
        msg=f'Point collision at {dsi}, {dsj}: ({vx[dsi]},{vy[dsi]}) =  ({vx[dsj]},{vy[dsj]}).'
        msg+='\nIt seems the mesh grid reordering has failed.'
        reporting.report_and_exit(msg)

    # is list of triangles unique
    s_vtri = set()
    clean = True
    for tri in triangs:
        stri = str(tri)
        if stri in s_vtri:
            print("        Double triangle at", stri)
            clean = False
        else:
            s_vtri.add(stri)
    if clean:
        print("  No doubled triangles found")
    else:
        print("  Found doubled triangles")



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
        #self.triang6p = None  #needed after triangularions are built?
        #self.triang1p = None
        self.interper_f = None

        #self.v_x6p = None
        #self.v_y6p = None
        #self.v_Fx6p = None
        #self.v_Fy6p = None
        #self.v_Fz6p = None


        #self.v_x = None  # final plot coordinate arrays in microns
        #self.v_y = None
        #self.m_X = None  # mesh grid versions of the same
        #self.m_Y = None

        self.xy_raW = {}

        self.xy_out = {}

    def cleanup(self):
        # Now that these are part of the sim object,
        # we need to get rid of them to allow pickling of sim files, as they contain C++ objects

        #del self.triang6p
        #del self.triang1p
        del self.interper_f

        #del self.v_x, self.v_y, self.m_X, self.m_Y
        #del self.v_x6p, self.v_y6p
        #del self.v_Fx6p, self.v_Fy6p, self.v_Fz6p

        self.xy_raw = {}
        self.xy_out = {}


        self.zero_arrays()


#    def make_picklable(self): #before this object can be copied or saved it needs internal C++ objects removed. They are rebuilt when needed.
#        if not self.triang6p is None: del self.triang6p._cpp_triangulation
#        if not self.triang1p is None: del self.triang1p._cpp_triangulation

    def _init_plot_params(self):
        self.plot_params = {'xlim_min': 0, 'xlim_max': 0, 'ylim_min': 0, 'ylim_max': 0,
                            'ticks': True, 'num_ticks': None,
                            'colorbar': True, 'contours': False, 'contour_lst': None,
                            'EM_AC': FieldType.EM_E,
                            'prefix': 'tmp', 'suffix': '',
                            'decorator': plotmodes.Decorator(),
                            'suppress_imimre': True,
                            'quiver_points': 30
                            }

    def update_plot_params(self, d_params):
        self.plot_params.update(d_params)

    # def set_plot_params(self,
    #                     xlim_min=0, xlim_max=0, ylim_min=0, ylim_max=0,
    #                     field_type=FieldType.EM_E,
    #                     quiver_points=30,
    #                     num_ticks=None, ticks=False, colorbar=True, contours=False, contour_lst=None,
    #                     suppress_imimre=True,
    #                     prefix='tmp', suffix='', decorator=plotting.Decorator(), ):

    #     #pf = numbat.NumBATApp().path_fields()
    #     #if prefix and not Path(pf).exists(): Path(pf).mkdir()  # TODO: shouldn't ned Path() wrapper

    #     #field_type = FieldType.AC if self.sim_result.is_AC() else FieldType.from_str(field_type)

    #     self.plot_params = {'xlim_min': xlim_min, 'xlim_max': xlim_max, 'ylim_min': ylim_min,
    #                         'ylim_max': ylim_max, 'ticks': ticks, 'num_ticks': num_ticks,
    #                         'colorbar': colorbar, 'contours': contours, 'contour_lst': contour_lst, 'EM_AC': field_type,
    #                         'prefix': prefix, 'suffix': suffix,
    #                         'decorator': decorator,
    #                         'suppress_imimre': suppress_imimre,
    #                         'quiver_points': quiver_points
    #                         }


    def interpolate_mode_i(self, ival, field_type):
        # self.v_Fx6p etc could propbably be made local to this function
        # construct the meshed field from fortran solution

        simres = self.sim_result

        fm = simres.fem_mesh

        # extract the field data at every node of every elt for the desired mode and field
        # v_Fxxx have length 6*n_msh_el

        v_Fx6p = np.zeros(6*fm.n_msh_el, dtype=np.complex128)
        v_Fy6p = np.zeros(6*fm.n_msh_el, dtype=np.complex128)
        v_Fz6p = np.zeros(6*fm.n_msh_el, dtype=np.complex128)

        fem_evecs = simres.fem_evecs_H if field_type == FieldType.EM_H else simres.fem_evecs

        i = 0
        for i_el in range(fm.n_msh_el):
            for i_node in range(6):  # TODO: make this one xyz array so we can broadcast
                #if field_type in (FieldType.EM_E, FieldType.AC):
                    v_Fx6p[i] = simres.fem_evecs[0, i_node, ival, i_el]
                    v_Fy6p[i] = simres.fem_evecs[1, i_node, ival, i_el]
                    v_Fz6p[i] = simres.fem_evecs[2, i_node, ival, i_el]
                #else:
                #    v_Fx6p[i] = simres.fem_evecs_H[0, i_node, ival, i_el]
                #    v_Fy6p[i] = simres.fem_evecs_H[1, i_node, ival, i_el]
                #    v_Fz6p[i] = simres.fem_evecs_H[2, i_node, ival, i_el]

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

    def _choose_plot_points(self, n_points):
        '''Picks actual data points for the plot grid based on requested resolution.'''
        self.setup_for_npoints = n_points

        fm = self.sim_result.fem_mesh

        x_min, x_max = np_min_max(fm.mesh_xy[0,:])
        y_min, y_max = np_min_max(fm.mesh_xy[1,:])

        area = abs((x_max-x_min)*(y_max-y_min))
        self.n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
        self.n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))

        # Now use the coords user would like to think in
        shiftx, shifty = self.sim_result.get_xyshift()
        #self.shiftx, self.shifty = shiftx, shifty  # TODO: get rid of these.

        # These are the actual x and y domains of the final plots
        v_x = np.linspace(x_min, x_max, self.n_pts_x)
        v_y = np.linspace(y_min, y_max, self.n_pts_y)
        m_X, m_Y = np.meshgrid(v_x, v_y)

        self.xy_raw = {'v_x': v_x, 'v_y': v_y, 'm_x': m_X, 'm_y': m_Y}

        v_x_out = (v_x + shiftx) / SI_um
        v_y_out = (v_y + shifty) / SI_um
        m_X_out, m_Y_out = np.meshgrid(v_x_out, v_y_out)

        self.xy_out = {'v_x': v_x_out, 'v_y': v_y_out, 'm_x': m_X_out, 'm_y': m_Y_out}

        print('''
  Structure has raw domain(x,y)   = [{0:.5f}, {1:.5f}] x [ {2:.5f}, {3:.5f}] (um),
                mapped to (x',y') = [{4:.5f}, {5:.5f}] x [ {6:.5f}, {7:.5f}] (um)
                    '''.format(
            v_x[0]/SI_um, v_x[-1]/SI_um, v_y[0]/SI_um, v_y[-1]/SI_um,
            v_x_out[0], v_x_out[-1], v_y_out[0], v_y_out[-1]))
        return shiftx, shifty

    def _save_triangulation_plots(self, triang1p, triang6p):
        fig, axs=plt.subplots(1,2)
        axs[0].triplot(triang1p, linewidth=.5)
        axs[1].triplot(triang6p, linewidth=.5)
        for ax in axs:
            ax.set_aspect(1.0)
            ax.scatter(self.mesh_xy[:,0],self.mesh_xy[:,1], s=2, c='red')

        pref = numbat.NumBATApp().outprefix()
        fname = pref + f'-{'ac' if self.sim_result.is_AC else 'em'}_triplots.png'
        save_and_close_figure(fig, fname)



    def setup_plot_grid(self, n_points=501):
        '''Define interpolation plotting grids for a nominal n_points**2 points distributed evenly amongst x and y.'''

        if self.setup_for_npoints == n_points:
            return  # only need to repeat if the grid density changes

        fm = self.sim_result.fem_mesh

        shiftx, shifty = self._choose_plot_points(n_points)

        # unrolling data for the interpolators
        # TODO: for EM, table_nod seems to be identical to the MailData one
        #       mesh_xy seems to be the same but with some fractional scaling.

        # Sim version is in fortran ordering
        # This version is in python ordering.  Eeek!
        self.table_nod = fm.table_nod.T
        self.mesh_xy = fm.mesh_xy.T  # is v_x, v_y  * d_in_m

        # dense triangulation with multiple points
        v_x6p = np.zeros(6*fm.n_msh_el)
        v_y6p = np.zeros(6*fm.n_msh_el)
        #self.v_Fx6p = np.zeros(6*sim.n_msh_el, dtype=np.complex128)
        #self.v_Fy6p = np.zeros(6*sim.n_msh_el, dtype=np.complex128)
        #self.v_Fz6p = np.zeros(6*sim.n_msh_el, dtype=np.complex128)
        v_triang6p = []

        # In table_nod
        # Nodes around a triangle element are numbered as corners: 0 1 2,  midpts: 3,4,5
        # This induces a 4-triangle sub-triangulation of each element, with clockwise vertices
        # (0 3 5), (1, 4, 3), (2, 5, 4),  (3, 4, 5)


        # create sub-triangles from combos of the element nodes
        for idx in range(0, 6*fm.n_msh_el, 6):
            triangles = [[idx+0, idx+3, idx+5],
                         [idx+1, idx+4, idx+3],
                         [idx+2, idx+5, idx+4],
                         [idx+3, idx+4, idx+5]]
            v_triang6p.extend(triangles)


        tabnod_py = self.table_nod -1  #  shift fortran to python indexing

        # Create vectors v_x6p, v_y6p which are unwrapped points at nodes of each element
        # i is the index for the coordinates FIND A BETTER NAME
        i = 0
        for i_el in range(fm.n_msh_el):
            for i_node in range(6):
                i_ex = tabnod_py[i_el, i_node]
                v_x6p[i] = self.mesh_xy[i_ex, 0]
                v_y6p[i] = self.mesh_xy[i_ex, 1]
                i += 1



        # Interpolate onto triangular grid - honest to FEM elements
        # dense triangulation with unique points
        v_triang1p = []
        #table_nod = self.table_nod
        for i_el in np.arange(fm.n_msh_el):
            triangles = [[tabnod_py[i_el, 0], tabnod_py[i_el, 3], tabnod_py[i_el, 5]],
                         [tabnod_py[i_el, 1], tabnod_py[i_el, 4], tabnod_py[i_el, 3]],
                         [tabnod_py[i_el, 2], tabnod_py[i_el, 5], tabnod_py[i_el, 4]],
                         [tabnod_py[i_el, 3], tabnod_py[i_el, 4], tabnod_py[i_el, 5]]]
            v_triang1p.extend(triangles)

        # TODO: There seems to be no difference between v_triang6p and v_triang1p.

        # This is for testing only. Normally turn off
        check_tris = False
        if check_tris:
            check_triangulation( self.mesh_xy[:,0], self.mesh_xy[:,1], self.v_triang1p)

        # triangulations:  x and y coords of all points, list of triangles defined by triples of indices of the points
        tri_triang6p = matplotlib.tri.Triangulation(v_x6p, v_y6p, v_triang6p)
        tri_triang1p = matplotlib.tri.Triangulation(self.mesh_xy[:, 0], self.mesh_xy[:, 1], v_triang1p)

        # l_triang1p = copy.deepcopy(self.triang1p)
        # l_triang6p = copy.deepcopy(self.triang6p)
        # self.triang6p = None
        # self.triang1p = None


        self._save_triangulation_plots(tri_triang1p, tri_triang6p)

        # building interpolators: triang1p for the finder, triang6p for the values
        # TODO: could be more efficient only interpolating the fields which are ultimately to be used?
        # create rectangular arrays corresponding to the v_x, v_y grids

        # There might be a cleaner way of doing this
        v_x_flat = self.xy_raw['m_x'].flatten('F') - shiftx
        v_y_flat = self.xy_raw['m_y'].flatten('F') - shifty
        finder = matplotlib.tri.TrapezoidMapTriFinder(tri_triang1p)
        self.interper_f = lambda x: matplotlib.tri.LinearTriInterpolator(
            tri_triang6p, x, trifinder=finder)(v_x_flat, v_y_flat).reshape(self.n_pts_x, self.n_pts_y)









class Mode:
    '''This is a base class for both EM and AC modes.'''

    def __init__(self, sim, m):
        self.mode_num = m
        self.sim_result = sim
        self.fracs = []  # fx, fy, ft, fz
        self.r0 = None  # centre of mass
        self.w2 = None  # second moment width
        self.r0_offset = (0.0, 0.0)
        self.extra_data = {}
        self.analysed = False
        self.interpolated = { FieldType.EM_E: False,  FieldType.EM_H: False, FieldType.AC: False }
        self.d_fields = {}
        self.clear_mode_plot_data()

    def get_mode_helper(self):
        return self.sim_result.get_mode_helper()

    def prepare_mode(self, n_points, field_type):

        mh = self.get_mode_helper()

        mh.setup_plot_grid(n_points=n_points)

        if self.is_AC():
            self.field_type = FieldType.AC
        else:
            self.field_type = field_type


        if not self.interpolated[field_type]:
            self.interpolate_mode(mh)
            self.interpolated[field_type] = True



    def plot_mode(self, comps, field_type=FieldType.EM_E, ax=None,
                  n_points=501, decorator=None):  # TODO get this random parameters hooked better into mode_helper.plot_params


        self.prepare_mode(n_points, field_type)


        mh = self.get_mode_helper()


        # FIX ME
        if not decorator is None:
            mh.plot_params['decorator'] = decorator
        elif mh.plot_params['decorator'] is None:
            # don't want to do this.
            mh.plot_params['decorator'] = plotmodes.Decorator()

        # Just for now
        # mh.plot_params['decorator'].set_singleplot_axes_property('axes.linewidth',.5)
        # mh.plot_params['quiver_points']=6
        # mh.plot_params['colorbar']=False
        # mh.plot_params['add_title']=False

        self._plot_me(mh, comps, field_type, ax)

        self.clear_mode_plot_data()

    def plot_mode_H(self, comps):  # plot magnetic field for EM modes
        self.plot_mode(comps, field_type=FieldType.EM_H)

    def plot_strain(self):
        if not self.sim_result.is_AC():
            print("Doing strain in an EM sim.!")
        print('doing strain')
        mh = self.get_mode_helper()
        mh.plot_strain_mode_i(self.mode_num)

    def clear_mode_plot_data(self):
        for k in self.d_fields.keys(): self.d_fields[k] = None

    def interpolate_mode(self, mode_helper):
        mh = mode_helper
        self.d_fields = mh.interpolate_mode_i(self.mode_num, self.field_type)

        if self.field_type == FieldType.EM_H:  # scale H fields by Z0 to get common units and amplitude with E
            for m_F in self.d_fields.values():   # Do this when they are first made
                m_F *= vacuum_impedance_Z0


    def _plot_me(self, mode_helper, comps, field_type, ax=None):

        # TODO: weirdly, we only ax != None when there is one component to plot
        if not ax is None and len(comps) != 1:
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
            # options are ['Ex', 'Hx', 'ux', 'Ey', 'Hy', 'uy', 'Ez', 'Hz', 'uz','Eabs', 'Habs', 'uabs', 'Et', 'Ht', 'ut']
            for comp in comps:
                cc = component_t(comp)
                plotmodes.plot_one_component(mh.xy_out, self.d_fields, mh.plot_params, self.mode_num, cc, ax)

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
            print('mode has not being analysed')
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

    def set_r0_offset(self, x0, y0):
        '''Sets the transverse position in the grid that is to be regarded as the origin for calculations of center-of-mass.

           This can be useful in aligning the FEM coordinate grid with a physically sensible place in the waveguide.

           :param float x0: *x* position of nominal origin.
           :param float y0: *y* position of nominal origin.
      '''
        self.r0_offset = (x0, y0)

    def analyse_mode(self, n_points=501, EM_field=FieldType.EM_E):
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

        self.prepare_mode(n_points, EM_field)

        self.analysed = True

        mh = self.get_mode_helper()
        v_x = mh.xy_out['v_x']
        v_y = mh.xy_out['v_y']

        mFs = self.d_fields

        m_Fx2 = mFs['Fxr']**2 + mFs['Fxi']**2
        m_Fy2 = mFs['Fyr']**2 + mFs['Fyi']**2
        m_Fz2 = mFs['Fzr']**2 + mFs['Fzi']**2
        m_Fall2 = m_Fx2 + m_Fy2 + m_Fz2

        s_fx = int2d(m_Fx2)
        s_fy = int2d(m_Fy2)
        s_fz = int2d(m_Fz2)

        s_f = s_fx+s_fy+s_fz
        f_x = s_fx/s_f
        f_y = s_fy/s_f
        f_t = f_x+f_y
        f_z = s_fz/s_f
        self.fracs = [f_x, f_y, f_t, f_z]


        [m_x, m_y] = np.meshgrid(v_x, v_y, indexing='ij')  # This is opposite to normal indexing to get image style ordering
        m_yud = np.flipud(m_y)  # Flipping upside down y to get sensible values for r0 position.

        m_xmod = m_x * m_Fall2  # could do this by broadcasting without meshgrid?
        m_ymod = m_yud * m_Fall2

        x0 = int2d(m_xmod)/s_f
        y0 = int2d(m_ymod)/s_f
        m_x2mod = np.power((m_x-x0), 2) * m_Fall2
        m_y2mod = np.power((m_yud-y0), 2) * m_Fall2
        w2x = sqrt(int2d(m_x2mod)/s_f)
        w2y = sqrt(int2d(m_y2mod)/s_f)
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
