# mode_calcs.py is a subroutine of NumBAT that contains methods to
# calculate the EM and Acoustic modes of a structure.

# Copyright (C) 2017 Bjorn Sturmberg, Kokou Dossou.

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


import numpy as np
import sys
import os
import copy
sys.path.append("../backend/")
from math import *

import matplotlib

import plotting
import integration
from nbtypes import *
from reporting import *

from fortran import NumBAT

VacCSpeed=299792458

def load_simulation(prefix):  # allow this function to be module leve rather than class-level
    return Simulation.load_simulation(prefix)

class ModePlotHelper(object): # helper class for plotting. factors common info from Simulation that each mode can draw on, but we only need to do once for each Sim

    def __init__(self, sim):# , field_type):
      self.sim=sim
      self.setup_for_npoints=0

      #self.field_type = field_type  # Enum.FieldType
      self.plot_params = {}

      self.init_arrays()
      self.set_plot_params(prefix='')


    def init_arrays(self):
      self.triang6p = None
      self.triang1p = None
      self.interper_f = None
      self.v_x = None
      self.v_y = None
      self.m_X = None
      self.m_Y = None
      self.v_x6p = None
      self.v_y6p = None
      self.v_Fx6p = None
      self.v_Fy6p = None
      self.v_Fz6p = None

    def cleanup(self):
      # Now that these are part of the sim object,
      # we need to get rid of them to allow pickling of sim files, as they contain C++ objects

      del self.triang6p
      del self.triang1p
      del self.interper_f

      del self.v_x, self.v_y, self.m_X, self.m_Y
      del self.v_x6p, self.v_y6p
      del self.v_Fx6p, self.v_Fy6p, self.v_Fz6p

      self.init_arrays()


#    def make_picklable(self): #before this object can be copied or saved it needs internal C++ objects removed. They are rebuilt when needed.
#        if not self.triang6p is None: del self.triang6p._cpp_triangulation
#        if not self.triang1p is None: del self.triang1p._cpp_triangulation

    def set_plot_params(self, #n_points=501, quiver_points=30, 
                  xlim_min=0, xlim_max=0, ylim_min=0, ylim_max=0,
                  EM_AC='EM_E', 
                  quiver_points=30,
                  num_ticks=None, ticks=False, colorbar=True, contours=False, contour_lst=None,
                  suppress_imimre=True, pdf_png='png', 
                  prefix='tmp', suffix='', decorator=plotting.Decorator(), ):
                  #modal_gains_PE=None,
                  #modal_gains_MB=None,
                  #modal_gains=None):


      if len(prefix) and not os.path.exists("%s-fields" % prefix): 
          os.mkdir("%s-fields" % prefix)

      if type(EM_AC) is type('a'): # aim to get rid of this
        if self.sim.is_AC():
            EM_AC=FieldType.AC
        else:
          try:
            EM_AC=FieldType.from_str(EM_AC)  # TODO:ugly that this changes from string to enum
          except:
            raise ValueError("EM_AC must be either 'AC', 'EM_E' or 'EM_H'.")


      self.plot_params={'xlim_min': xlim_min, 'xlim_max': xlim_max, 'ylim_min': ylim_min, 
                 'ylim_max': ylim_max, 'ticks': ticks, 'num_ticks':num_ticks,
                  'colorbar':colorbar, 'contours':contours, 'contour_lst':contour_lst, 'EM_AC':EM_AC,
                  'prefix': prefix, 'suffix': suffix, 'pdf_png': pdf_png, 
                 # 'modal_gain':modal_gain, 
                  'decorator': decorator,
                  'suppress_imimre':suppress_imimre,
              # 'n_pts_x': n_pts_x, 'n_pts_y': n_pts_y, 
               'quiver_points': quiver_points 
               }

    def interpolate_mode_i(self, ival, field_type):
      # self.v_Fx6p etc could propbably be made local to this function
      # construct the meshed field from fortran solution
      sim = self.sim
      i = 0
      for i_el in np.arange(sim.n_msh_el):
          for i_node in np.arange(6):
              if field_type == FieldType.EM_E or field_type == FieldType.AC :
                  self.v_Fx6p[i] = sim.sol1[0,i_node,ival,i_el]
                  self.v_Fy6p[i] = sim.sol1[1,i_node,ival,i_el]
                  self.v_Fz6p[i] = sim.sol1[2,i_node,ival,i_el]
              if field_type == FieldType.EM_H:
                  self.v_Fx6p[i] = sim.sol1_H[0,i_node,ival,i_el]
                  self.v_Fy6p[i] = sim.sol1_H[1,i_node,ival,i_el]
                  self.v_Fz6p[i] = sim.sol1_H[2,i_node,ival,i_el]
              i += 1

      self.v_F6p = np.sqrt(np.abs(self.v_Fx6p)**2 + np.abs(self.v_Fy6p)**2 + np.abs(self.v_Fz6p)**2)

      # Construct rectangular  interpolated fields.  
      # Always need these ones.
      m_ReFx=self.interper_f(self.v_Fx6p.real)
      m_ReFy=self.interper_f(self.v_Fy6p.real)
      m_ImFz=self.interper_f(self.v_Fz6p.imag)
      m_AbsF=self.interper_f(self.v_F6p)

      # often not needed for plotting, but are used for measuring fractions. (Could fix taht?)
      m_ImFx=self.interper_f(self.v_Fx6p.imag)
      m_ImFy=self.interper_f(self.v_Fy6p.imag)
      m_ReFz=self.interper_f(self.v_Fz6p.real)

      return (m_ReFx, m_ImFx, m_ReFy, m_ImFy, m_ReFz, m_ImFz, m_AbsF)

    def setup_plot_grid(self, n_points=501, #quiver_points=30, 
            #xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None): # these zlim_ points are not actually used here
            ):

      if self.setup_for_npoints == n_points: return   #only need to repeat if the grid density changes
      self.setup_for_npoints = n_points

      sim=self.sim

      # field mapping
      x_tmp = []
      y_tmp = []
      for i in np.arange(sim.n_msh_pts):
          x_tmp.append(sim.mesh_xy[0,i])
          y_tmp.append(sim.mesh_xy[1,i])
      x_min = np.min(x_tmp); 
      x_max = np.max(x_tmp)
      y_min = np.min(y_tmp); 
      y_max = np.max(y_tmp)
      area = abs((x_max-x_min)*(y_max-y_min))
      self.n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
      self.n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))

      # unrolling data for the interpolators
      self.table_nod = sim.table_nod.T
      self.mesh_xy = sim.mesh_xy.T

      # dense triangulation with multiple points
      self.v_x6p = np.zeros(6*sim.n_msh_el)
      self.v_y6p = np.zeros(6*sim.n_msh_el)
      self.v_Fx6p = np.zeros(6*sim.n_msh_el, dtype=np.complex128)
      self.v_Fy6p = np.zeros(6*sim.n_msh_el, dtype=np.complex128)
      self.v_Fz6p = np.zeros(6*sim.n_msh_el, dtype=np.complex128)
      self.v_triang6p = []


      #i = 0  # what is this counting?
      for i_el in np.arange(sim.n_msh_el):
          # triangles
          idx = np.arange(6*i_el, 6*(i_el+1))
          triangles = [[idx[0], idx[3], idx[5]],
                       [idx[1], idx[4], idx[3]],
                       [idx[2], idx[5], idx[4]],
                       [idx[3], idx[4], idx[5]]]
          self.v_triang6p.extend(triangles)



      # Create vectors v_x6p, v_y6p which are unwrapped points at nodes of each element
      # i is the index for the coordinates FIND A BETTER NAME
      i = 0  
      for i_el in np.arange(sim.n_msh_el):
          for i_node in np.arange(6):
              i_ex = self.table_nod[i_el, i_node]-1
              # values
              self.v_x6p[i] = self.mesh_xy[i_ex, 0]
              self.v_y6p[i] = self.mesh_xy[i_ex, 1]    # Fact that this is mesh_xy and not y_arr seems to be right
              i += 1


      ### Interpolate onto triangular grid - honest to FEM elements
      # dense triangulation with unique points
      self.v_triang1p = []
      for i_el in np.arange(sim.n_msh_el):
          # triangles
          table_nod = self.table_nod
          triangles = [[table_nod[i_el,0]-1,table_nod[i_el,3]-1,table_nod[i_el,5]-1],
                       [table_nod[i_el,1]-1,table_nod[i_el,4]-1,table_nod[i_el,3]-1],
                       [table_nod[i_el,2]-1,table_nod[i_el,5]-1,table_nod[i_el,4]-1],
                       [table_nod[i_el,3]-1,table_nod[i_el,4]-1,table_nod[i_el,5]-1]]
          self.v_triang1p.extend(triangles)


      # triangulations:  x and y coords of all points, list of triangles defined by triples of indices of the points
      self.triang6p = matplotlib.tri.Triangulation(self.v_x6p, self.v_y6p, self.v_triang6p)
      self.triang1p = matplotlib.tri.Triangulation(self.mesh_xy[:,0], self.mesh_xy[:,1], self.v_triang1p)

      # Now use the coords user would like to think in
      shiftx, shifty = self.sim.get_xyshift()  # Odd that we are doing this at ModeHelper level?

      self.v_x = np.linspace(x_min, x_max, self.n_pts_x) + shiftx
      self.v_y = np.linspace(y_min, y_max, self.n_pts_y) + shifty
      self.m_X, self.m_Y = np.meshgrid(self.v_x, self.v_y)

      #TODO: this would be much more useful *before* the FEM calc
      print('''
  Structure has raw domain(x,y)   = [{0:.5f}, {1:.5f}] x [ {2:.5f}, {3:.5f}] (um), 
                mapped to (x',y') = [{4:.5f}, {5:.5f}] x [ {6:.5f}, {7:.5f}] (um)
                    '''.format( 
                        1e6*(self.v_x[0]-shiftx), 
                        1e6*(self.v_x[-1]-shiftx), 
                        1e6*(self.v_y[0]-shifty), 1e6*(self.v_y[-1]-shifty),
          1e6*self.v_x[0], 1e6*self.v_x[-1], 1e6*self.v_y[0], 1e6*self.v_y[-1]))

      # building interpolators: triang1p for the finder, triang6p for the values
      #TODO: could be more efficient only interpolating the fields which are ultimately to be used?
      # create rectangular arrays corresponding to the v_x, v_y grids

      v_x_flat=self.m_X.flatten('F') -shiftx  # There might be a cleaner way of doing this
      v_y_flat=self.m_Y.flatten('F') -shifty
      finder = matplotlib.tri.TrapezoidMapTriFinder(self.triang1p)
      self.interper_f = lambda x: matplotlib.tri.LinearTriInterpolator(
           self.triang6p, x, trifinder=finder)(
                   v_x_flat, v_y_flat).reshape(self.n_pts_x, self.n_pts_y)


    def plot_strain_mode_i(self, ival):
             #TODO: this interpolation looks very old. Can we get strain directly from fortran?
     ### Interpolate onto rectangular Cartesian grid

     m_Fx = self.m_ReFx + 1j*self.m_ImFx
     m_Fy = self.m_ReFy + 1j*self.m_ImFy
     m_Fz = self.m_ReFz + 1j*self.m_ImFz
     dx=self.v_x[1]-self.v_x[0]
     dy=self.v_y[1]-self.v_y[0]


     print('finding gradients') 
     #TODO: Check that the axis choice checks out
     del_x_Fx = np.gradient(m_Fx, dx, axis=0)
     del_y_Fx = np.gradient(m_Fx, dy, axis=1)
     del_x_Fy = np.gradient(m_Fy, dx, axis=0)
     del_y_Fy = np.gradient(m_Fy, dy, axis=1)
     del_x_Fz = np.gradient(m_Fz, dx, axis=0)
     del_y_Fz = np.gradient(m_Fz, dy, axis=1)
     del_z_Fx = 1j*self.sim.q_AC*m_Fx
     del_z_Fy = 1j*self.sim.q_AC*m_Fy
     del_z_Fz = 1j*self.sim.q_AC*m_Fz

     
     return

     self.v_x=np.linspace(x_min, x_max, self.n_pts_x)
     # For now, get these avlues from v_x already figured out earlier.
     x_min = self.v_x[0]
     x_max = self.v_x[-1]
     n_pts_x = len(self.v_x)
     y_min = self.v_y[0]
     y_max = self.v_y[-1]
     n_pts_y = len(self.v_y)

     xy = list(zip(self.v_x6p, self.v_y6p))

     # This seems to be equivalent  to taking grid_x = self.m_Y, grid_y = self.m_X
     # CONFIRM!
     #grid_x, grid_y = np.mgrid[x_min:x_max:n_pts_x*1j, y_min:y_max:n_pts_y*1j]  #OLD CODE
     grid_x, grid_y = self.m_Y, self.m_X                                         #NEW CODE

     m_ReFx = interpolate.griddata(xy, v_Fx6p.real, (grid_x, grid_y), method='linear')
     m_ReFy = interpolate.griddata(xy, v_Fy6p.real, (grid_x, grid_y), method='linear')
     m_ReFz = interpolate.griddata(xy, v_Fz6p.real, (grid_x, grid_y), method='linear')
     m_ImFx = interpolate.griddata(xy, v_Fx6p.imag, (grid_x, grid_y), method='linear')
     m_ImFy = interpolate.griddata(xy, v_Fy6p.imag, (grid_x, grid_y), method='linear')
     m_ImFz = interpolate.griddata(xy, v_Fz6p.imag, (grid_x, grid_y), method='linear')
     m_AbsF = interpolate.griddata(xy, v_F6p.real, (grid_x, grid_y), method='linear')

     dx = grid_x[-1,0] - grid_x[-2,0]
     dy = grid_y[0,-1] - grid_y[0,-2]

     m_Fx = m_ReFx + 1j*m_ImFx
     m_Fy = m_ReFy + 1j*m_ImFy
     m_Fz = m_ReFz + 1j*m_ImFz
     m_Fx = m_Fx.reshape(n_pts_x,n_pts_y)
     m_Fy = m_Fy.reshape(n_pts_x,n_pts_y)
     m_Fz = m_Fz.reshape(n_pts_x,n_pts_y)
     m_AbsF = m_AbsF.reshape(n_pts_x,n_pts_y)

     m_ReFx = np.real(m_Fx)
     m_ReFy = np.real(m_Fy)
     m_ReFz = np.real(m_Fz)
     m_ImFx = np.imag(m_Fx)
     m_ImFy = np.imag(m_Fy)
     m_ImFz = np.imag(m_Fz)

     del_x_Fx = np.gradient(m_Fx, dx, axis=0)
     del_y_Fx = np.gradient(m_Fx, dy, axis=1)
     del_x_Fy = np.gradient(m_Fy, dx, axis=0)
     del_y_Fy = np.gradient(m_Fy, dy, axis=1)
     del_x_Fz = np.gradient(m_Fz, dx, axis=0)
     del_y_Fz = np.gradient(m_Fz, dy, axis=1)
     del_z_Fx = 1j*sim_wguide.q_AC*m_Fx
     del_z_Fy = 1j*sim_wguide.q_AC*m_Fy
     del_z_Fz = 1j*sim_wguide.q_AC*m_Fz

     # Flip y order as imshow has origin at top left
     del_mat = np.array([del_x_Ex[:,::-1].real, del_x_Ey[:,::-1].real, del_x_Ez[:,::-1].real, del_x_Ex[:,::-1].imag, del_x_Ey[:,::-1].imag, del_x_Ez[:,::-1].imag, del_y_Ex[:,::-1].real, del_y_Ey[:,::-1].real, del_y_Ez[:,::-1].real, del_y_Ex[:,::-1].imag, del_y_Ey[:,::-1].imag, del_y_Ez[:,::-1].imag, del_z_Ex[:,::-1].real, del_z_Ey[:,::-1].real, del_z_Ez[:,::-1].real, del_z_Ex[:,::-1].imag, del_z_Ey[:,::-1].imag, del_z_Ez[:,::-1].imag])
     v_labels = ["Re($S_{xx}$)","Re($S_{xy}$)","Re($S_{xz}$)", "Im($S_{xx}$)","Im($S_{xy}$)","Im($S_{xz}$)","Re($S_{yx}$)","Re($S_{yy}$)","Re($S_{yz}$)","Im($S_{yx}$)","Im($S_{yy}$)","Im($S_{yz}$)","Re($S_{zx}$)","Re($S_{zy}$)","Re($S_{zz}$)","Im($S_{zx}$)","Im($S_{zy}$)","Im($S_{zz}$)"]

     # stress field plots
     plt.clf()
     fig = plt.figure(figsize=(15,30))
     for i_p,plot in enumerate(del_mat):
         ax = plt.subplot(6,3,i_p+1)
         im = plt.imshow(plot.T);
         # no ticks
         plt.xticks([])
         plt.yticks([])
         # limits
         if xlim_min >0 :
             ax.set_xlim(xlim_min*n_points,(1-xlim_max)*n_points)
         if ylim_min >0 :
             ax.set_ylim((1-ylim_min)*n_points,ylim_max*n_points)
         # titles
         plt.title(v_labels[i_p],fontsize=decorator.get_font_size('subplot_title'))
         # colorbar
         divider = make_axes_locatable(ax)
         cax = divider.append_axes("right", size="5%", pad=0.1)
         cbar = plt.colorbar(im, cax=cax, format='%.2e')
         if num_ticks:
             cbarticks = np.linspace(np.min(plot), np.max(plot), num=num_ticks)                
         elif ylim_min != 0:
             if xlim_min/ylim_min > 3:
                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=3)
             if xlim_min/ylim_min > 1.5:
                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=5)
             else:
                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
         else:
             cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
         cbar.set_ticks(cbarlabels)
         cbarlabels = ['%.2f' %t for t in cbarlabels]
         cbar.set_ticklabels(cbarlabels)
         if contours:
             if contour_lst:
                 cbarticks = contour_lst
             if np.max(np.abs(plot[~np.isnan(plot)])) > plot_threshold:
                 CS2 = ax.contour(m_X, m_Y, plot.T, levels=cbarticks, colors=colors[::-1], linewidths=(1.5,))
             cbar.add_lines(CS2)
         cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))
     fig.set_tight_layout(True)
     n_str = ''
     if np.imag(sim_wguide.Eig_values[ival]) < 0:
         k_str = r'$\Omega/2\pi = %(re_k)f %(im_k)f i$ GHz'% \
             {'re_k' : np.real(sim_wguide.Eig_values[ival]*1e-9),
             'im_k' : np.imag(sim_wguide.Eig_values[ival]*1e-9)}
     else:
         k_str = r'$\Omega/2\pi = %(re_k)f + %(im_k)f i$ GHz'% \
             {'re_k' : np.real(sim_wguide.Eig_values[ival]*1e-9),
             'im_k' : np.imag(sim_wguide.Eig_values[ival]*1e-9)}
     plt.suptitle('Mode #' + str(ival) + '   ' + k_str + '   ' + n_str, fontsize=decorator.get_font_size('title'))

     if pdf_png=='png':
         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.png' %
             {'pre' : prefix, 's' : EM_AC, 'i' : ival, 'add' : suffix})
     elif pdf_png=='pdf':
         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.pdf' %
             {'pre' : prefix, 's' : EM_AC, 'i' : ival, 'add' : suffix}, bbox_inches='tight')
     if not keep_plots_open: plt.close()


      



class Mode(object):
  '''This is a base class for both EM and AC modes.''' 
  def __init__(self, sim, m):
    self.mode_num=m
    self.sim=sim
    self.fracs=[]  # fx, fy, ft, fz
    self.r0=None  # centre of mass
    self.w2=None  # second moment width
    self.r0_offset=(0.0, 0.0)
    self.extra_data={}
    self.analysed = False
    self.interpolated = False
    self.clear_mode_plot_data()

  def _get_mode_helper(self):
      return self.sim.get_mode_helper()

  def analyse_mode(self, n_points=501, EM_field=FieldType.EM_E):
      mh = self._get_mode_helper()
      mh.setup_plot_grid(n_points=n_points)  # TODO: make all this cleaner. this should happen once for all modes
      if not self.interpolated:

          #self.field_type=FieldType.EM_E  #REMOVE THIS REPEAT WITH plot_mode
          if self.is_AC(): 
             self.field_type=FieldType.AC
          else:
             self.field_type=EM_field

          self.interpolate_mode(mh)
      self._analyse_mode(mh.v_x, mh.v_y, 
              self.m_ReFx, self.m_ReFy, self.m_ReFz, self.m_ImFx, self.m_ImFy, self.m_ImFz, self.m_AbsF)

  def plot_mode(self, comps, EM_field=FieldType.EM_E, ax=None, 
          n_points=501, decorator=None): #TODO get this random parameters hooked better into mode_helper.plot_params

      self.field_type=FieldType.EM_E
      if self.is_AC(): 
          self.field_type=FieldType.AC
      else:
          self.field_type=EM_field

      mh = self._get_mode_helper()

      mh.setup_plot_grid(n_points=n_points)  # TODO: make all this cleaner. this should happen once for all modes

      #mh.set_plot_params() # can delete shortly

      #print(mh.plot_params)
      #print('dec', mh.plot_params['decorator'])

      #FIX ME
      if not decorator is None: 
          mh.plot_params['decorator']=decorator 
      elif mh.plot_params['decorator'] is None:
          mh.plot_params['decorator']=plotting.Decorator() # don't want to do this.


      ## Just for now
      #mh.plot_params['decorator'].set_singleplot_axes_property('axes.linewidth',.5)
      #mh.plot_params['quiver_points']=6
      #mh.plot_params['colorbar']=False
      #mh.plot_params['add_title']=False


      self.interpolate_mode(mh)
      self._plot_me(mh, comps, self.field_type, ax)
      self.clear_mode_plot_data()

  def plot_mode_H(self, comps): # plot magnetic field for EM modes
      self.plot_mode(comps, EM_field=FieldType.EM_H)

  def plot_strain(self):
     if not self.sim.is_AC():
         print("Doing strain in an EM sim.!")
     print('doing strain')
     mode_helper = self._get_mode_helper()
     mode_helper.plot_strain_mode_i(self.mode_num)

  def clear_mode_plot_data(self):
      self.m_ReFx = None
      self.m_ImFx = None
      self.m_ReFy = None
      self.m_ImFy = None
      self.m_ReFz = None
      self.m_ImFz = None
      self.m_AbsF = None


  def interpolate_mode(self, mode_helper):
      self.interpolated = True
      sim=self.sim
      ival=self.mode_num

      mh = mode_helper
      
      (self.m_ReFx, self.m_ImFx, self.m_ReFy, self.m_ImFy, 
              self.m_ReFz, self.m_ImFz, self.m_AbsF)=mh.interpolate_mode_i(ival, self.field_type)

      if self.field_type == FieldType.EM_H:  # scale H fields by Z0 to get common units and amplitude with E
          mu0 = 1.25663706212e-6
          eps0 =  8.8541878128e-12
          Z0 = sqrt(mu0/eps0)
          self.m_ReFx *= Z0
          self.m_ReFy *= Z0
          self.m_ReFz *= Z0
          self.m_ImFx *= Z0
          self.m_ImFy *= Z0
          self.m_ImFz *= Z0
          self.m_AbsF *= Z0


  def _plot_me(self, mode_helper, comps, field_type, ax=None):
      v_plots = {'Fxr':self.m_ReFx, 'Fyr':self.m_ReFy, 'Fzi':self.m_ImFz, 
              'Fxi':self.m_ImFx, 'Fyi':self.m_ImFy, 'Fzr':self.m_ReFz, 'Fabs':self.m_AbsF}
     

      # TODO: weirdly, we only ax != None when there is one component to plot
      if not ax is None and len(comps) !=1:
          print('\nError: when providing an axis to plot on, must specify exactly one modal component.')
          return

      mh = mode_helper
      decorator=mh.plot_params['decorator']

      decorator._set_for_multi()
      mh.plot_params['EM_AC']=field_type # TODO this is a kludgy way of doing this. send it through separately

      if ax is None: # can't do multiplots on a provided axis (would need a provided figure)
        plotting.plot_all_components(mh.v_x, mh.v_y, mh.m_X, mh.m_Y, v_plots, 
          mh.plot_params, self.sim, self.mode_num) 


      if len(comps):
        decorator._set_for_single()
        for comp in comps: #options are ['Ex', 'Hx', 'ux', 'Ey', 'Hy', 'uy', 'Ez', 'Hz', 'uz','Eabs', 'Habs', 'uabs', 'Et', 'Ht', 'ut'] 
          cc=component_t(comp)
          plotting.plot_one_component(mh.m_X, mh.m_Y, v_plots, mh.plot_params, self.sim, self.mode_num, cc, ax)

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
    if not self.analysed: print('mode has not being analysed')
    return self.fracs

  def __str__(self): 
    '''String representation of the mode.'''
    s='Abstract mode class'
    return s

  def is_poln_ex(self):
    '''Returns true if mode is predominantly x-polarised (ie if *fx*>0.7).
       
       :rtype: bool
       '''
    polthresh=.7
    return self.fracs[0]>polthresh

  def is_poln_ey(self):
    '''Returns true if mode is predominantly y-polarised (ie if *fy*>0.7).
       
       :rtype: bool
       '''
    polthresh=.7
    return self.fracs[1]>polthresh

  def is_poln_indeterminate(self):
    '''Returns true if transverse polarisation is neither predominantly *x* or *y* oriented.

       :rtype: bool
       '''
    return not (self.is_poln_ex() or self.is_poln_ey())

  def is_EM(self):
    '''Returns true if the mode is an electromagnetic mode.

       :rtype: bool
       '''
    return self.sim.is_EM()

  def is_AC(self):
    '''Returns true if the mode is an acoustic mode.

       :rtype: bool
       '''
    return self.sim.is_AC()

  def center_of_mass(self): 
    '''Returns the centre of mass of the mode relative to the specified origin.

       :rtype: float
    '''
    return self.r0-self.r0_offset

  def second_moment_widths(self): 
    '''Returns the second moment widths :math:`(w_x, w_y, \sqrt{w_x^2+w_y^2})` of the mode relative to the specified origin.

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
    '''Returns the combined second moment width :math:`\sqrt{w_x^2+w_y^2}`.

       :rtype: float
       '''
    return self.w2[2]

  def set_r0_offset(self, x0, y0):
    '''Sets the transverse position in the grid that is to be regarded as the origin for calculations of center-of-mass.
  
       This can be useful in aligning the FEM coordinate grid with a physically sensible place in the waveguide.

       :param float x0: *x* position of nominal origin.
       :param float y0: *y* position of nominal origin.
  '''
    self.r0_offset=(x0, y0)

  def _analyse_mode(self, v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf):
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
    self.analysed = True
    s_fx=np.sum(np.sum(m_Refx*m_Refx+m_Imfx*m_Imfx))
    s_fy=np.sum(np.sum(m_Refy*m_Refy+m_Imfy*m_Imfy)) 
    s_fz=np.sum(np.sum(m_Refz*m_Refz+m_Imfz*m_Imfz)) 
    s_f=s_fx+s_fy+s_fz
    f_x=s_fx/s_f
    f_y=s_fy/s_f
    f_t=f_x+f_y
    f_z=s_fz/s_f
    self.fracs=[f_x, f_y, f_t, f_z]
    [m_x, m_y]=np.meshgrid(v_x, v_y, indexing='ij') # TODO: IS THIS RIGHT?
    #print('vs', v_x[0], v_x[1], v_x[-1], v_y[0], v_y[1], v_y[-1])
    #print('sh', len(v_x), len(v_y), m_Refx.shape, m_x.shape)
    #print('m_xa', m_x[0,:5])
    #print('m_xb', m_x[:5,0])
    #print('m_xb', m_x[-5:,0])
    #print('m_ya', m_y[0,:5])
    #print('m_yb', m_y[0,-5:])
    #print('m_yc', m_y[:5,0])


    m_mod= m_Refx*m_Refx+m_Imfx*m_Imfx+m_Refy*m_Refy+m_Imfy*m_Imfy+m_Refz*m_Refz+m_Imfz*m_Imfz
    m_xmod= m_x * m_mod  # could do this by broadcasting without meshgrid?
    m_yud = np.flipud(m_y)
    m_ymod= m_yud * m_mod   # Flipping upside down y to get sensible values for r0 position. 
                            # Should this happen much closer to the extraction from FEM?
    x0=np.sum(np.sum(m_xmod))/s_f
    y0=np.sum(np.sum(m_ymod))/s_f
    m_x2mod= np.power((m_x-x0),2) * m_mod
    m_y2mod= np.power((m_yud-y0),2) * m_mod
    w2x=sqrt(np.sum(np.sum(m_x2mod))/s_f)
    w2y=sqrt(np.sum(np.sum(m_y2mod))/s_f)
    w2=sqrt(w2x*w2x+w2y*w2y)
#print ('sums', s_f, np.sum(np.sum(m_mod)), w2x, w2y, w2)
    self.r0=np.array([x0, y0])
    self.w2=np.array([w2x, w2y, w2])

class ModeEM(Mode):
  '''Class representing a single electromagnetic (EM) mode.'''
  def __init__(self, sim, m):
    super().__init__(sim, m)

  def __str__(self): 
    s='EM mode # {0}'.format(self.mode_num)
    return s

  def _analyse_mode(self, v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf):
    super()._analyse_mode(v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf)


class ModeAC(Mode):
  '''Class representing a single acoustic (AC) mode.'''
  def __init__(self, sim, m):
    super().__init__(sim, m)

    self.gain={}  # { (EM_p_i, EM_s_j): gain}
    self.gain_PE={}
    self.gain_MB={}

  def __str__(self): 
    s='AC mode # {0}'.format(self.mode_num)
    return s

  def _analyse_mode(self, v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf):
    super()._analyse_mode(v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf)

    

class Simulation(object):
    '''Class for calculating the electromagnetic and/or acoustic modes of a ``Struct`` object.
    '''

    def __init__(self, structure, num_modes=20, wl_nm=1, n_eff=None, shift_Hz=None, 
                 q_AC=None, EM_sim=None, Stokes=False, 
                 calc_EM_mode_energy=False, calc_AC_mode_power=False, debug=False):
        '''Sets up the problem for the mode calculation at a given optical wavelength `wl_nm` or acoustic wavenumber `q_AC`.

           For electromagnetic problems, the tool solves for the effective index or wavenumber at a given wavelength.
           For acoustic problems, the tool solves for the acoustic frequency at a given wavenumber.

             :param Simulation structure: The waveguide structure to be solved.
             :param int num_modes: The number of modes to be found.
             :param float wl_nm: For electromagnetic problems, the vacuum wavelength in nanometers.
             :param float n_eff: For electromagnetic problems, an estimated effective index to begin the eigenvalue search.
             :param float shift_Hz: For acoustic problems, an estimated frequency offset to begin the eigenvalue search.
             :param float q_AC: For acoustic problems, the acoustic wavevector of the mode.
             :param float EM_sim: For acoustic problems, the results of a previously solved EM problem to speed calculations. 
             :param bool calc_EM_mode_energy: For electromagnetic problems, whether to calculate the optical mode energy.
             :param bool calc_AC_mode_power: For acoustic problems, whether to calculate the acoustic mode power.
          '''

        self.structure = structure
        self.mode_plot_helper = None
        self.sim_type = SimType.EM

        self.lambda_m = wl_nm*1e-9
        self.n_eff = n_eff
        self.shift_Hz = shift_Hz

        self.q_AC = q_AC
        self.Omega_AC = None
        self.EM_sim = EM_sim

        self.n_modes = num_modes
        self.Stokes = Stokes
        self.mode_pol = None
        self.k_0 = 2 * np.pi / self.lambda_m
        # just off normal incidence to avoid degeneracies
        self.k_perp = np.array([1e-16, 1e-16])
        speed_c = 299792458
        self.omega_EM = 2*np.pi*speed_c/self.lambda_m # Angular freq in units of rad/s
        self.calc_EM_mode_energy = calc_EM_mode_energy
        self.calc_AC_mode_power = calc_AC_mode_power

        self.EM_mode_energy = None
        self.EM_mode_power = None

        self.AC_mode_energy = None
        self.AC_mode_power = None

        self.debug = debug
        self.sym_reps = None
        self.point_group = PointGroup.Unknown
        self.Q_method = QAcMethod.NotSet
        self.ac_alpha_t = None   # temporal acoustic loss [1/s]
        self.ac_linewidth = None   # acoustic linewidth [Hz]
        self.ac_Qmech = None   # acoustic mechanical Q [dimless]

        self.mode_set=[]
        self.r0_offset = [0, 0] # passed to modes when created 


        self.n_msh_pts = 0   # number of points in .msh mesh file
        self.n_msh_el = 0    # number of elements in .msh mesh file
        self.table_nod = None
        self.type_el = None
        self.type_nod = None
        self.mesh_xy = None
        self.ls_material = None

        # Takes list of all material refractive indices
        # Discards any that are zero
        # Set up mapping table for refractive indices
        # (Why we need this is mystery)
        # el_conv_table_n maps the number of the material to the position in the nonzero v_refindexn
        # el_conv_table_n[ith material] = index into v_refindexn  of non-zero refractive indices
        # Except for zero index materials, 
        #  it will always be {1:1, 2:2, 3:3, .., num_mats:num_mats}
        # (MJS: Not sure about the counting from 1, possibly needed for fortran?)
        v_refindexn = []
        v_refindexn_tmp = np.array([m.refindex_n for m in self.structure.d_materials.values()])
        self.el_conv_table_n = {}
        i = 1; j = 1
        for n in v_refindexn_tmp:
            if n != 0:
                v_refindexn.append(n)
                self.el_conv_table_n[i] = j
                j += 1
            i += 1
        self.v_refindexn = np.array(v_refindexn)
        v_refindexn = None

        if self.structure.loss is False:
            self.v_refindexn = self.v_refindexn.real

    def get_xyshift(self):
        if self.is_EM():
            return self.structure.shift_em_x, self.structure.shift_em_y
        else:
            return self.structure.shift_ac_x, self.structure.shift_ac_y

    def clean_for_save(self):
        if self.mode_plot_helper is not None:
            self.mode_plot_helper.cleanup()

    def save_simulation(self, prefix): 
        self.clean_for_save()

        # Acoustic sims can contain EM sims which must also be clean for saving
        if self.EM_sim is not None: self.EM_sim.clean_for_save() 

        np.savez(prefix, simulation = self)  

    @staticmethod
    def load_simulation(prefix): 
        npzfile = np.load(prefix+'.npz', allow_pickle=True)
        return npzfile['simulation'].tolist()

    def get_mode_helper(self):
        if self.mode_plot_helper is None:
          self.mode_plot_helper=ModePlotHelper(self)

        return self.mode_plot_helper

    def is_EM(self): 
      '''Returns true if the solver is setup for an electromagnetic problem.'''
      return self.sim_type == SimType.EM

    def is_AC(self): 
      '''Returns true if the solver is setup for an acoustic problem.'''
      return self.sim_type == SimType.AC
    
    def analyse_modes(self, n_points=501):
        modes = self.get_modes()
        for m in modes: m.analyse_mode(n_points=n_points)

    def get_modes(self):
      '''Returns an array of class `Mode` containing the solved electromagnetic or acoustic modes.
         
         :rtype: numarray(Mode)
         '''
      if not len(self.mode_set):
        for m in range(self.n_modes):
          if self.is_EM():
            mode=ModeEM(self, m)
          else:
            mode=ModeAC(self, m)
          mode.set_r0_offset(self.r0_offset[0], self.r0_offset[1])  # awkward and specific to do this here, but might have already been set in the Simulation object befores modes are created
          self.mode_set.append(mode)
        
      return self.mode_set

    def set_r0_offset(self, rx, ry): # this is clumsy and only works if called after modes have been calced.
      self.r0_offset = [rx, ry]
      for m in self.get_modes(): m.set_r0_offset(rx, ry)

    def symmetry_classification(self, m):
      '''If the point group of the structure has been specified, returns the symmetry class of the given mode.
         
         :param int m: Index of the mode of interest.
         :rtype: PointGroup
         '''
      if self.point_group == PointGroup.Unknown: return ''
      return '{0}:{1}'.format(self.point_group.name, self.sym_reps[m].name)

    def neff(self, m): 
      ''' Returns the effective index of EM mode `m`.
  
      :param int m: Index of the mode of interest.
      :rtype: float
      ''' 
      assert(self.is_EM())
      return np.real(self.Eig_values[m]*self.lambda_m/(2*np.pi))

    def ngroup_EM_available(self): 
      '''Returns true if a measure of the electromagnetic group index is available.'''
      return not(self.EM_mode_energy is None or self.EM_mode_power is None)

    def vgroup_AC_available(self): 
      '''Returns true if a measure of the acoustic group velocity is available.'''
      return not(self.AC_mode_energy is None or self.AC_mode_power is None)

    def ngroup_EM(self, m):
      '''Returns the group index of electromagnetic mode `m`, if available, otherwise returns zero with a warning message.
         
         :param int m: Index of the mode of interest.
         :return: Group index of the mode.
         :rtype: float
         '''

      if not self.ngroup_EM_available():
        print('''EM group index requires calculation of mode energy and mode power when calculating EM modes. 
               Set calc_EM_mode_energy=True and calc_AC_mode_power=True in call to Simulation''')
        ng = 0

      if abs(self.EM_mode_energy[m]) > 0.0:
        vg= np.real(self.EM_mode_power[m]/self.EM_mode_energy[m])
        ng=VacCSpeed/vg
      else:
        ng = 0

      return ng

    def ngroup_EM_all(self):
      '''Returns a numpy array of the group index of all electromagnetic modes, if available, 
         otherwise returns a zero numarray with a warning message.
         
         :return: numpy array of  index of the mode.
         :rtype: array(float)
         '''
      if not self.ngroup_EM_available():
        print('''EM group index requires calculation of mode energy and mode power when calculating EM modes. 
               Set calc_EM_mode_energy=True in call to calc_EM_modes''')
        return np.zeros(len(self.Eig_values), dtype=float)
      vg= np.real(self.EM_mode_power/self.EM_mode_energy)
      ng=VacCSpeed/vg
      return ng

    def kz_EM(self, m): 
      '''Returns the wavevector in 1/m of electromagnetic mode `m`.
         
         :param int m: Index of the mode of interest.
         :return: Wavevector k in 1/m.
         :rtype: float
         '''
      assert(self.is_EM())
      return self.Eig_values[m]

    def nu_AC(self, m): 
      '''Returns the frequency in Hz of acoustic mode `m`.
         
         :param int m: Index of the mode of interest.
         :return: Frequency :math:`\\nu` in Hz
         :rtype: float
         '''
      assert(self.is_AC())
      return self.Eig_values[m]

    def Omega_AC(self, m): 
      '''Returns the frequency in 1/s of acoustic mode `m`.
         
         :param int m: Index of the mode of interest.
         :return: Angular requency :math:`\\Omega` in Hz
         :rtype: float
         '''
      return self.Eig_values[m]*2*np.pi

    def neff_all(self): 
      ''' Return an array of the effective index of all electromagnetic modes.

         :return: numpy array of effective indices
         :rtype: array(float)
          '''
      assert(self.is_EM())
      return np.real(self.Eig_values*self.lambda_m/(2*np.pi))

    def kz_EM_all(self): 
      ''' Return an array of the wavevector in 1/m of all electromagnetic modes.

         :return: numpy array of wavevectors in 1/m
         :rtype: array(float)
          '''
      assert(self.is_EM())
      return self.Eig_values

    def nu_AC_all(self): 
      ''' Return an array of the frequency in Hz of all acoustic modes.

         :return: numpy array of frequencies in Hz
         :rtype: array(float)
         '''
      assert(self.is_AC())
      return self.Eig_values

    def Omega_AC_all(self): 
      ''' Return an array of the angular frequency in 1/s of all acoustic modes.

         :return: numpy array of angular frequencies in 1/s
         :rtype: array(float)
         '''
      assert(self.is_AC())
      return self.Eig_values*2*np.pi

    def vp_AC(self, m): 
      """ Return phase velocity of all AC modes in m/s"""
      ''' Return the phase velocity in m/s of acoustic mode `m`.

         :return: Phase velocity of acoustic mode `m` in m/s
         :rtype: float
         '''
      assert(self.is_AC())
      return np.pi*2*np.real(self.Eig_values[m])/self.q_AC

    def vp_AC_all(self): 
      ''' Return an array of the phase velocity in m/s of all acoustic modes.

         :return: numpy array of elastic phase velocities in m/s
         :rtype: array(float)
         '''
      assert(self.is_AC())
      return np.pi*2*np.real(self.Eig_values)/self.q_AC

    def vg_AC(self, m): 
      """ Return group velocity of AC mode m in m/s"""
      if self.AC_mode_energy is None or self.AC_mode_power is None:
        print('''AC group velocity requires calculation of mode energy and mode power when calculating AC modes. 
               Set calc_AC_mode_power=True in call to calc_AC_modes''')
        return 0
      vg= np.real(self.AC_mode_power[m]/self.AC_mode_energy[m])
      return vg

    def vg_AC_all(self): 
      assert(self.is_AC())
      """ Return group velocity of all AC modes in m/s"""
      if self.AC_mode_energy is None or self.AC_mode_power is None:
        print('''AC group velocity requires calculation of mode energy and mode power when calculating AC modes. 
               Set calc_AC_mode_power=True in call to calc_AC_modes''')
        return np.zeros(len(self.Eig_values), dtype=float)
      vg= np.real(self.AC_mode_power/self.AC_mode_energy)
      return vg



    def alpha_t_AC(self, m): 
      assert(self.is_AC())
      return self.ac_alpha_t[m]

    def alpha_t_AC_all(self): 
      assert(self.is_AC())
      return self.ac_alpha_t

    def alpha_s_AC(self, m): # spatial loss [1/m]  #TODO: surely this should be the group velocity for scaling between spatial and temporal decays #Which is a problem because it requires knowing vg
      assert(self.is_AC())
      print('alpha_s:')
      print(self.alpha_t_AC(m))
      print(self.vg_AC(m))
      return self.alpha_t_AC(m)/self.vg_AC(m)

    def alpha_s_AC_all(self): # spatial loss [1/m]
      assert(self.is_AC())
      return self.alpha_t_AC_all/self.vg_AC_all

    def Qmech_AC(self, m): 
      assert(self.is_AC())
      return self.ac_Qmech[m]

    def Qmech_AC_all(self):
      assert(self.is_AC())
      return self.ac_Qmech

    def linewidth_AC(self, m):
      assert(self.is_AC())
      return self.ac_linewidth[m]

    def linewidth_AC_all(self):
      assert(self.is_AC())
      return self.ac_linewidth

    def analyse_symmetries(self, ptgrp):
      self.point_group=ptgrp
      symlist = integration.symmetries(self)
      self.sym_reps = []
      if ptgrp == PointGroup.C2V:
        for m, sym in enumerate(symlist):
          if sym == (1,1,1):     self.sym_reps.append(SymRep.A)
          elif sym == (-1,1,-1): self.sym_reps.append(SymRep.B1)
          elif sym == (1,-1,-1): self.sym_reps.append(SymRep.B2)
          elif sym == (-1,-1,1): self.sym_reps.append(SymRep.B3)
          else: 
            print('Warning: Unknown symmetry pattern', sym)
            self.sym_reps.append(SymRep.Unknown)


      else:
        print("unknown symmetry properties in mode_calcs")

    def calc_acoustic_losses(self, fixed_Q=None): # TODO: make sure this is not done more than once on the same Simulation
      alpha = None
      if fixed_Q is None: 
        self.Q_method=QAcMethod.Intrinsic

        # Calc alpha (loss) Eq. 45
        print("Acoustic loss calc")
        nnodes=6 # TODO: is this right?
        try:
            if self.EM_sim.structure.inc_shape in self.EM_sim.structure.linear_element_shapes:
                alpha = NumBAT.ac_alpha_int_v2(self.n_modes,
                    self.n_msh_el, self.n_msh_pts, nnodes,
                    self.table_nod, self.type_el, self.mesh_xy,
                    self.structure.n_typ_el_AC, self.structure.eta_tensor,
                    self.q_AC, self.Omega_AC, self.sol1,
                    # sim_AC.AC_mode_power) # appropriate for alpha in [1/m]
                    self.AC_mode_energy) # appropriate for alpha in [1/s]
            else:
                if self.EM_sim.structure.inc_shape not in self.EM_sim.structure.curvilinear_element_shapes:
                    print("Warning: ac_alpha_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
                Fortran_debug=0
                overlap=np.zeros(self.n_modes, dtype=complex)  # not sure why this is needed by ac_alpha_int
                alpha = NumBAT.ac_alpha_int(self.n_modes,
                    self.n_msh_el, self.n_msh_pts, nnodes,
                    self.table_nod, self.type_el, self.mesh_xy,
                    self.structure.n_typ_el_AC, self.structure.eta_tensor,
                    self.q_AC, self.Omega_AC, self.sol1,
                    # sim_AC.AC_mode_power, Fortran_debug) # appropriate for alpha in [1/m]
                    self.AC_mode_energy, Fortran_debug) # appropriate for alpha in [1/s]
        except KeyboardInterrupt:
            print("\n\n Routine ac_alpha_int interrupted by keyboard.\n\n")
        self.ac_alpha_t = np.real(alpha)
        # Q_factors = 0.5*(q_AC/alpha)*np.ones(n_modes) # appropriate for alpha in [1/m]
        self.ac_Qmech = 0.5*(np.real(self.Omega_AC)/self.ac_alpha_t)*np.ones(self.n_modes) # appropriate for alpha in [1/s]
      else:
        self.Q_method=QAcMethod.Fixed
        # factor of a 1/2 because alpha is for power!
        # alpha [1/m] = Omega_AC/(2*vg*fixed_Q) = q_AC/fixed_Q
        # alpha [1/s] = vg * alpha [1/m]
        # alpha [1/s] = Omega_AC/(2*fixed_Q)
        # alpha = 0.5*(q_AC/fixed_Q)*np.ones(n_modes) # appropriate for alpha in [1/m]
        self.ac_Qmech = fixed_Q*np.ones(self.n_modes)
        self.ac_alpha_t = 0.5*(np.real(self.Omega_AC)/fixed_Q)*np.ones(self.n_modes) # appropriate for alpha in [1/s]

      self.ac_linewidth = self.ac_alpha_t/np.pi # SBS linewidth of each resonance in [Hz]   #TODO: not sure about the 1/pi.  
                                                                                            #If linewdith should be amplitude rate in Hz, wouldn't it be
                                                                                            # alpha/(2 * 2pi)  since alpha is a power decay rate

    def calc_EM_modes(self):
        """ Run a Fortran FEM calculation to find the optical modes.

        Returns a ``Simulation`` object that has these key values:

        Eig_values: a 1d array of Eigenvalues (propagation constants) in [1/m]

        sol1: the associated Eigenvectors, ie. the fields, stored as [field comp, node nu on element, Eig value, el nu]

        EM_mode_power: the power in the optical modes. Note this power is negative for modes travelling in the negative
                       z-direction, eg the Stokes wave in backward SBS.
        """
        self.sim_type = SimType.EM

        tstruc=self.structure
        self.d_in_m = tstruc.unitcell_x*1e-9   # TODO: don't think fortran really needs this. Why does it not care about y dimension?


        if self.n_modes < 20:
            self.n_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")


        # Parameters that control how FEM routine runs
        self.E_H_field = 1  # Selected formulation (1=E-Field, 2=H-Field)
        bnd_cdn_i = 2  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        itermax = 30  # Maximum number of iterations for convergence
        EM_FEM_debug = self.debug  # Fortran routines will display & save add. info

        print(' Boundary conditions: %s'% {0:'Dirichlet', 1:'Neumann', 2:'Periodic'}[bnd_cdn_i])
        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        shift_ksqr = self.n_eff**2 * self.k_0**2

        if EM_FEM_debug == 1:
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")
            if not os.path.exists("Output"):
                os.mkdir("Output")

        
        # mesh sizes is at top of mail file    #TODO: would these values be better as part of tstruc?
        self.n_msh_pts, self.n_msh_el = map(int, tstruc.get_mail_data()[0].split())
            
        print('\n Structure has {0} mesh points and {1} mesh elements.'.format(self.n_msh_pts, self.n_msh_el))

        # Size of Fortran's complex superarray (scales with mesh)
        #int_max, cmplx_max, real_max = NumBAT.array_size(self.n_msh_el, self.n_modes)
        #if EM_FEM_debug == 1:
        #  print("Mesh calculated: %d nodes."%self.n_msh_el)

        EM_FEM_debug = 0
#        try:
        fort_err = 0

        resm = NumBAT.calc_em_modes( self.lambda_m, self.d_in_m, tstruc.mesh_file, 
            self.n_modes, self.n_msh_pts, self.n_msh_el, 
            tstruc.n_typ_el, self.v_refindexn,
            self.k_perp, shift_ksqr, bnd_cdn_i, itermax, self.E_H_field, EM_FEM_debug,
            # TODO: these are all obselete, remove.
            tstruc.plotting_fields, tstruc.plot_real, tstruc.plot_imag, tstruc.plot_abs) 
    #        cmplx_max, real_max, int_max)

        self.Eig_values, self.sol1, self.mode_pol, self.table_nod, \
        self.type_el, self.type_nod, self.mesh_xy, self.ls_material, fort_err, fort_mesg = resm

        if fort_err != 0:
            fort_mesg = str(fort_mesg, 'utf-8')  # fort_mesg comes back as a byte string.
            report_and_exit('Fortran error in solving for electromagnetic modes: \n'
                            ' NumBAT Fortran error code = %d. \n Message: \n %s'%(fort_err, fort_mesg))

#        except KeyboardInterrupt:
#            print("\n\n FEM routine calc_EM_modes",\
#            "interrupted by keyboard.\n\n")


        # if not tstruc.plot_field_conc:
        #     self.mode_pol = None

        # if tstruc.plotting_fields != 1:
        #     self.sol1 = None
        #     self.v_refindexn = None
        #     self.E_H_field = None
        #     self.table_nod = None
        #     self.type_el = None
        #     self.mesh_xy = None
        #     self.n_msh_pts = None
        #     self.n_msh_el = None

        #if tstruc.plt_mesh:
        #    print("Suppressed inefficient matplotlib plotting of mesh...")
            #plotting.plot_msh(self.mesh_xy, prefix=tstruc.mesh_file, suffix='_EM')


### Calc unnormalised power in each EM mode Kokou equiv. of Eq. 8.
        try:
            print("  Calculating EM mode powers...")
            nnodes = 6
            if tstruc.inc_shape in tstruc.linear_element_shapes:
            ## Integration using analytically evaluated basis function integrals. Fast.
                self.EM_mode_power = NumBAT.em_mode_energy_int_v2_ez(
                    self.k_0, self.n_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod,
                    self.mesh_xy, self.Eig_values, self.sol1)
            else:
                if tstruc.inc_shape not in tstruc.curvilinear_element_shapes:
                    print("Warning: em_mode_energy_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
            # Integration by quadrature. Slowest.
                self.EM_mode_power = NumBAT.em_mode_energy_int_ez(
                    self.k_0, self.n_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod,
                    self.mesh_xy, self.Eig_values, self.sol1)
            # Bring Kokou's def into line with CW formulation.
            self.EM_mode_power = 2.0*self.EM_mode_power

        except KeyboardInterrupt:
            print("\n\n FEM routine EM_mode_energy_int interrupted by keyboard.\n\n")


### Calc energy (not power) in each EM mode - PRA Eq. 6.
        if self.calc_EM_mode_energy is True:
            print("Calculating EM mode energies...")
            try:
                nnodes = 6
                # import time
                # start = time.time()
                if tstruc.inc_shape in tstruc.linear_element_shapes:
                # # Semi-analytic integration. Fastest!
                # else:
                #     if tstruc.inc_shape not in tstruc.curvilinear_element_shapes:
                #         print("Warning: em_mode_e_energy_int - not sure if mesh contains curvi-linear elements", 
                #             "\n using slow quadrature integration by default.\n\n")
                # # Integration by quadrature. Slowest.
                    self.EM_mode_energy = NumBAT.em_mode_e_energy_int(
                        self.n_modes, self.n_msh_el, self.n_msh_pts, nnodes,
                        self.table_nod, self.type_el, tstruc.n_typ_el, self.v_refindexn,
                        self.mesh_xy, self.sol1)
                else:
                  print("\n\n FEM routine em_mode_e_energy_int is not implemented for this structure. Can't find group index. \n\n")
                  self.EM_mode_energy=np.zeros(self.n_modes, dtype=float)
                  
            except KeyboardInterrupt:
                print("\n\n FEM routine em_mode_e_energy_int interrupted by keyboard.\n\n")

        # This group velocity calc is not accurate in the presence of dispersion!
        # self.group_velocity_EM = self.EM_mode_power/self.EM_mode_power_energy

        # If considering a the backwards propagating Stokes field.
        if self.Stokes == True:
            self.Eig_values = -1*self.Eig_values
            self.sol1 = np.conj(self.sol1)

        ### Not necessary because EM FEM mesh always normalised in area to unity.
        # print area
        # x_tmp = []
        # y_tmp = []
        # for i in np.arange(self.n_msh_pts):
        #     x_tmp.append(self.mesh_xy[0,i])
        #     y_tmp.append(self.mesh_xy[1,i])
        # x_min = np.min(x_tmp); x_max=np.max(x_tmp)
        # y_min = np.min(y_tmp); y_max=np.max(y_tmp)
        # area = abs((x_max-x_min)*(y_max-y_min))
        # print area
        # self.EM_mode_power = self.EM_mode_power*area


#def calc_EM_mode_energy(self):  # these require extraction of numerical props from Fortran. Clean that up first.
#      assert(self.is_EM())
#      if not self.EM_mode_energy is None: return  # nothing to do

#    def calc_EM_mode_power(self):
#      assert(self.is_EM())
#      if not self.EM_mode_power is None: return  # nothing to do

#    def calc_EM_mode_power_and_energy(self):
#      assert(self.is_EM())
#      self.calc_EM_mode_power()
#      self.calc_EM_mode_energy()


    def calc_AC_modes(self, bcs=None):
        """ Run a Fortran FEM calculation to find the acoustic modes.

        Returns a ``Simulation`` object that has these key values:

        Eig_values: a 1d array of Eigenvalues (frequencies) in [1/s]

        sol1: the associated Eigenvectors, ie. the fields, stored as
               [field comp, node num on element, Eig value, el num]

        AC_mode_energy: the elastic power in the acoutic modes.
        """
        self.sim_type = SimType.AC

        self.d_in_m = self.structure.inc_a_x*1e-9

        # Build a table of materials only containing elastic properties referencing the original full list in Struct
        # Needed for restricting meshes to elastic only for example
        # el_conv_table maps the number of the acoustically active material in original material 
        # list to the position in the list of acoustically active materials
        # eg [vacuum, silicon, glass, vacuum, chalc] ->  {2:1,3:2,5:3}
        el_conv_table = {}
        i = 1; j = 1
        for mat in self.structure.d_materials.values():
            if mat.has_elastic_properties():
                el_conv_table[i] = j
                j += 1
            i += 1

        self.typ_el_AC = {}
        for k,v in el_conv_table.items():
          self.typ_el_AC[self.el_conv_table_n[k]] = v  # now keeps its own rather than take from EM_sim which might not exist

        #print('el_conv_table_n EM', self.el_conv_table_n, self.v_refindexn)
        #print('el_conv_table, AC', el_conv_table)
        #print('typ_el_AC', self.typ_el_AC)

        if self.n_modes < 20:
            self.n_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")

        # Parameters that control how FEM routine runs
        bnd_cdn_i = 0  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        if bcs == 'Open': 
            print('Attempting open elastic boundary conditions.')
            icond=1  # TODO: DO THIS ACTUILLY WORK?

        itermax = 30  # Maximum number of iterations for convergence
        AC_FEM_debug = 0  # Fortran routines will display & save add. info
        ARPACK_tol = 1e-10  # ARPACK accuracy (0.0 for machine precision)

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        if self.shift_Hz is None:
            # For AC problem shift is a frequency; [shift] = s^-1.
            v_list = []
            for el in range(self.structure.n_typ_el_AC):
                # Using acoustic velocity of longitudinal mode pg 215 Auld vol 1.
                v_list.append(np.sqrt(self.structure.c_tensor[0,0][el]/self.structure.rho[el]))
                # # Using acoustic velocity of shear mode pg 215 Auld vol 1.
                # v_list.append(np.sqrt(self.structure.c_tensor[3,3][el]/self.structure.rho[el]))
            AC_velocity = np.real(v_list).min()
            shift = np.real(AC_velocity*self.q_AC/(2.*np.pi))
            # print "shift", shift
            shift_nu = 0.9*shift
            # print "shift", shift
        else:
            shift_nu = self.shift_Hz

        # Take existing msh from EM FEM and manipulate mesh to exclude vacuum areas.
        if self.EM_sim:
            suplied_geo_flag = 1
            n_msh_el = self.EM_sim.n_msh_el
            n_msh_pts = self.EM_sim.n_msh_pts
            type_el = self.EM_sim.type_el
            type_nod = self.EM_sim.type_nod
            table_nod = self.EM_sim.table_nod
            mesh_xy = self.EM_sim.mesh_xy
            n_el_kept = 0
            n_msh_pts_AC = 0
            type_el_AC = []
            table_nod_AC_tmp = np.zeros(np.shape(table_nod))
            el_convert_tbl = {}
            el_convert_tbl_inv = {}
            node_convert_tbl = {}
            #if self.structure.plot_mesh: #TODO turn this back on
            #    plotting.plot_msh(mesh_xy, prefix=self.structure.mesh_file, suffix='_AC-orig')

            for el in range(n_msh_el):
                # print type_el[el]
                if type_el[el] in self.typ_el_AC:
                    # print "in", type_el[el]
                    type_el_AC.append(self.typ_el_AC[type_el[el]])
                    el_convert_tbl[n_el_kept] = el
                    el_convert_tbl_inv[el] = n_el_kept
                    for i in range(6):
                        # Leaves node numbering untouched
                        table_nod_AC_tmp[i][n_el_kept] = table_nod[i][el]
                    n_el_kept += 1
            n_msh_el_AC = n_el_kept
            # Find unique nodes
            node_lst_tmp = []
            for el in range(n_msh_el_AC):
                for i in range(6):
                    node_lst_tmp.append(table_nod_AC_tmp[i][el])
            unique_nodes = list(set(node_lst_tmp))
            n_msh_pts_AC = len(unique_nodes)
            unique_nodes = [int(j) for j in unique_nodes]
            # Mapping unique nodes to start from zero
            for i in range(n_msh_pts_AC):
                node_convert_tbl[unique_nodes[i]] = i
            # Creating finalised table_nod.
            table_nod_AC = []
            for i in range(6):
                el_tbl = []
                for el in range(n_msh_el_AC):
                    # Note table_nod needs to be adjust back to fortran indexing
                    el_tbl.append(node_convert_tbl[table_nod_AC_tmp[i][el]]+1)
                table_nod_AC.append(el_tbl)
            # Find the coordinates of chosen nodes.
            mesh_xy_AC = np.zeros((2,n_msh_pts_AC))
            for node in unique_nodes:
                # Note mesh_xy needs to be adjust back to fortran indexing
                mesh_xy_AC[0,node_convert_tbl[node]] = (mesh_xy[0,node-1])
                mesh_xy_AC[1,node_convert_tbl[node]] = (mesh_xy[1,node-1])

            self.el_convert_tbl = el_convert_tbl
            self.el_convert_tbl_inv = el_convert_tbl_inv
            self.node_convert_tbl = node_convert_tbl


            ### AC FEM uses Neumann B.C.s so type_nod is totally irrelevant!
            # # Find nodes on boundaries of materials
            # node_array = -1*np.ones(n_msh_pts)
            # interface_nodes = []
            # for el in range(n_msh_el):
            #     for i in range(6):
            #         node = table_nod[i][el]
            #         # Check if first time seen this node
            #         if node_array[node - 1] == -1: # adjust to python indexing
            #             node_array[node - 1] = type_el[el]
            #         else:
            #             if node_array[node - 1] != type_el[el]:
            #                 interface_nodes.append(node)
            # interface_nodes = list(set(interface_nodes))
            type_nod_AC = np.zeros(n_msh_pts_AC)

            # import matplotlib
            # matplotlib.use('pdf')
            # import matplotlib.pyplot as plt
            # plt.clf()
            # plt.figure(figsize=(13,13))
            # ax = plt.subplot(1,1,1)
            # for node in unique_nodes:
            #     if node in interface_nodes:
            #         type_nod_AC[node_convert_tbl[node]] = bnd_cdn_i
            #         plt.plot(mesh_xy_AC[0,node_convert_tbl[node]], mesh_xy_AC[1,node_convert_tbl[node]], 'ok')
            # ax.set_aspect('equal')
            # plt.savefig('boundary.pdf', bbox_inches='tight')
            self.n_msh_pts = n_msh_pts_AC
            self.n_msh_el = n_msh_el_AC
        # Default, indicates to use geometry subroutine in FEM routine.
        else: # No EM mesh data supplied
            suplied_geo_flag = 0
            with open(self.structure.mesh_file) as f:
                self.n_msh_pts, self.n_msh_el = [int(i) for i in f.readline().split()]
            table_nod_AC = np.zeros((6, self.n_msh_el))
            type_el_AC = np.zeros(self.n_msh_el)
            mesh_xy_AC = np.zeros((2,self.n_msh_pts))
            type_nod_AC = np.zeros(self.n_msh_pts)

        if AC_FEM_debug == 1:
            print('shift_nu', shift_nu)
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Output"):
                os.mkdir("Output")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")

        # Size of Fortran's complex superarray (scales with mesh)
        #int_max, cmplx_max, real_max = NumBAT.array_size(self.n_msh_el, self.n_modes)

        print('\n Structure has {0} mesh points and {1} mesh elements.'.format(self.n_msh_pts, self.n_msh_el))

        try:
            fort_err = 0

            tstruc=self.structure
            show_mem_est = False

            resm = NumBAT.calc_ac_modes(
                self.q_AC, self.n_modes,
                AC_FEM_debug, show_mem_est, tstruc.mesh_file, self.n_msh_pts,
                self.n_msh_el, tstruc.n_typ_el_AC, tstruc.c_tensor, tstruc.rho,
                self.d_in_m, shift_nu, bnd_cdn_i, itermax, ARPACK_tol,
                tstruc.plotting_fields,
                suplied_geo_flag, type_nod_AC, 
                tstruc.symmetry_flag, table_nod_AC, type_el_AC, mesh_xy_AC)

            table_nod_out, type_el_out, mesh_xy_out, \
            self.Eig_values, self.sol1, self.mode_pol, fort_err, fort_mesg = resm

            if fort_err != 0:
                fort_mesg = str(fort_mesg, 'utf-8')  # fort_mesg comes back as a byte string.
                report_and_exit('Fortran error in solving for acoustic modes: \n'
                                ' NumBAT Fortran error code = %d. \n Message: \n %s'%(fort_err, fort_mesg))

            # FEM Eigenvalue is frequency, rather than angular frequency Omega
            self.Omega_AC = self.Eig_values*2*np.pi

        except KeyboardInterrupt:
            print("\n\n FEM routine calc_AC_modes",\
            "interrupted by keyboard.\n\n")

        # Retrieve the material properties of each mesh point.
        self.ls_material = NumBAT.array_material_ac(
                self.n_msh_el, self.structure.n_typ_el_AC, type_el_AC,
             self.structure.rho, self.structure.c_tensor, 
             self.structure.p_tensor, self.structure.eta_tensor)

        #if self.structure.plt_mesh:  #TODO re-enable
        #    plotting.plot_msh(mesh_xy_AC, prefix=self.structure.mesh_file, suffix='_AC-in')
        #    plotting.plot_msh(mesh_xy_out, prefix=self.structure.mesh_file, suffix='_AC-out')

        # if self.EM_sim is None:
        #     table_nod_out = None
        #     type_el_out = None
        #     mesh_xy_out = None
        #     self.table_nod = table_nod_AC
        #     self.type_el = type_el_AC
        #     self.mesh_xy = mesh_xy_AC
        # else:
        self.table_nod = table_nod_out
        self.type_el = type_el_out
        self.mesh_xy = mesh_xy_out

### Calc unnormalised power in each AC mode - PRA Eq. 18.
        if self.calc_AC_mode_power is True:
            try:
                nnodes = 6
                if self.structure.inc_shape in self.structure.linear_element_shapes:
                # Semi-analytic integration following KD 9/9/16 notes. Fastest!
                    self.AC_mode_power = NumBAT.ac_mode_power_int_v4(
                        self.n_modes, self.n_msh_el, self.n_msh_pts,
                        nnodes, self.table_nod, self.type_el, self.mesh_xy,
                        self.structure.n_typ_el_AC, self.structure.c_tensor,
                        self.q_AC, self.Omega_AC, self.sol1)
                else:
                    if self.structure.inc_shape not in self.structure.curvilinear_element_shapes:
                        print("Warning: ac_mode_power_int - not sure if mesh contains curvi-linear elements", 
                            "\n using slow quadrature integration by default.\n\n")
                # Integration by quadrature. Slowest.
                    self.AC_mode_power = NumBAT.ac_mode_power_int(
                        self.n_modes, self.n_msh_el, self.n_msh_pts,
                        nnodes, self.table_nod, self.type_el, self.mesh_xy,
                        self.structure.n_typ_el_AC, self.structure.c_tensor_z,
                        self.q_AC, self.Omega_AC, self.sol1, AC_FEM_debug)
            except KeyboardInterrupt:
                print("\n\n FEM routine AC_mode_energy_int interrupted by keyboard.\n\n")


### Calc unnormalised elastic energy in each AC mode - PRA Eq. 16.
        try:
            nnodes = 6
            if self.structure.inc_shape in self.structure.linear_element_shapes:
            # Semi-analytic integration. Fastest!
                self.AC_mode_energy= NumBAT.ac_mode_elastic_energy_int_v4(
                    self.n_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod, self.type_el, self.mesh_xy,
                    self.structure.n_typ_el_AC, self.structure.rho,
                    self.Omega_AC, self.sol1)
            else:
                if self.structure.inc_shape not in self.structure.curvilinear_element_shapes:
                    print("Warning: ac_mode_elastic_energy_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
            # Integration by quadrature. Slowest.
                self.AC_mode_energy= NumBAT.ac_mode_elastic_energy_int(
                    self.n_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod, self.type_el, self.mesh_xy,
                    self.structure.n_typ_el_AC, self.structure.rho,
                    self.Omega_AC, self.sol1, AC_FEM_debug)
        except KeyboardInterrupt:
            print("\n\n FEM routine AC_mode_elastic_energy_int interrupted by keyboard.\n\n")


def bkwd_Stokes_modes(EM_sim):
    """ Defines the backward travelling Stokes waves as the conjugate
        of the forward travelling pump waves.

    Returns a ``Simulation`` object that has these key values:

    Eig_values: a 1d array of Eigenvalues (propagation constants) in [1/m]

    sol1: the associated Eigenvectors, ie. the fields, stored as
           [field comp, node nu on element, Eig value, el nu]

    EM_mode_power: the power in the Stokes modes. Note this power is negative because the modes 
                   are travelling in the negative z-direction.
    """

    EM_sim.clean_for_save()

    Stokes_modes = copy.deepcopy(EM_sim)
    Stokes_modes.sol1 = np.conj(Stokes_modes.sol1)
    Stokes_modes.Eig_values = -1.0*Stokes_modes.Eig_values
    Stokes_modes.EM_mode_power = -1.0*Stokes_modes.EM_mode_power
    return Stokes_modes



def fwd_Stokes_modes(EM_sim):
    """ Defines the forward travelling Stokes waves as a copy
        of the forward travelling pump waves.

    Returns a ``Simulation`` object that has these key values:

    """

    EM_sim.clean_for_save()

    Stokes_modes = copy.deepcopy(EM_sim)
    return Stokes_modes
