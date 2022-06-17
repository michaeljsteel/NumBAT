# materials.py is a subroutine of NumBAT that defines Material objects,
# these represent dispersive lossy refractive indices and possess
# methods to interpolate n from tabulated data.

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


import os
import numpy as np
import sys
import traceback
from scipy.interpolate import interp1d
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
from math import sqrt

import json
import re

from nbtypes import *
from reporting import *



class VoigtTensor4(object):
  def __init__(self, sym, src_dict, src_file):
    self.mat=np.zeros([7,7],dtype=float)  # unit indexing
    self.sym = sym  # eg 'c', 'p', 'eta'
    self.d=src_dict
    self.json_file=src_file

  def read(self, m, n, optional=False):
    elt='{0}_{1}{2}'.format(self.sym, m, n)

    if elt not in self.d: 
      if not optional:
        report_and_exit('Failed to read required tensor element {0} in material data file {1}'.format(
          elt, self.json_file))
      else:
        return False

    self.mat[m,n]=self.d[elt]
    return True

  def load_isotropic(self):
    self.read(1,1) 
    self.read(1,2) 
    self.read(4,4)
    self.mat[2,2]= self.mat[1,1]
    self.mat[3,3]= self.mat[1,1]
    self.mat[5,5]= self.mat[4,4]
    self.mat[6,6]= self.mat[4,4]
    self.mat[2,1]= self.mat[1,2]
    self.mat[2,3]= self.mat[1,2]
    self.mat[1,3]= self.mat[1,2]
    self.mat[3,1]= self.mat[1,2]
    self.mat[3,2]= self.mat[1,2]

#  def copy(self, tm, tn, fm, fn):
#    self.mat[tm,tn]=self.mat[fm,fn]

#  def copy(self, tm, tn, fm, fn):
#    self.mat[tm,tn]=self.mat[fm,fn]

#def value(self, m, n): return self.mat[m-1, n-1]

  def __getitem__(self, k):
    return self.mat[k[0], k[1]]

  def __setitem__(self, k, v):
    self.mat[k[0], k[1]]=v

  def rotate(self, theta, rotation_axis):
    rmat=rotate_tensor(self.mat[1:,1:], theta, rotation_axis)
    self.mat[1:,1:]=rmat

  def __str__(self):
    s='\nVoigt tensor {0}:\n'.format(self.sym)
    s+=str(self.mat[1:,1:])
    return s

  def dump(self):
    print('\nVoigt tensor {0}'.format(self.sym))
    print(self.mat[1:,1:])

class Material(object):
    """Class representing a waveguide material. 
      
      Materials include the following properties and corresponding units:
          -  Refractive index []
          -  Density [kg/m3]
          -  Stiffness tensor component [Pa]
          -  Photoelastic tensor component []
          -  Acoustic loss tensor component [Pa s]

    """

# class level variables
    _data_loc=''
    _materials={}

    @classmethod
    def _set_file_locations(cls):
      this_dir= os.path.dirname(os.path.realpath(__file__))
      Material._data_loc= os.path.join(this_dir, "material_data", "")

    @classmethod
    def _get_material(cls, s):
      if not len(Material._materials):
        Material._set_file_locations()
        for f in os.listdir(Material._data_loc):
          if f.endswith(".json"):
            Material._materials[f[:-5]] = Material(f)
      try:
        mat= Material._materials[s]
      except KeyError:
        report_and_exit('Material data {0} file not found.'.format(s))

      return mat 

    def __init__(self, data_file):

        if not len(Material._data_loc): self.__class__._set_file_locations()

        self.json_file=data_file
        try:
            self.load_data_file(Material._data_loc, self.json_file)
        except FileNotFoundError:
            report_and_exit('Material data {0} file not found.'.format(self.json_file))



    def __str__(self):
      s='''
          Material: {0}
          File: {1}
          Source: {2}
          Date: {3}
          '''.format(
            self.chemical,
            self.file_name,
            self.author,
            self.date)
      if len(self.comment): s+='Comment: '+self.comment
      return s

    def full_str(self):
      s=self.__str__()
      s+= str(self.c_tensor)
      s+= str(self.eta_tensor)
      s+= str(self.p_tensor)
      return s

    def vac_longitudinal(self):
      assert(not self.anisotropic)
      # lame lambda = c_12
      # lame mu = c_44
      #  v = sqrt(c_11/rho)
      #    =sqrt((c12 + 2c44)/rho)
      #    =sqrt((lambda + 2mu)/rho)
      return sqrt(self.c_tensor[1,1]/self.s)


    def vac_shear(self):
      assert(not self.anisotropic)
      #  v = sqrt(c_44/rho)
      #    =sqrt((mu)/rho)
      return sqrt(self.c_tensor[4,4]/self.s)


    def has_elastic_properties(self):
        return self.s is not None

    def load_data_file(self, dataloc, data_file, alt_path=''):  
        """
        Load data from json file.
        
        Args:
            data_file  (str): name of data file located in NumBAT/backend/material_data
            
            alt_path  (str): non standard path to data_file
        
        """

        with open(dataloc+data_file,'r') as fin:
            s_in = ''.join(fin.readlines())
            s_in = re.sub(r'//.*\n','\n', s_in)

            try:
              self._params = json.loads(s_in)
            except Exception as err:
              traceback.print_exc()
              report_and_exit('JSON parsing error: {0} for file {1}'.format(err, self.json_file))

            self.file_name = self._params['file_name']  # Name of this file, will be used as identifier
            self.chemical = self._params['chemical']  # Chemical composition
            self.author = self._params['author']  # Author of data
            self.date = self._params['date']  # Year of data publication/measurement
            self.institution = self._params['institution']  # Source institution
            self.doi = self._params['doi']  # doi or, failing that, the http address

            self.comment = self._params.get('comment', '')  # general comment for any purpose

            Re_n = self._params['Re_n']  # Real part of refractive index []
            Im_n = self._params['Im_n']  # Imaginary part of refractive index []
            self.n = (Re_n + 1j*Im_n)  # Complex refractive index []
            self.s = self._params['s']  # Density [kg/m3]

#self.c_11 = self._params['c_11']  # Stiffness tensor component [Pa]
#            self.c_12 = self._params['c_12']  # Stiffness tensor component [Pa]
#            self.c_44 = self._params['c_44']  # Stiffness tensor component [Pa]
#            self.p_11 = self._params['p_11']  # Photoelastic tensor component []
#            self.p_12 = self._params['p_12']  # Photoelastic tensor component []
#            self.p_44 = self._params['p_44']  # Photoelastic tensor component []
#            self.eta_11 = self._params['eta_11']  # Acoustic loss tensor component [Pa s]
#            self.eta_12 = self._params['eta_12']  # Acoustic loss tensor component [Pa s]
#            self.eta_44 = self._params['eta_44']  # Acoustic loss tensor component [Pa s]

            self.crystal = CrystalGroup.Isotropic

            self.c_tensor=VoigtTensor4('c', self._params, self.json_file)
            self.eta_tensor=VoigtTensor4('eta', self._params, self.json_file)
            self.p_tensor=VoigtTensor4('p', self._params, self.json_file)

            self.c_tensor.load_isotropic()
            self.p_tensor.load_isotropic()
            self.eta_tensor.load_isotropic()

            self.load_tensors()


    def load_cubic_crystal(self): #(don't really need this as isotropic materials are the same)

      try:
        self.c_tensor.read(1,1)
        self.c_tensor.read(1,2)
        self.c_tensor[1,3] = self.c_tensor[1,2]
        self.c_tensor[2,1] = self.c_tensor[1,2]
        self.c_tensor[2,2] = self.c_tensor[1,1]
        self.c_tensor[2,3] = self.c_tensor[1,2]
        self.c_tensor[3,1] = self.c_tensor[1,2]
        self.c_tensor[3,2] = self.c_tensor[1,2]
        self.c_tensor[3,3] = self.c_tensor[1,1]
        self.c_tensor.read(4,4)
        self.c_tensor[5,5] = self.c_tensor[4,4]
        self.c_tensor[6,6] = self.c_tensor[4,4]

        self.eta_tensor.read(1,1)
        self.eta_tensor.read(1,2)
        self.eta_tensor[1,3] = self.eta_tensor[1,2]
        self.eta_tensor[2,1] = self.eta_tensor[1,2]
        self.eta_tensor[2,2] = self.eta_tensor[1,1]
        self.eta_tensor[2,3] = self.eta_tensor[1,2]
        self.eta_tensor[3,1] = self.eta_tensor[1,2]
        self.eta_tensor[3,2] = self.eta_tensor[1,2]
        self.eta_tensor[3,3] = self.eta_tensor[1,1]
        self.eta_tensor.read(4,4)
        self.eta_tensor[5,5] = self.eta_tensor[4,4]
        self.eta_tensor[6,6] = self.eta_tensor[4,4]

        self.p_tensor.read(1,1)
        self.p_tensor.read(1,2)

        self.p_tensor[1,3] = self.p_tensor[1,2]
        self.p_tensor[2,1] = self.p_tensor[1,2]
        self.p_tensor[2,2] = self.p_tensor[1,1]
        self.p_tensor[2,3] = self.p_tensor[1,2]
        self.p_tensor[3,1] = self.p_tensor[1,2]
        self.p_tensor[3,2] = self.p_tensor[1,2]
        self.p_tensor[3,3] = self.p_tensor[1,1]
        self.p_tensor.read(4,4)

        # According to Powell, for Oh group, these are distinct elements, but no one seems to quote them
        if not self.p_tensor.read(5,5, optional=True): self.p_tensor[5,5] = self.p_tensor[4,4]
        if not self.p_tensor.read(6,6, optional=True): self.p_tensor[6,6] = self.p_tensor[4,4]

      except:
        report_and_exit('Failed to load cubic crystal class in material data file {0}'.format(self.json_file))

    def load_trigonal_crystal(self):
      try:
        for (i,j) in [(1,1), (1,2), (1,3), (1,4), (3,3), (4,4), (6,6)]:
          self.c_tensor.read(i,j)

        self.c_tensor[2,1] = self.c_tensor[1,2]
        self.c_tensor[2,2] = self.c_tensor[1,1]
        self.c_tensor[2,3] = self.c_tensor[1,3]
        self.c_tensor[2,4] = -self.c_tensor[1,4]

        self.c_tensor[3,1] = self.c_tensor[1,3]
        self.c_tensor[3,2] = self.c_tensor[1,3]
        self.c_tensor[3,3] = self.c_tensor[3,3]

        self.c_tensor[4,1] = self.c_tensor[1,4]
        self.c_tensor[4,2] = -self.c_tensor[1,4]

        self.c_tensor[5,5] = self.c_tensor[4,4]
        self.c_tensor[5,6] = self.c_tensor[1,4]
        self.c_tensor[6,5] = self.c_tensor[1,4]

        for (i,j) in [(1,1), (1,2), (1,3), (1,4), (3,3), (4,4), (6,6)]:
          self.eta_tensor.read(i,j)

        self.eta_tensor[2,1] = self.eta_tensor[1,2]
        self.eta_tensor[2,2] = self.eta_tensor[1,1]
        self.eta_tensor[2,3] = self.eta_tensor[1,3]
        self.eta_tensor[2,4] = -self.eta_tensor[1,4]

        self.eta_tensor[3,1] = self.eta_tensor[1,3]
        self.eta_tensor[3,2] = self.eta_tensor[1,3]
        self.eta_tensor[3,3] = self.eta_tensor[3,3]

        self.eta_tensor[4,1] = self.eta_tensor[1,4]
        self.eta_tensor[4,2] = -self.eta_tensor[1,4]

        self.eta_tensor[5,5] = self.eta_tensor[4,4]
        self.eta_tensor[5,6] = self.eta_tensor[1,4]
        self.eta_tensor[6,5] = self.eta_tensor[1,4]

        #TODO: confirm correct symmetry properties for p. Using trigonal = C3v from Powell
        self.p_tensor.read(1,1)
        self.p_tensor.read(1,2)
        self.p_tensor.read(1,3)
        self.p_tensor.read(1,4)

        self.p_tensor[2,1] = self.p_tensor[1,2]
        self.p_tensor[2,2] = self.p_tensor[1,1]
        self.p_tensor[2,3] = self.p_tensor[1,3]
        self.p_tensor[2,4] = -self.p_tensor[1,4]

        self.p_tensor.read(3,1)
        self.p_tensor[3,2] = self.p_tensor[3,1]
        self.p_tensor.read(3,3)

        self.p_tensor.read(4,1)
        self.p_tensor[4,2] = -self.p_tensor[4,1]
        self.p_tensor.read(4,4)

        self.p_tensor[5,5] = -self.p_tensor[4,4]
        self.p_tensor[5,6] = 2*self.p_tensor[4,1]
        self.p_tensor[6,5] = -self.p_tensor[1,4]
        self.p_tensor[6,6] = self.p_tensor[1,1]- self.p_tensor[1,2]

      except:
        report_and_exit('Failed to load trigonal crystal class in material data file {0}'.format(self.json_file))
          
    def load_general_crystal(self):
      try:  # full anisotropic tensor components
          for i in range(1,7):
            for j in range(1,7):
              self.c_tensor.read(i,j)
              self.p_tensor.read(i,j)
              self.eta_tensor.read(i,j)

      except KeyError:
        report_and_exit('Failed to load anisotropic crystal class in material data file {0}'.format(self.json_file))
    


    def set_refractive_index(self, nr, ni=0.0):
      self.n = nr + 1j*ni

    def load_tensors(self): # not do this unless symmetry is off?

      self.anisotropic = False

      self.crystal=CrystalGroup.Unknown

      if not 'crystal_class' in self._params: return

      try:
        self.crystal = CrystalGroup[self._params['crystal_class']]
      except ValueError:
        print('Unknown crystal class in material data file')
        sys.exit(1)


      if self.crystal==CrystalGroup.Isotropic: return

      self.anisotropic = True

      if self.crystal==CrystalGroup.Trigonal:
        self.load_trigonal_crystal()
        return
      elif self.crystal==CrystalGroup.Cubic:
        self.load_cubic_crystal()
        return
        
      elif self.crystal==CrystalGroup.GeneralAnisotropic:
        self.load_general_crystal()
        return




    def rotate_axis(self, theta, rotate_axis, save_rotated_tensors=False):
        """ Rotate crystal axis by theta radians.

            Args:
                theta  (float): Angle to rotate by in radians.

                rotate_axis  (str): Axis around which to rotate.

            Keyword Args:
                save_rotated_tensors  (bool): Save rotated tensors to csv.

            Returns:
                ``Material`` object with rotated tensor values.
        """
     
        self.c_tensor.rotate(theta, rotate_axis)
        self.p_tensor.rotate(theta, rotate_axis)
        self.eta_tensor.rotate(theta, rotate_axis)
        if save_rotated_tensors:
            np.savetxt('rotated_c_tensor.csv', self.c_tensor.mat, delimiter=',')
            np.savetxt('rotated_p_tensor.csv', self.p_tensor.mat, delimiter=',')
            np.savetxt('rotated_eta_tensor.csv', self.eta_tensor.mat, delimiter=',')


# Array that converts between 4th rank tensors in terms of x,y,z and Voigt notation
#               [[xx,xy,xz], [yx,yy,yz], [zx,zy,zz]]
to_Voigt = np.array([[0,5,4], [5,1,3], [4,3,2]]) 



def rotation_matrix_sum(i, j, k, l, tensor_orig, mat_R):
    """
    Inner loop of rotation matrix summation.
    """
    tensor_prime_comp = 0
    for q in range(3):
        for r in range(3):
            V1 = to_Voigt[q,r]
            for s in range(3):
                for t in range(3):
                    V2 = to_Voigt[s,t]
                    tensor_prime_comp += mat_R[i,q] * mat_R[j,r] * mat_R[k,s] * mat_R[l,t] * tensor_orig[V1,V2]

    return tensor_prime_comp


def rotate_tensor(tensor_orig, theta, rotation_axis):
    """
    Rotate all acoustic material tensor by theta radians around chosen
    rotation_axis.

    Args:
        tensor_orig  (array): Tensor to be rotated.

        theta  (float): Angle to rotate by in radians.

        rotation_axis  (str): Axis around which to rotate.
    """
    if rotation_axis == 'x-axis':
        mat_R = np.array([[1,0,0], [0,np.cos(theta),-np.sin(theta)], [0,np.sin(theta),np.cos(theta)]])
    if rotation_axis == 'y-axis':
        mat_R = np.array([[np.cos(theta),0,np.sin(theta)], [0,1,0], [-np.sin(theta),0,np.cos(theta)]])
    if rotation_axis == 'z-axis':
        mat_R = np.array([[np.cos(theta),-np.sin(theta),0], [np.sin(theta),np.cos(theta),0], [0,0,1]])

    tensor_prime = np.zeros((6,6))
    for i in range(3):
        for j in range(3):
            V1 = to_Voigt[i,j]
            for k in range(3):
                for l in range(3):
                    V2 = to_Voigt[k,l]
                    tensor_prime[V1,V2] = rotation_matrix_sum(i,j,k,l,tensor_orig,mat_R)

    return tensor_prime


def isotropic_stiffness(E, v):
    """
    Calculate the stiffness matrix components of isotropic 
    materials, given the two free parameters.

    Ref: www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm

    Args:
        E  (float): Youngs modulus

        v  (float): Poisson ratio
    """
    c_11 = E*(1-v)/((1+v)*(1-2*v))
    c_12 = E*(v)/((1+v)*(1-2*v))
    c_44 = (E*(1-2*v)/((1+v)*(1-2*v)))/2

    return c_11, c_12, c_44


#g_materials={}


def get_material(s):
  return Material._get_material(s)
#  global g_materials
#  if not len(g_materials):
#this_dir= os.path.dirname(os.path.realpath(__file__))
#    data_loc= os.path.join(this_dir, "material_data", "")
#  for f in os.listdir(data_loc):
#    if f.endswith(".json"):
##g_materials[f[:-5]] = Material(data_loc, f[:-5])
#        print('loading mat', f)
#        g_materials[f[:-5]] = Material(f)

#  return g_materials[s]

# This code is deprecated and will be removed. Use get_material() instead.
#print('loading materials dict')
#materials_dict = {}
#this_directory = os.path.dirname(os.path.realpath(__file__))
#data_location = os.path.join(this_directory, "material_data", "")
#for file in os.listdir(data_location):
#    if file.endswith(".json"):
#        materials_dict[file[:-5]] = Material(data_location, file[:-5])
##print('cone loading materials dict')
