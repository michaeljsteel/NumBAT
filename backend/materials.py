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
import matplotlib.pyplot as plt
from math import sqrt

import json
import re

from nbtypes import *
from reporting import *

class BadMaterialFileError(Exception):
    pass


class VoigtTensor4(object):
  '''A class for representing rank 4 tensors in the compact Voigt representation.'''

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
    rmat = _rotate_Voigt_tensor(self.mat[1:,1:], theta, rotation_axis)
    self.mat[1:,1:] = rmat

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
    def _make_material(cls, s):
      if not len(Material._materials):  # first time through, we load all the materials
        Material._set_file_locations()
        for f in os.listdir(Material._data_loc):
            if f.endswith(".json"):
                try:
                    new_mat = Material(f)
                except BadMaterialFileError as err:
                    report_and_exit(str(err))

                if new_mat.material_name in Material._materials: # two mats have the same name  TODO: file_name is actually a mat_name
                    raise report_and_exit(f"Material file {f} has the same name as an existing material {new_mat.material_name}.")
            #Material._materials[f[:-5]] 
                Material._materials[new_mat.material_name] = new_mat

      try:
        mat= Material._materials[s]
      except KeyError:
        report_and_exit(f'Material {s} not found in material_data folder.')

      return mat 

    def __init__(self, data_file):

        if not len(Material._data_loc): self.__class__._set_file_locations()

        self.json_file=data_file
        try:
            self._load_data_file(Material._data_loc, self.json_file)
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
            self.material_name,
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

    def elastic_properties(self):
        '''Returns a string containing key elastic properties of the material.'''
        try:
            s =  'Material:       {0}'.format(self.material_name)
            s+='\nDensity:        {0:.3f}'.format(self.rho)
            s+='\nVelocity long.: {0:.3f}'.format( self.Vac_longitudinal() )
            s+='\nVelocity shear: {0:.3f}'.format( self.Vac_shear())
        except:
            s='Unknown/undefined elastic parameters in material '+self.material_name
        return s

    def Vac_longitudinal(self):
      '''For an isotropic material, returns the longitudinal (P-wave) elastic phase velocity.'''
      assert(not self.anisotropic)
      # lame lambda = c_12
      # lame mu = c_44
      #  v = sqrt(c_11/rho)
      #    =sqrt((c12 + 2c44)/rho)
      #    =sqrt((lambda + 2mu)/rho)
      if not self.rho or self.rho == 0: # Catch vacuum cases
          return 0.
      else: 
          return sqrt(self.c_tensor[1,1]/self.rho)


    def Vac_shear(self):
      '''For an isotropic material, returns the shear (S-wave) elastic phase velocity.'''
      assert(not self.anisotropic)
      #  v = sqrt(c_44/rho)
      #    =sqrt((mu)/rho)
      if not self.rho or self.rho == 0: # Catch vacuum cases
          return 0.
      else: 
          return sqrt(self.c_tensor[4,4]/self.rho)


    def has_elastic_properties(self):
        '''Returns true if the material has at least some elastic properties defined.'''
        return self.rho is not None

    def _load_data_file(self, dataloc, data_file, alt_path=''):  
        """
        Load material data from json file.
        
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

            self.material_name = self._params.get('material_name', 'NOFILENAME')  # Name of this file, will be used as identifier and must be present
            if self.material_name == 'NOFILENAME':
                raise BadMaterialFileError(f"Material file {data_file} has no 'material_name' field.")

            self.format = self._params.get('format', 'NOFORMAT')  
            if self.format == 'NOFORMAT':
                raise BadMaterialFileError(f"Material file {data_file} has no 'format' field.")

            if self.format != 'NumBATMaterial-fmt-2.0':
                raise BadMaterialFileError(f"Material file {data_file} must be in format 'NumBATMaterial-Fmt-2.0'.")


            self.chemical = self._params['chemical']  # Chemical composition
            self.author = self._params['author']  # Author of data
            self.date = self._params['date']  # Year of data publication/measurement
            self.institution = self._params['institution']  # Source institution
            self.doi = self._params['doi']  # doi or, failing that, the http address

            self.comment = self._params.get('comment', '')  # general comment for any purpose

            Re_n = self._params['Re_n']  # Real part of refractive index []
            Im_n = self._params['Im_n']  # Imaginary part of refractive index []
            self.refindex_n = (Re_n + 1j*Im_n)  # Complex refractive index []
            self.rho = self._params['s']  # Density [kg/m3]

#self.c_11 = self._params['c_11']  # Stiffness tensor component [Pa]
#            self.c_12 = self._params['c_12']  # Stiffness tensor component [Pa]
#            self.c_44 = self._params['c_44']  # Stiffness tensor component [Pa]
#            self.p_11 = self._params['p_11']  # Photoelastic tensor component []
#            self.p_12 = self._params['p_12']  # Photoelastic tensor component []
#            self.p_44 = self._params['p_44']  # Photoelastic tensor component []
#            self.eta_11 = self._params['eta_11']  # Acoustic loss tensor component [Pa s]
#            self.eta_12 = self._params['eta_12']  # Acoustic loss tensor component [Pa s]
#            self.eta_44 = self._params['eta_44']  # Acoustic loss tensor component [Pa s]


            if not 'crystal_class' in self._params:
                raise BadMaterialFileError(f"Material file {data_file} has no 'crystal_class' field.")
            try:
                self.crystal = CrystalGroup[self._params['crystal_class']]
            except ValueError:
                print('Unknown crystal class in material data file')
                raise BadMaterialFileError(f"Unknown crystal class in material data file {data_file}")

            self.crystal = CrystalGroup.Isotropic

            self.c_tensor=VoigtTensor4('c', self._params, self.json_file)
            self.eta_tensor=VoigtTensor4('eta', self._params, self.json_file)
            self.p_tensor=VoigtTensor4('p', self._params, self.json_file)

            self.c_tensor.load_isotropic()
            self.p_tensor.load_isotropic()
            self.eta_tensor.load_isotropic()

            self.load_tensors()
            


    def is_vacuum(self): 
        '''Returns True if the material is the vacuum.'''
        return self.chemical == 'Vacuum'

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
      self.refindex_n = nr + 1j*ni

    def load_tensors(self): # not do this unless symmetry is off?

      self.anisotropic = False

      self.crystal=CrystalGroup.Unknown
      
      if not 'crystal_class' in self._params: return

      try:
        self.crystal = CrystalGroup[self._params['crystal_class']]
      except ValueError:
        print('Unknown crystal class in material data file')
        sys.exit(1)  #TODO: exit properly


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



def rotate_tensor_elt(i, j, k, l, T_pqrs, mat_R):
    '''
    Calculates the element ijkl of the rotated tensor T' from the original 
    rank-4 tensor T_PQ in 6x6 Voigt notation under the rotation specified by the 3x3 matrix R.
    '''

    Tp_ijkl = 0

    for q in range(3):
        for r in range(3):
            V1 = to_Voigt[q,r]
            for s in range(3):
                for t in range(3):
                    V2 = to_Voigt[s,t]
                    Tp_ijkl += mat_R[i,q] * mat_R[j,r] * mat_R[k,s] * mat_R[l,t] * T_pqrs[V1,V2]

    return Tp_ijkl 


def _rotate_Voigt_tensor(T_PQ, theta, rotation_axis):
    """
    Rotate an acoustic material tensor by theta radians around a specified rotation_axis.
    T_PQ is a rank-4 tensor expressed in 6x6 Voigt notation.

    The complete operation in 3x3x3x3 notation is 
    T'_ijkl  = sum_{pqrs} R_ip R_jq R_kr R_ls  T_pqrs.

    The result T'_ijkl is returned in Voigt format T'_PQ.

    Args:
        T_PQ  (array): Tensor to be rotated.

        theta  (float): Angle to rotate by in radians.

        rotation_axis  (str): Axis around which to rotate.
    """

    if type(rotation_axis) == type('axis'):
        if rotation_axis == 'x-axis':
            mat_R = np.array([[1,0,0], [0,np.cos(theta),-np.sin(theta)], [0,np.sin(theta),np.cos(theta)]])
        elif rotation_axis == 'y-axis':
            mat_R = np.array([[np.cos(theta),0,np.sin(theta)], [0,1,0], [-np.sin(theta),0,np.cos(theta)]])
        elif rotation_axis == 'z-axis':
            mat_R = np.array([[np.cos(theta),-np.sin(theta),0], [np.sin(theta),np.cos(theta),0], [0,0,1]])
    else:
        emsg = f"Can't convert {rotation_axis} to a 3-element unit vector."
        try:
            if isinstance(rotation_axis, list): # try to convert to numpy
                uvec = np.array(rotation_axis)
            elif isinstance(rotation_axis, np.ndarray): 
                uvec = rotation_axis
            else:
                report_and_exit(emsg)
        except:
            report_and_exit(emsg)

        # uvec is now a numpy array of some length
        if len(uvec) != 3 or np.abs(uvec) == 0.0 :    
            report_and_exit(emsg)

        # normlise u 
        uvec = uvec / np.abs(uvec)
        ct = math.cos(theta)
        st = math.sin(theta)
        omct = 1-ct
        ux,uy,uz = uvec[:]

        mat_R =  np.array([
            [ct + ux**2*omct,  ux*uy*omct-uz*st, ux*uz*omct+uy*st ],
            [uy*ux*omct+uz*st, ct + uy**2*omct,  uy*uz*omct-ux*st ],
            [uz*ux*omct-uy*st, uz*uy*omct+ux*st, ct+uz**2*omct]
            ])


    Tp_PQ = np.zeros((6,6))
    for i in range(3):
        for j in range(3):
            V1 = to_Voigt[i,j]
            for k in range(3):
                for l in range(3):
                    V2 = to_Voigt[k,l]
                    Tp_PQ[V1,V2] = rotate_tensor_elt(i,j,k,l, T_PQ, mat_R)

    return Tp_PQ


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


def make_material(s): 
    return Material._make_material(s)

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
