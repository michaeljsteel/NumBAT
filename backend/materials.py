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


import sys
import os
import traceback
import math
import json
import re

import numpy as np
import numpy.linalg
import tempfile
import subprocess
import scipy.linalg

import numbattools

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors


from nbtypes import CrystalGroup
from reporting import report_and_exit


class BadMaterialFileError(Exception):
    pass


# TODO: this class should not know about files. Pull the reading from files out.


# TODO: VoigtTensors check for symmetry consistency.  Seems to be not quite symmetric, 
# a mix of symmetric and skew-symmetric. Is that right?

class VoigtTensor4(object):
    '''A class for representing rank 4 tensors in the compact Voigt representation.'''

    def __init__(self, sym, src_dict=None, src_file=None):
        self.mat = np.zeros([7, 7], dtype=float)  # unit indexing
        self.sym = sym  # eg 'c', 'p', 'eta'

        self.d = src_dict
        self.json_file = src_file

    def as_zerobase_matrix(self):
        return self.mat[1:, 1:]
    
    def read(self, m, n, optional=False):
        elt = f'{self.sym}_{m}{n}'

        if self.d is None or self.json_file is None:
            report_and_exit('Voigt Tensor' + self.sym +
                            'is not setup for file reading.')

        if elt not in self.d:
            if not optional:
                report_and_exit(
                    f'Failed to read required tensor element {elt} in material data file {self.json_file}')
            else:
                return False

        self.mat[m, n] = self.d[elt]
        return True

    def load_isotropic(self):
        self.read(1, 1)
        self.read(1, 2)
        self.read(4, 4)
        self.set_isotropic(self.mat[1, 1], self.mat[1, 2], self.mat[4, 4])

    #  def copy(self, tm, tn, fm, fn):
    #    self.mat[tm,tn]=self.mat[fm,fn]

    #  def copy(self, tm, tn, fm, fn):
    #    self.mat[tm,tn]=self.mat[fm,fn]

    # def value(self, m, n): return self.mat[m-1, n-1]
    def set_isotropic(self, m11, m12, m44):
        self.mat[1, 1] = m11
        self.mat[1, 2] = m12
        self.mat[4, 4] = m44

        self.mat[2, 2] = self.mat[1, 1]
        self.mat[3, 3] = self.mat[1, 1]
        self.mat[5, 5] = self.mat[4, 4]
        self.mat[6, 6] = self.mat[4, 4]
        self.mat[2, 1] = self.mat[1, 2]
        self.mat[2, 3] = self.mat[1, 2]
        self.mat[1, 3] = self.mat[1, 2]
        self.mat[3, 1] = self.mat[1, 2]
        self.mat[3, 2] = self.mat[1, 2]

    def __getitem__(self, k):
        return self.mat[k[0], k[1]]

    def __setitem__(self, k, v):
        self.mat[k[0], k[1]] = v

    def rotate(self, theta, rotation_axis):
        '''Rotates the crystal axes an angle theta in radians around the axis defined by
        the rotation_axis expressed as the 3-tuple (nx, ny, nz).


        CHECK ME!!: This is based on NumBAT's orientation being x horizontal, y vertical, and z along the waveguide (out of the page for a right hand set).
        ''' 
        rmat = _rotate_Voigt_tensor(self.mat[1:, 1:], theta, rotation_axis)
        self.mat[1:, 1:] = rmat

    def __str__(self):
        s = f'\nVoigt tensor {self.sym}:\n'
        s += str(self.mat[1:, 1:])
        return s

    def dump(self):
        print(f'\nVoigt tensor {self.sym}')
        print(self.mat[1:, 1:])


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
    _data_loc = ''
    _materials = {}

    @classmethod
    def _set_file_locations(cls):
        this_dir = os.path.dirname(os.path.realpath(__file__))
        Material._data_loc = os.path.join(this_dir, "material_data", "")

    @classmethod
    def _load_materials_library(cls):
        Material._set_file_locations()
        for f in os.listdir(Material._data_loc):
            if f.endswith(".json"):
                try:
                    new_mat = Material(f)
                except BadMaterialFileError as err:
                    report_and_exit(str(err))

                # two mats have the same name  TODO: file_name is actually a mat_name
                if new_mat.material_name in Material._materials:
                    report_and_exit(
                        f"Material file {f} has the same name as an existing material {new_mat.material_name}.")

                Material._materials[new_mat.material_name] = new_mat

    @classmethod
    def _make_material(cls, s):
        if not Material._materials:  # first time through, we load all the materials
            Material._load_materials_library()

        try:
            mat = Material._materials[s]
        except KeyError:
            report_and_exit(
                f'Material {s} not found in material_data folder.\nEither the material file is missing or the name field in the material file has been incorrectly specified.')

        return mat

    # This should never be called directly by users, who should call make_material()

    def __init__(self, data_file):

        self.json_file = data_file


        # members to fill out amongst others
        self._crystal_axes=[]  # a,b,c crystal axes according to standard conventions
        
        try:
            self._load_data_file(self.json_file)

        except FileNotFoundError:
            report_and_exit(
                f'Material data {self.json_file} file not found.')

    def __str__(self):
        s = f'''
          Material: {self.chemical}
          File: {self.material_name}
          Source: {self.author}
          Date: {self.date}
          '''
        if len(self.comment):
            s += 'Comment: '+self.comment
        return s

    def full_str(self):
        s = str(self)
        s += str(self.c_tensor)
        s += str(self.eta_tensor)
        s += str(self.p_tensor)
        return s

    def elastic_properties(self):
        '''Returns a string containing key elastic properties of the material.'''
        try:
            s = f'Material:       {self.material_name}'
            s += f'\nDensity:        {self.rho:.3f} kg/m^3'

            if self.is_isotropic():
                s += f'\nc11:            {self.c_tensor.mat[1, 1]*1e-9:.3f} GPa'
                s += f'\nc12:            {self.c_tensor.mat[1, 2]*1e-9:.3f} GPa'
                s += f'\nc44:            {self.c_tensor.mat[4, 4]*1e-9:.3f} GPa'
                s += f"\nYoung's mod E:  {self.EYoung*1e-9:.3f} GPa"
                s += f'\nPoisson ratio:  {self.nuPoisson:.3f}'
                s += f'\nVelocity long.: {self.Vac_longitudinal():.3f} m/s'
                s += f'\nVelocity shear: {self.Vac_shear():.3f} m/s'
            else:
                s += '\nStiffness c:' + str(self.c_tensor)

        except Exception:
            s = 'Unknown/undefined elastic parameters in material '+self.material_name
        return s

    def Vac_longitudinal(self):
        '''For an isotropic material, returns the longitudinal (P-wave) elastic phase velocity.'''
        assert (self.is_isotropic())
        # lame lambda = c_12
        # lame mu = c_44
        #  v = math.sqrt(c_11/rho)
        #    =math.sqrt((c12 + 2c44)/rho)
        #    =math.sqrt((lambda + 2mu)/rho)
        if not self.rho or self.rho == 0:  # Catch vacuum cases
            return 0.
        else:
            return math.sqrt(self.c_tensor[1, 1]/self.rho)

    def Vac_shear(self):
        '''For an isotropic material, returns the shear (S-wave) elastic phase velocity.'''
        assert (self.is_isotropic())
        #  v = math.sqrt(c_44/rho)
        #    =math.sqrt((mu)/rho)
        if not self.rho or self.rho == 0:  # Catch vacuum cases
            return 0.
        else:
            return math.sqrt(self.c_tensor[4, 4]/self.rho)

    def has_elastic_properties(self):
        '''Returns true if the material has at least some elastic properties defined.'''
        return self.rho is not None

    def _load_data_file(self, data_file):
        """
        Load material data from json file.

        Args:
            data_file  (str): name of data file located in NumBAT/backend/material_data


        """

        with open(Material._data_loc + data_file, 'r') as fin:
            s_in = ''.join(fin.readlines())
            s_in = re.sub(r'//.*\n', '\n', s_in)

            try:
                self._params = json.loads(s_in)
            except Exception as err:
                traceback.print_exc()
                report_and_exit(
                    f'JSON parsing error: {err} for file {self.json_file}')

            # Name of this file, will be used as identifier and must be present
            self.material_name = self._params.get(
                'material_name', 'NOFILENAME')
            if self.material_name == 'NOFILENAME':
                raise BadMaterialFileError(
                    f"Material file {data_file} has no 'material_name' field.")

            self.format = self._params.get('format', 'NOFORMAT')
            if self.format == 'NOFORMAT':
                raise BadMaterialFileError(
                    f"Material file {data_file} has no 'format' field.")

            if self.format != 'NumBATMaterial-fmt-2.0':
                raise BadMaterialFileError(
                    f"Material file {data_file} must be in format 'NumBATMaterial-Fmt-2.0'.")

            self.chemical = self._params['chemical']  # Chemical composition
            self.author = self._params['author']  # Author of data
            # Year of data publication/measurement
            self.date = self._params['date']
            # Source institution
            self.institution = self._params['institution']
            # doi or, failing that, the http address
            self.doi = self._params['doi']

            # general comment for any purpose
            self.comment = self._params.get('comment', '')

            Re_n = self._params['Re_n']  # Real part of refractive index []
            # Imaginary part of refractive index []
            Im_n = self._params['Im_n']
            self.refindex_n = (Re_n + 1j*Im_n)  # Complex refractive index []
            self.rho = self._params['s']  # Density [kg/m3]
            self.EYoung = None
            self.nuPoisson = None


            #self.crystal_avec = np.array([1.0,0.0,0.0])
            #self.crystal_bvec = np.array([0.0,1.0,0.0])
            #self.crystal_cvec = np.array([0.0,0.0,1.0])


# self.c_11 = self._params['c_11']  # Stiffness tensor component [Pa]
#            self.c_12 = self._params['c_12']  # Stiffness tensor component [Pa]
#            self.c_44 = self._params['c_44']  # Stiffness tensor component [Pa]
#            self.p_11 = self._params['p_11']  # Photoelastic tensor component []
#            self.p_12 = self._params['p_12']  # Photoelastic tensor component []
#            self.p_44 = self._params['p_44']  # Photoelastic tensor component []
#            self.eta_11 = self._params['eta_11']  # Acoustic loss tensor component [Pa s]
#            self.eta_12 = self._params['eta_12']  # Acoustic loss tensor component [Pa s]
#            self.eta_44 = self._params['eta_44']  # Acoustic loss tensor component [Pa s]

            if not 'crystal_class' in self._params:
                raise BadMaterialFileError(
                    f"Material file {data_file} has no 'crystal_class' field.")
            try:
                self.crystal = CrystalGroup[self._params['crystal_class']]
            except ValueError as exc:
                print('Unknown crystal class in material data file')
                raise BadMaterialFileError(
                    f"Unknown crystal class in material data file {data_file}") from exc

            if self.crystal == CrystalGroup.Isotropic:
                self.load_crystal_isotropic()
                
            else:
                self.c_tensor = VoigtTensor4('c', self._params, self.json_file)
                self.eta_tensor = VoigtTensor4(
                    'eta', self._params, self.json_file)
                self.p_tensor = VoigtTensor4('p', self._params, self.json_file)
                #self.load_tensors()
                self.load_crystal_anisotropic()
            self._store_original_tensors()

    def _store_original_tensors(self):
        self._c_tensor_orig = self.c_tensor
        self._p_tensor_orig = self.p_tensor
        self._eta_tensor_orig  = self.eta_tensor

        
    def is_vacuum(self):
        '''Returns True if the material is the vacuum.'''
        return self.chemical == 'Vacuum'

    # (don't really need this as isotropic materials are the same)
    def load_crystal_cubic(self):

        # plain cartesian axes
        self.set_crystal_axes(np.array([1.0,0.,0.]),np.array([0.,1.0,0.]),np.array([0.0,0.,1.]))
        
        try:
            self.c_tensor.read(1, 1)
            self.c_tensor.read(1, 2)
            self.c_tensor[1, 3] = self.c_tensor[1, 2]
            self.c_tensor[2, 1] = self.c_tensor[1, 2]
            self.c_tensor[2, 2] = self.c_tensor[1, 1]
            self.c_tensor[2, 3] = self.c_tensor[1, 2]
            self.c_tensor[3, 1] = self.c_tensor[1, 2]
            self.c_tensor[3, 2] = self.c_tensor[1, 2]
            self.c_tensor[3, 3] = self.c_tensor[1, 1]
            self.c_tensor.read(4, 4)
            self.c_tensor[5, 5] = self.c_tensor[4, 4]
            self.c_tensor[6, 6] = self.c_tensor[4, 4]

            self.eta_tensor.read(1, 1)
            self.eta_tensor.read(1, 2)
            self.eta_tensor[1, 3] = self.eta_tensor[1, 2]
            self.eta_tensor[2, 1] = self.eta_tensor[1, 2]
            self.eta_tensor[2, 2] = self.eta_tensor[1, 1]
            self.eta_tensor[2, 3] = self.eta_tensor[1, 2]
            self.eta_tensor[3, 1] = self.eta_tensor[1, 2]
            self.eta_tensor[3, 2] = self.eta_tensor[1, 2]
            self.eta_tensor[3, 3] = self.eta_tensor[1, 1]
            self.eta_tensor.read(4, 4)
            self.eta_tensor[5, 5] = self.eta_tensor[4, 4]
            self.eta_tensor[6, 6] = self.eta_tensor[4, 4]

            self.p_tensor.read(1, 1)
            self.p_tensor.read(1, 2)

            self.p_tensor[1, 3] = self.p_tensor[1, 2]
            self.p_tensor[2, 1] = self.p_tensor[1, 2]
            self.p_tensor[2, 2] = self.p_tensor[1, 1]
            self.p_tensor[2, 3] = self.p_tensor[1, 2]
            self.p_tensor[3, 1] = self.p_tensor[1, 2]
            self.p_tensor[3, 2] = self.p_tensor[1, 2]
            self.p_tensor[3, 3] = self.p_tensor[1, 1]
            self.p_tensor.read(4, 4)

            # According to Powell, for Oh group, these are distinct elements, but no one seems to quote them
            if not self.p_tensor.read(5, 5, optional=True):
                self.p_tensor[5, 5] = self.p_tensor[4, 4]
            if not self.p_tensor.read(6, 6, optional=True):
                self.p_tensor[6, 6] = self.p_tensor[4, 4]

        except Exception:
            report_and_exit(
                f'Failed to load cubic crystal class in material data file {self.json_file}')

    def load_crystal_trigonal(self):
        # Good source for these rules is the supp info of doi:10.1364/JOSAB.482656 (Gustavo surface paper)

        self.set_crystal_axes(np.array([1.0,0.,0.]),np.array([0.,1.0,0.]),np.array([0.0,0.,1.]))
                                
        try:
            for lintens in [self.c_tensor, self.eta_tensor]:
                for (i, j) in [(1, 1), (1, 2), (1, 3), (1, 4), (3, 3), (4, 4) ]:
                    lintens.read(i, j)

            lintens[2, 1] = lintens[1, 2]
            lintens[2, 2] = lintens[1, 1]
            lintens[2, 3] = lintens[1, 3]
            lintens[2, 4] = -lintens[1, 4]

            lintens[3, 1] = lintens[1, 3]
            lintens[3, 2] = lintens[1, 3]

            lintens[4, 1] = lintens[1, 4]
            lintens[4, 2] = -lintens[1, 4]

            lintens[5, 5] = lintens[4, 4]
            lintens[5, 6] = lintens[1, 4]
            lintens[6, 5] = lintens[1, 4]
            lintens[6, 6] = (lintens[1, 1]-lintens[1,2])/2.0


            # TODO: confirm correct symmetry properties for p. 
            # PreviouslyuUsing trigonal = C3v from Powell, now the paper above
            self.p_tensor.read(1, 1)
            self.p_tensor.read(1, 2)
            self.p_tensor.read(1, 3)
            self.p_tensor.read(1, 4)
            self.p_tensor.read(3, 1)
            self.p_tensor.read(3, 3)
            self.p_tensor.read(4, 1)
            self.p_tensor.read(4, 4)

            self.p_tensor[2, 1] = self.p_tensor[1, 2]
            self.p_tensor[2, 2] = self.p_tensor[1, 1]
            self.p_tensor[2, 3] = self.p_tensor[1, 3]
            self.p_tensor[2, 4] = -self.p_tensor[1, 4]

            self.p_tensor[3, 2] = self.p_tensor[3, 1]

            self.p_tensor[4, 2] = -self.p_tensor[4, 1]

            self.p_tensor[5, 5] = self.p_tensor[4, 4]
            self.p_tensor[5, 6] = self.p_tensor[4, 1]
            self.p_tensor[6, 5] = self.p_tensor[1, 4]
            self.p_tensor[6, 6] = (self.p_tensor[1, 1] - self.p_tensor[1, 2])/2

        except Exception:
            report_and_exit(
                f'Failed to load trigonal crystal class in material data file {self.json_file}')

    def load_crystal_general(self):
        try:  # full anisotropic tensor components
            for i in range(1, 7):
                for j in range(1, 7):
                    self.c_tensor.read(i, j)
                    self.p_tensor.read(i, j)
                    self.eta_tensor.read(i, j)

        except KeyError:
            report_and_exit(
                'Failed to load anisotropic crystal class in material data file {self.json_file}')

    def set_refractive_index(self, nr, ni=0.0):
        self.refindex_n = nr + 1j*ni

    def is_isotropic(self): return not self._anisotropic


    def rotate_axis(self, theta, rotation_axis, save_rotated_tensors=False):
        self.rotate(theta, rotation_axis, save_rotated_tensors)
        
    def rotate(self, theta, rot_axis_spec, save_rotated_tensors=False):
        """ Rotate crystal axis by theta radians.

            Args:
                theta  (float): Angle to rotate by in radians.

                rotate_axis  (str): Axis around which to rotate.

            Keyword Args:
                save_rotated_tensors  (bool): Save rotated tensors to csv.

            Returns:
                ``Material`` object with rotated tensor values.
        """

        rotation_axis = parse_rotation_axis(rot_axis_spec)
        self.c_tensor.rotate(theta, rotation_axis)
        self.p_tensor.rotate(theta, rotation_axis)
        self.eta_tensor.rotate(theta, rotation_axis)

        caxes = self._crystal_axes.copy()
        self.set_crystal_axes(
            _rotate_3vector(caxes[0], theta, rotation_axis),
            _rotate_3vector(caxes[1], theta, rotation_axis),
            _rotate_3vector(caxes[2], theta, rotation_axis)
        )

        if save_rotated_tensors:
            np.savetxt('rotated_c_tensor.csv',
                       self.c_tensor.mat, delimiter=',')
            np.savetxt('rotated_p_tensor.csv',
                       self.p_tensor.mat, delimiter=',')
            np.savetxt('rotated_eta_tensor.csv',
                       self.eta_tensor.mat, delimiter=',')

    def reset_orientation(self):  #restore orientation to original axes in spec file.
        self.c_tensor = self._c_tensor_orig
        self.p_tensor = self._p_tensor_orig
        self.eta_tensor  = self._eta_tensor_orig

        self.set_crystal_axes(np.array([1.0,0.,0.]),np.array([0.,1.0,0.]),np.array([0.0,0.,1.]))

        #self.crystal_avec = np.array([1.0,0.0,0.0])
        #self.crystal_bvec = np.array([0.0,1.0,0.0])
        #self.crystal_cvec = np.array([0.0,0.0,1,0])

        
    def set_orientation(self, label):  # rotate original crystal to specific named-orientation, eg x-cut, y-cut. '111' etc.
        self.reset_orientation()

        try:
            ocode = self._params[f'orientation_{label.lower()}']
        except KeyError:
            report_and_exit(f'Orientation "{label}" is not defined for material {self.material_name}.')

        if ocode == 'ident':  # native orientation is the desired one
            return 
        
        try:
            ux, uy, uz, rot = map(float, ocode.split(','))
        except:
            report_and_exit(f"Can't parse crystal orientation code {ocode} for material {self.material_name}.")
        rot_axis = np.array((ux, uy, uz))
        theta = rot*np.pi/180

        self.rotate(theta, rot_axis)
        
    def plot_bulk_dispersion(self, pref):
        '''Draw slowness curve in the horizontal (x-z) plane for the crystal axes current orientation.

        Solving the Christoffel equation
        D C D^T u = -\rho v_p^2 u, for eigenvalue v_p and eigengector u.
        C is the Voigt form stiffness.
        D = [
            [κx  0   0   0  κz  κy  ]
            [0   κy  0   κz 0   κx  ]
            [0   0   κz  κy κx  0]
        ] where κ=(cos phi, 0, sin phi)
        '''
        #kx0, kx1 = -5.,5.
        #kz0, kz1 = -5.,5.
        npts = 500
        v_kphi = np.linspace(0.,np.pi*2,npts)
        v_vl = np.zeros(npts)
        v_vs1  = np.zeros(npts)
        v_vs2  = np.zeros(npts)

        with open(pref+'-bulkdisp.dat', 'w') as fout:

            fout.write('#phi    kapx    kapz    vl         vs1        vs2        vlx     vly      vlz     vs1x    vs1y     vs1z    vs2x    vs2y     vs2z\n')


            kapcomp = np.zeros(3)
            for ik,kphi in enumerate(v_kphi):
                kapx = np.cos(kphi)
                kapz = np.sin(kphi)
                kapy = 0.0
                vkap = np.array([kapx, kapy, kapz])
                mD = np.array([
                    [kapx, 0,    0,    0,    kapz, kapy],
                    [0,    kapy, 0,    kapz, 0,    kapx],
                    [0,    0,    kapz, kapy, kapx, 0]
                ])
                
                mLHS = np.matmul(np.matmul(mD, self.c_tensor.as_zerobase_matrix()), mD.T) / self.rho

                fout.write(f'{kphi:.4f}  {vkap[0]:+.4f}  {vkap[2]:+.4f}  ')

                #Solve and normalise
                evals, evecs = scipy.linalg.eig(mLHS)
                for i in range(3): evecs[:,i] /=np.linalg.norm(evecs[:,i]) # TODO: make a oneliner.
                vels = np.sqrt(np.real(evals))
                
                # identify the longitudinal mode (dot product with kappa = \pm 1)
                for i in range(3): kapcomp[i] = np.dot(vkap, evecs[:,i]) # TODO: make a oneliner

                if numbattools.almost_zero(1-np.abs(kapcomp[0])):
                    ivs = [0,1,2]  # indices of eigsols vl, vs1, vs2
                elif numbattools.almost_zero(1-np.abs(kapcomp[1])):
                    ivs = [1,0,2]
                else:
                    ivs = [2,0,1]
                    
                for iv in ivs: fout.write(f'{vels[iv]:.4f}  ')
                for iv in ivs: fout.write(f'{evecs[0,iv]:.4f}  {evecs[1,iv]:.4f}   {evecs[2,iv]:.4f}  ')
                fout.write(f'{kapcomp[0]:.4f}  {kapcomp[1]:.4f} {kapcomp[2]:.4f}')

                fout.write('\n')

                # Store vels for plotting in km/s
                v_vl[ik] = vels[ivs[0]] * 0.001
                v_vs1[ik] = vels[ivs[1]]* 0.001
                v_vs2[ik] = vels[ivs[2]]* 0.001
                
        #fig, ax = plt.subplots(subplot_kw={'projection':'polar'})
        fig, ax = plt.subplots()
        ax.plot(np.cos(v_kphi)/v_vl, np.sin(v_kphi)/v_vl, 'g', lw=1, markersize=1, label=r'$v_l$')
        ax.plot(np.cos(v_kphi)/v_vs1, np.sin(v_kphi)/v_vs1, lw=1, color='brown', label=r'$v_{s,i}$')
        ax.plot(np.cos(v_kphi)/v_vs2, np.sin(v_kphi)/v_vs2, lw=1, color='brown')
        #ax.set_rticks([])
        ax.set_xlabel(r'$1/v_x$ [s/km]')
        ax.set_ylabel(r'$1/v_z$ [s/km]')
        #ax.grid(False)
        ax.set_aspect(1.0)
        ax.axhline(0, c='gray')
        ax.axvline(0, c='gray')
        ax.legend(loc='upper right', frameon=False, fontsize=16)

        #decorator.add_extra_axes_commands()
        # add label
        ax.text(0.05, 0.95, self.material_name, fontsize=14, style='italic',
                transform=ax.transAxes)
        
        plt.savefig(pref+'-bulkdisp.png')
                    
    def make_crystal_plot(self, pref):
        '''Build crystal coordinates diagram using call to external asymptote application.'''
        
        fn = tempfile.NamedTemporaryFile(suffix='.asy', mode='w+t', delete=False)

        asy_cmds = get_asy_crystal_axes(self._crystal_axes)
        fn.write(asy_cmds)
        fn.close()

        # run .asy 
        subprocess.run(['asy', fn.name, '-o', f'{pref}crystal'])

    def set_crystal_axes(self, va, vb, vc):
        self._crystal_axes = [va, vb, vc]
        
    def load_crystal_isotropic(self):
        # ordinary Cartesian axes for the crystal axes
        self.set_crystal_axes(
            np.array([1.0,0.,0.]),
             np.array([0.,1.0,0.]),
             np.array([0.0,0.,1.]))
            
        
        self._anisotropic = False

        # Try to read isotropic from stiffness and then from Young's modulus and Poisson ratio
        if 'c_11' in self._params and 'c_12' in self._params and 'c_44' in self._params:
            self.c_tensor = VoigtTensor4(
                'c', self._params, self.json_file)
            self.c_tensor.load_isotropic()
            mu = self.c_tensor.mat[4, 4]
            lam = self.c_tensor.mat[1, 2]
            r = lam/mu
            self.nuPoisson = 0.5*r/(1+r)
            self.EYoung = 2*mu*(1+self.nuPoisson)

        elif 'EYoung' in self._params and 'nuPoisson' in self._params:
            self.EYoung = self._params['EYoung']
            self.nuPoisson = self._params['nuPoisson']
            c44 = 0.5*self.EYoung/(1+self.nuPoisson)
            c12 = self.EYoung*self.nuPoisson / \
                ((1+self.nuPoisson) * (1-2*self.nuPoisson))
            c11 = c12+2*c44
            self.c_tensor = VoigtTensor4('c')
            self.c_tensor.set_isotropic(c11, c12, c44)
        else:
            report_and_exit(
                'Broken isotropic material file:' + self.json_file)

        self.eta_tensor = VoigtTensor4(
            'eta', self._params, self.json_file)
        self.p_tensor = VoigtTensor4('p', self._params, self.json_file)

        self.p_tensor.load_isotropic()
        self.eta_tensor.load_isotropic()

    def load_crystal_anisotropic(self):  # not do this unless symmetry is off?

        self.c_tensor = VoigtTensor4('c', self._params, self.json_file)
        self.eta_tensor = VoigtTensor4('eta', self._params, self.json_file)
        self.p_tensor = VoigtTensor4('p', self._params, self.json_file)
                
        self._anisotropic = True

        #TODO: change to match/case
        if self.crystal == CrystalGroup.Trigonal:
            self.load_crystal_trigonal()
        elif self.crystal == CrystalGroup.Cubic:
            self.load_crystal_cubic()
        elif self.crystal == CrystalGroup.GeneralAnisotropic:
            self.load_crystal_general()

# Array that converts between 4th rank tensors in terms of x,y,z and Voigt notation
#               [[xx,xy,xz], [yx,yy,yz], [zx,zy,zz]]
to_Voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])


def rotate_tensor_elt(i, j, k, l, T_pqrs, mat_R):
    '''
    Calculates the element ijkl of the rotated tensor T' from the original
    rank-4 tensor T_PQ in 6x6 Voigt notation under the rotation specified by the 3x3 matrix R.
    '''

    Tp_ijkl = 0

    for q in range(3):
        for r in range(3):
            V1 = to_Voigt[q, r]
            for s in range(3):
                for t in range(3):
                    V2 = to_Voigt[s, t]
                    Tp_ijkl += mat_R[i, q] * mat_R[j, r] * \
                        mat_R[k, s] * mat_R[l, t] * T_pqrs[V1, V2]

    return Tp_ijkl

def parse_rotation_axis(rot_axis_spec):
    '''Convert one of several forms - string, numpy 3vec -  to a standard 3vec'''

    if isinstance(rot_axis_spec, str):
        ral = rot_axis_spec.lower()
        if ral in ('x', 'x-axis'):
            rot_axis = np.array([1.,0.,0.])
        elif ral in ('y', 'y-axis'):
            rot_axis = np.array([0.,1.,0.])
        elif ral in ('z', 'z-axis'):
            rot_axis = np.array([0.,0.,1.])
        else:
            report_and_exit(f"Can't convert {rot_axis_spec} to a 3-element unit vector.")
    else: # should be a numeric 3 vector

        try:
            if isinstance(rot_axis_spec, (tuple,list)):  # try to convert to numpy
                rot_axis = np.array(rot_axis_spec)
            elif isinstance(rot_axis_spec, np.ndarray):
                rot_axis = rot_axis_spec
            else:
                report_and_exit(emsg)
        except Exception:
            report_and_exit(emsg)
        if len(rot_axis) !=3:
            report_and_exit(f'Rotation axis {rot_axis} must have length 3.')

    return rot_axis

def _make_rotation_matrix(theta, rotation_axis): 
    """
    Return the SO(3) matrix corresponding to a rotation of theta radians the specified rotation_axis.
    """
    if isinstance(rotation_axis, str):
        if rotation_axis.lower() in ('x', 'x-axis'):
            mat_R = np.array([[1, 0, 0], [0, np.cos(
                theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
        elif rotation_axis.lower()  == ('y', 'y-axis'):
            mat_R = np.array([[np.cos(theta), 0, np.sin(theta)], [
                             0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
        elif rotation_axis.lower()  == ('z', 'z-axis'):
            mat_R = np.array([[np.cos(theta), -np.sin(theta), 0],
                             [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    else:
        emsg = f"Can't convert {rotation_axis} to a 3-element unit vector."
        try:
            if isinstance(rotation_axis, (tuple,list)):  # try to convert to numpy
                uvec = np.array(rotation_axis)
            elif isinstance(rotation_axis, np.ndarray):
                uvec = rotation_axis
            else:
                report_and_exit(emsg)
        except Exception:
            report_and_exit(emsg)

        # uvec is now a numpy array of some length
        nvec = np.linalg.norm(uvec)
        if len(uvec) != 3 or nvec == 0.0:
            report_and_exit(emsg)

        # normlise u
        uvec = uvec / nvec
        ct = math.cos(theta)
        st = math.sin(theta)
        omct = 1-ct
        ux, uy, uz = uvec[:]

        mat_R = np.array([
            [ct + ux**2*omct,  ux*uy*omct-uz*st, ux*uz*omct+uy*st],
            [uy*ux*omct+uz*st, ct + uy**2*omct,  uy*uz*omct-ux*st],
            [uz*ux*omct-uy*st, uz*uy*omct+ux*st, ct+uz**2*omct]
        ])

    return mat_R

def _rotate_3vector(vec, theta, rotation_axis):
    mat_R = _make_rotation_matrix(theta, rotation_axis)

    vrot = 0*vec
    for i in range(3):
        vrot[i] = mat_R[i,0]*vec[0] + mat_R[i,1]*vec[1] + mat_R[i,2]*vec[2]

    return vrot
        
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

    mat_R = _make_rotation_matrix(theta, rotation_axis)
    
    
    Tp_PQ = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            V1 = to_Voigt[i, j]
            for k in range(3):
                for l in range(3):
                    V2 = to_Voigt[k, l]
                    Tp_PQ[V1, V2] = rotate_tensor_elt(i, j, k, l, T_PQ, mat_R)

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


def make_material(s):
    return Material._make_material(s)


def get_asy_crystal_axes(crystal_axes):

    (va, vb, vc) = crystal_axes
    s_avec = '('+','.join(map(str, va))+')'
    s_bvec = '('+','.join(map(str, vb))+')'
    s_cvec = '('+','.join(map(str, vc))+')'
    
    s1='''
settings.outformat='png';
settings.render=8;
import three;
import graph3;

size(2cm,0);
defaultpen(fontsize(7pt));
defaultpen(.2);

real axlen=1.25;
int arrsize=3;
real blen=.5;

//currentprojection=orthographic(1,1,1);
currentprojection=oblique;


draw(O--2X, black, Arrow3(arrsize), L=Label("$\hat{x}$", position=EndPoint));
draw(O--2Y, black, Arrow3(arrsize), L=Label("$\hat{y}$", position=EndPoint));
draw(O--3Z, black, Arrow3(arrsize), L=Label("$\hat{z}$", position=EndPoint));

draw(O-- -2X, gray);
draw(O-- -2Y, gray);
draw(O-- -2Z, gray);


//label("$\hat{x}$", 3X*1.1);
//label("$\hat{y}$", 3Y*1.1);
//label("$\hat{z}$", 3Z*1.1);

draw(box((-1,-.5,-2)*blen,(1,.5,2)*blen),blue);
'''

    s2=f'''triple avec={s_avec};
triple bvec={s_bvec};
triple cvec={s_cvec};
'''

    s3 = '''triple corig=(0,.5,2)*blen;
draw(corig--avec+corig, red, Arrow3(arrsize), L=Label("$c_x$", position=EndPoint));
draw(corig--bvec+corig, red, Arrow3(arrsize), L=Label("$c_y$", position=EndPoint));
draw(corig--cvec+corig, red, Arrow3(arrsize), L=Label("$c_z$", position=EndPoint));

triple k0=(1,-1,-1);
triple k1=k0+(0,0,2);

draw(k0--k1,green, Arrow3(arrsize), L=Label("$k$"));
'''


    return s1 + s2 + s3
