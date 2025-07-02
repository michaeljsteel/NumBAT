# Copyright (C) 2017-2025  Michael Steel.

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

import math
import numpy as np
import numpy.linalg as npla
import copy

from nbtypes import unit_x, unit_y, unit_z
import reporting
import numbattools as nbtools

# Array that converts between 4th rank tensors in terms of x,y,z and Voigt notation
#               [[xx,xy,xz], [yx,yy,yz], [zx,zy,zz]]
to_Voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
from_Voigt = [
        None,
        (0,0),
        (1,1),
        (2,2),
        (1,2),
        (0,2),
        (0,1),
        ]



# These functions are indexed from 0.
# Note factors of 0.5 in the strain but not the stress

def strain_6col_to_3mat(Sv):
    return np.array([
        [   Sv[0], .5*Sv[5], .5*Sv[4]],
        [.5*Sv[5],    Sv[1], .5*Sv[3]],
        [.5*Sv[4], .5*Sv[3],    Sv[2]]])

def stress_6col_to_3mat(Tv):
    return np.array([
        [Tv[0], Tv[5], Tv[4]],
        [Tv[5], Tv[1], Tv[3]],
        [Tv[4], Tv[3], Tv[2]]])

def strain_3mat_to_6col(Sm):
    return np.array([Sm[0,0], Sm[1,1], Sm[2,2], 2*Sm[1,2], 2*Sm[0,2], 2*Sm[0,1] ])

def stress_3mat_to_6col(Tm):
    return np.array([Tm[0,0], Tm[1,1], Tm[2,2], Tm[1,2], Tm[0,2], Tm[0,1] ])


def kvec_to_symmetric_gradient(kvec):
    '''Calculates 6x3 symmetric gradient operator for plane wave of wavevector kvec=(kx,ky,kz).

    See Auld I, eq 1.53
    '''
    kx, ky, kz = kvec

    nabla_Ij = np.array([
        [kx,   0.0, 0.0],
        [0.0,  ky,  0.0],
        [0.0,  0.0, kz],
        [0,    kz,  ky],
        [kz,   0.0, kx],
        [ky,   kx, 0.0]
    ])
    return nabla_Ij

def rotate_2tensor_elt(i, j, T_pq, mat_R):
    '''
    Calculates the element ij of the rotated tensor Tp from the original
    rank-2 tensor T_pq in zero-indexed 3x3 notation under the rotation specified by the 3x3 matrix R.
    '''

    Tp_ij = 0

    for q in range(3):
        for r in range(3):
            Tp_ij += mat_R[i, q] * mat_R[j, r] * T_pq[q,r]

    return Tp_ij

def rotate_3tensor_elt(i, j, k, T_pqr, mat_R):
    '''
    Calculates the element ijk of the rotated tensor Tp from the original
    rank-3 tensor T_pqr in zero-indexed 3x3x3 notation under the rotation specified by the 3x3 matrix R.
    '''

    Tp_ijk = 0

    for q in range(3):
        for r in range(3):
            for s in range(3):
                    Tp_ijk += mat_R[i, q] * mat_R[j, r] * mat_R[k, s] * T_pqr[q,r,s]

    return Tp_ijk


def rotate_Voigt_4tensor_elt(i, j, k, l, T_PQ, mat_R):
    '''
    Calculates the element ijkl of the rotated tensor Tp from the original
    rank-4 tensor T_PQ in 6x6 Voigt notation under the rotation specified by the 3x3 matrix R.
    '''

    Tp_ijkl = 0

    for q in range(3):
        for r in range(3):
            V_I = to_Voigt[q, r]
            for s in range(3):
                for t in range(3):
                    V_J = to_Voigt[s, t]
                    Tp_ijkl += mat_R[i, q] * mat_R[j, r] * \
                        mat_R[k, s] * mat_R[l, t] * T_PQ[V_I, V_J]

    return Tp_ijkl


def parse_rotation_axis(rot_axis_spec):
    '''Convert one of several forms - string, numpy 3vec -  to a standard unit 3vec'''

    if isinstance(rot_axis_spec, str):
        ral = rot_axis_spec.lower()
        if ral in ('x', 'x-axis'):
            rot_axis = unit_x
        elif ral in ('y', 'y-axis'):
            rot_axis = unit_y
        elif ral in ('z', 'z-axis'):
            rot_axis = unit_z
        else:
            reporting.report_and_exit(
                f"Can't convert {rot_axis_spec} to a 3-element unit vector.")
    else:  # should be a numeric 3 vector
        emsg = f"Can't convert {rot_axis_spec} to a 3-element unit vector."
        try:
            if isinstance(rot_axis_spec, (tuple, list)):  # try to convert to numpy
                rot_axis = np.array(rot_axis_spec)
            elif isinstance(rot_axis_spec, np.ndarray):
                rot_axis = rot_axis_spec
            else:
                reporting.report_and_exit(emsg)
        except Exception:
            reporting.report_and_exit(emsg)
        if len(rot_axis) != 3:
            reporting.report_and_exit(
                f'Rotation axis {rot_axis} must have length 3.')

    nvec = np.linalg.norm(rot_axis)
    if nbtools.almost_zero(nvec):
        reporting.report_and_exit(f'Rotation axis {rot_axis} has zero length.')

    return rot_axis/nvec


def make_rotation_matrix(rot_axis_spec, theta):
    """
    Return the SO(3) matrix corresponding to a rotation of theta radians the specified rotation_axis.

    See https://en.wikipedia.org/wiki/Rotation_matrix
    """

    uvec = parse_rotation_axis(rot_axis_spec)

    ct = math.cos(theta)
    st = math.sin(theta)
    omct = 1-ct
    ux, uy, uz = uvec[:]

    mat_R = np.array([
        [ct + ux**2*omct,  ux*uy*omct-uz*st, ux*uz*omct+uy*st],
        [uy*ux*omct+uz*st, ct + uy**2*omct,  uy*uz*omct-ux*st],
        [uz*ux*omct-uy*st, uz*uy*omct+ux*st, ct+uz**2*omct]
    ])

    reporting.assertion(nbtools.almost_unity(
        np.linalg.det(mat_R)), 'Rotation matrix has unit determinant.')

    return mat_R


def rotate_3vector(vec3, mat_R):


    #vrot = 0*vec
    #for i in range(3):
    #    vrot[i] = mat_R[i, 0]*vec[0] + mat_R[i, 1]*vec[1] + mat_R[i, 2]*vec[2]
    #return vrot

    return np.matmul(mat_R, vec3)


def _rotate_2tensor(T_ij, mat_R):
    """
    Rotate a material tensor by theta radians around a specified rotation_axis.
    T_ij is a rank-2 tensor expressed in 3x3 zero-indexed standard notation.

    The complete operation in 3x3 notation is
    T'_ij  = sum_{pqr} R_ip R_jq R_kr

    Args:
        T_ij  (array): Tensor to be rotated.

        theta  (float): Angle to rotate by in radians.

        rotation_axis  (str): Axis around which to rotate.
    """

    Tp_ij = np.zeros((3, 3), dtype=T_ij.dtype)

    for i in range(3):
        for j in range(3):
                Tp_ij[i,j] = rotate_2tensor_elt(i, j, T_ij, mat_R)

    return Tp_ij

def _rotate_3tensor(T_ijk, mat_R):
    """
    Rotate a material tensor by theta radians around a specified rotation_axis.
    T_ijk is a rank-3 tensor expressed in 3x3x3 zero-indexed standard notation.

    The complete operation in 3x3x3 notation is
    T'_ijk  = sum_{pqr} R_ip R_jq R_kr T_pqr.

    Args:
        T_ijk  (array): Tensor to be rotated.

        theta  (float): Angle to rotate by in radians.

        rotation_axis  (str): Axis around which to rotate.
    """

    Tp_ijk = np.zeros((3, 3, 3), dtype=T_ijk.dtype)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                Tp_ijk[i,j,k] = rotate_3tensor_elt(i, j, k, T_ijk, mat_R)

    return Tp_ijk

def _rotate_Voigt_4tensor(T_PQ, mat_R):
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

    Tp_PQ = np.zeros((6, 6), dtype=T_PQ.dtype)
    for i in range(3):
        for j in range(3):
            V1 = to_Voigt[i, j]
            for k in range(3):
                for l in range(3):
                    V2 = to_Voigt[k, l]
                    Tp_PQ[V1, V2] = rotate_Voigt_4tensor_elt(i, j, k, l, T_PQ, mat_R)

    return Tp_PQ


def Voigt3_iJ_to_ijk(mat_iJ, fac2mul = False):
    "Takes Voigt indexed (0..2) x (1..6) iJ matrix to 3x3 (0..2,0..2,0..2) indexed matrix"
    T_ijk = np.zeros([3,3,3], dtype=mat_iJ.dtype)

    #fac2mul=False
    for i in range(3):
        for j in range(3):
            for k in range(3):
                J = to_Voigt[j,k] + 1  # +1 because toVoigt runs from 0..5
                if J>=4 and fac2mul:
                    T_ijk[i,j,k] = 2*mat_iJ[i,J]
                else:
                    T_ijk[i,j,k] = mat_iJ[i,J]
    return T_ijk

def Voigt3_ijk_to_iJ(mat_ijk, fac2mul = False):
    "Takes 3x3 (0..2,0..2,0..2) indexed matrix to Voigt indexed (0..2) x (1..6) iJ matrix"
    T_iJ = np.zeros([3,7], dtype=mat_ijk.dtype)

    #fac2mul=False
    for i in range(3):
        for J in range(1,7):
            (j,k) = from_Voigt[J]
            if J>=4 and fac2mul:
                T_iJ[i,J] = mat_ijk[i,j,k] / 2
            else:
                T_iJ[i,J] = mat_ijk[i,j,k]
    return T_iJ

class VoigtTensor3_iJ(object):
    '''A class for representing rank 3 tensors with iJ indexing
      (regular indexing in the first slot, Voigt in the last two slots)
    '''

    def __init__(self, json_symbol, physical_name='', physical_symbol='',
                 unit=None, transforms_with_factor_2 = False):

        self.mat = np.zeros([3, 7], dtype=float)  # unit indexing in the last two indices

        self._json_symbol = json_symbol  # eg 'piezo_d_IJ'
        self._physical_name=physical_name
        self._physical_symbol=physical_symbol
        self.unit=unit
        self._transforms_with_factor_2 = transforms_with_factor_2  # converts to _ijk form with stiffness facs

    def rotate(self, mat_R):
        with np.printoptions(precision=4):
            #print('rot3 a\n', self.mat)
            T_ijk = Voigt3_iJ_to_ijk(self.mat, self._transforms_with_factor_2)
            Tp_ijk = _rotate_3tensor(T_ijk, mat_R)
            self.mat = Voigt3_ijk_to_iJ(Tp_ijk, self._transforms_with_factor_2)

            #print('rot3 b\n', T_ijk)
            #print('rot3 c\n', Tp_ijk)
            #print('rot3 d\n', self.mat)
    #def to_Tijk(self):
    #    return Voigt3_iJ_to_ijk(self.mat, self._transforms_with_factor_2)

    #def from_Tijk(self, mat_ijk):
    #    self.mat = Voigt3_ijk_to_iJ(mat_ijk, self._transforms_with_factor_2)

    def as_transpose_iJ(self): # data expressed as dT_Ij[1:7][0:3]
        return self.mat.T

    def __str__(self):

        sh = f'\n  {self._physical_name} {self._physical_symbol}'
        if self.unit is not None:
            sh += f', unit: {self.unit[0]}.'
        else:
            sh += '.'

        sh += '\n    Voigt 3-tensor:\n'
        with np.printoptions(precision=4):
            if self.unit is not None:
                sd = str(self.mat[0:, 1:]/self.unit[1])
            else:
                sd = str(self.mat[0:, 1:])

        return sh + nbtools.indent_string(sd, indent=4)


    def set_from_00matrix(self, mat00):
        self.mat[:, 1:7] = mat00

    def set_from_json(self, json_data):

        for ci,xyz in enumerate(('x', 'y', 'z')):

            for J in range(1,7):
                elt = f'{self._json_symbol}_{xyz}{J}'
                if elt in json_data:
                    self.mat[ci][J] = json_data[elt]

    def set_elt_from_json(self, s_iJ, json_data):

        elt = f'{self._json_symbol}_{s_iJ}'
        e_i = {'x':0, 'y':1, 'z':2}[s_iJ[0]]
        e_J = int(s_iJ[1])
        self.mat[e_i, e_J] =  json_data[elt]

    def elt_iJ(self, s_iJ):
        e_i = {'x':0, 'y':1, 'z':2}[s_iJ[0]]
        e_J = int(s_iJ[1])
        return self.mat[e_i, e_J]

    def set_elt_iJ(self, s_iJ, val):
        e_i = {'x':0, 'y':1, 'z':2}[s_iJ[0]]
        e_J = int(s_iJ[1])
        self.mat[e_i, e_J] = val



    def value(self):
        '''Returns copy of Voigt matrix indexed as m[0..2, 0..5].'''
        return self.mat[:, 1:].copy()

    def refvalue(self):
        '''Returns reference to internal Voigt matrix indexed as m[0..2, 0..5].'''
        return self.mat[:, 1:]

    def copy(self):
        return copy.deepcopy(self)

class VoigtTensor4(object):
    '''A class for representing rank 4 tensors in the compact Voigt representation.

    Internally uses 1-based indexing like the notation. Externally it passes zero-based values.'''

    def __init__(self, src_dict, json_symbol, physical_name='', physical_symbol='', unit=None):
        """unit is a symbol, scale value tuple:  eg ('GPa', 1.e9)"""

        self.mat = np.zeros([7, 7], dtype=float)  # unit indexing

        self.physical_name = physical_name
        self.json_symbol = json_symbol  # eg 'c', 'p', 'eta'

        self._json_data = src_dict
        self.physical_symbol=physical_symbol
        self.unit=unit

    def __str__(self):

        sh = f'\n  {self.physical_name} {self.physical_symbol}'
        if self.unit is not None:
            sh += f', unit: {self.unit[0]}.'
        else:
            sh += '., unit: dimensionless.'

        sh += '\n    Voigt 4-tensor:\n'
        with np.printoptions(precision=4):
            if self.unit is not None:
                sd = str(self.mat[1:, 1:]/self.unit[1])
            else:
                sd = str(self.mat[1:, 1:])

        return sh + nbtools.indent_string(sd, indent=4)

    # Allow direct indexing of Voigt tensor in [(i,j)] form

    def __getitem__(self, k):
        return self.mat[k[0], k[1]]

    def __setitem__(self, k, v):
        self.mat[k[0], k[1]] = v

    def elt_s_IJ(self, s_IJ):
        eI=int(s_IJ[0])
        eJ=int(s_IJ[1])
        assert(eI>=1 and eI<=6 and eJ>=1 and eJ<=6)
        return self.mat[eI,eJ]

    def set_elt_s_IJ(self, s_IJ, val):
        eI=int(s_IJ[0])
        eJ=int(s_IJ[1])
        assert(eI>=1 and eI<=6 and eJ>=1 and eJ<=6)
        self.mat[eI,eJ] = val

    def dump_rawdata(self):
        print(f'\nVoigt tensor {self.physical_name}, tensor {self.json_symbol}')
        print(self.mat[1:, 1:])

    def value(self):
        '''Returns copy of Voigt matrix indexed as m[0..5, 0..5].'''
        return self.mat[1:, 1:].copy()

    def refvalue(self):
        '''Returns reference to internal Voigt matrix indexed as m[0..5, 0..5].'''
        return self.mat[1:, 1:]


    def set_from_00matrix(self, mat00):
        if self.mat.dtype != mat00.dtype:
            self.mat = np.zeros([7, 7], dtype=mat00.dtype)
        self.mat[1:7, 1:7] = mat00

    def copy(self):
        return copy.deepcopy(self)

    def read_from_json(self, m, n, optional=False):
        '''Looks for data in the _json_data dict in form c_12, eta_23, etc.'''

        elt = f'{self.json_symbol}_{m}{n}'

        if elt not in self._json_data:
            if not optional:
                reporting.report_and_exit(
                    f'Failed to read required tensor element {elt} for material {self.physical_name}')
            else:
                return False

        self.mat[m, n] = self._json_data[elt]
        return True

    def load_isotropic_from_json(self):
        self.read_from_json(1, 1)
        self.read_from_json(1, 2)
        self.read_from_json(4, 4)
        self.make_isotropic_tensor(self.mat[1, 1], self.mat[1, 2], self.mat[4, 4])

    def make_isotropic_tensor(self, m11, m12, m44):
        '''Build Voigt matrix from 3 parameters for isotropic geometry.
        (Actually, only two are independent.)'''

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

    def check_symmetries(self):
        # Check matrix is symmetric and positive definite

        rtol = 1e-12
        tol = rtol * np.abs(self.mat).max()
        tmat = self.mat - self.mat.T
        mat_is_sym = nbtools.almost_zero(np.linalg.norm(tmat), tol)
        reporting.assertion(
            mat_is_sym, f'Material matrix {self.physical_name}-{self.physical_symbol} is symmetric.\n' + str(self.mat))

    def rotate(self, mat_R):
        '''Rotates the crystal according to the SO(3) matrix mat_R.
        '''

        rot_tensor = _rotate_Voigt_4tensor(self.mat[1:, 1:], mat_R)
        self.mat[1:, 1:] = rot_tensor

    def inverse(self, physical_name='', physical_symbol='', unit=None):
        """Return a new Voigt-4 tensor that is the inverse of this one."""

        m_data = self.value()
        m_dinv = npla.inv(m_data)

        vt_inv = self.copy()
        vt_inv.set_from_00matrix(m_dinv)

        if physical_name:
            vt_inv.physical_name=physical_name
        if physical_symbol:
            vt_inv.physical_symbol=physical_symbol
        if unit:
            vt_inv.unit=unit

        return vt_inv
