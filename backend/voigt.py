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
import numpy.random as nprand
import copy

from nbtypes import unit_x, unit_y, unit_z
import reporting
import numbattools as nbtools

# Array that converts between 4th rank tensors in terms of x,y,z and Voigt notation
#               [[xx,xy,xz], [yx,yy,yz], [zx,zy,zz]]

# Makes a Voigt3 index in [1..6] from (i,j) in [0..2]x[0..2]
to_Voigt3_index = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]]) + 1
"""Array mapping (i, j) indices in [0..2]x[0..2] to Voigt notation indices [1..6]."""

#used by Voigt4, J is 0..5 !
to_Voigt_zerobase = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
"""Array mapping (i, j) indices in [0..2]x[0..2] to zero-based Voigt indices [0..5]."""

# Sends a Voigt index in [1..6] to (i,j) in [0..2]x[0..2]
from_Voigt = [
        None,
        (0,0),
        (1,1),
        (2,2),
        (1,2),
        (0,2),
        (0,1),
        ]
"""List mapping Voigt indices [1..6] to (i, j) pairs in [0..2]x[0..2]."""


# These functions are indexed from 0.
# Note factors of 0.5 in the strain but not the stress

def strain_6col_to_3mat(Sv):
    """
    Convert a 6-element Voigt strain vector to a 3x3 symmetric strain matrix.

    Args:
        Sv (array-like): 6-element strain vector in Voigt notation.

    Returns:
        np.ndarray: 3x3 symmetric strain matrix.
    """
    return np.array([
        [   Sv[0], .5*Sv[5], .5*Sv[4]],
        [.5*Sv[5],    Sv[1], .5*Sv[3]],
        [.5*Sv[4], .5*Sv[3],    Sv[2]]])

def stress_6col_to_3mat(Tv):
    """
    Convert a 6-element Voigt stress vector to a 3x3 symmetric stress matrix.

    Args:
        Tv (array-like): 6-element stress vector in Voigt notation.

    Returns:
        np.ndarray: 3x3 symmetric stress matrix.
    """
    return np.array([
        [Tv[0], Tv[5], Tv[4]],
        [Tv[5], Tv[1], Tv[3]],
        [Tv[4], Tv[3], Tv[2]]])

def strain_3mat_to_6col(Sm):
    """
    Convert a 3x3 symmetric strain matrix to a 6-element Voigt strain vector.

    Args:
        Sm (np.ndarray): 3x3 symmetric strain matrix.

    Returns:
        np.ndarray: 6-element strain vector in Voigt notation.
    """
    return np.array([Sm[0,0], Sm[1,1], Sm[2,2], 2*Sm[1,2], 2*Sm[0,2], 2*Sm[0,1] ])

def stress_3mat_to_6col(Tm):
    """
    Convert a 3x3 symmetric stress matrix to a 6-element Voigt stress vector.

    Args:
        Tm (np.ndarray): 3x3 symmetric stress matrix.

    Returns:
        np.ndarray: 6-element stress vector in Voigt notation.
    """
    return np.array([Tm[0,0], Tm[1,1], Tm[2,2], Tm[1,2], Tm[0,2], Tm[0,1] ])


def kvec_to_symmetric_gradient(kvec):
    """
    Calculate the 6x3 symmetric gradient operator for a plane wave with wavevector kvec.

    Args:
        kvec (array-like): Wavevector (kx, ky, kz).

    Returns:
        np.ndarray: 6x3 symmetric gradient operator matrix.

    Reference:
        Auld I, eq 1.53
    """
    kx, ky, kz = kvec

    nabla_IJ = np.array([
        [kx,   0.0, 0.0],
        [0.0,  ky,  0.0],
        [0.0,  0.0, kz],
        [0,    kz,  ky],
        [kz,   0.0, kx],
        [ky,   kx, 0.0]
    ])
    return nabla_IJ

# def rotate_2tensor_elt(i, j, T_pq, mat_R):
#     """
#     Calculate the element (i, j) of a rotated rank-2 tensor under a given rotation.

#     Args:
#         i (int): Row index of the rotated tensor.
#         j (int): Column index of the rotated tensor.
#         T_pq (np.ndarray): Original 3x3 tensor.
#         mat_R (np.ndarray): 3x3 rotation matrix.

#     Returns:
#         float: Rotated tensor element Tp_ij.
#     """

#     Tp_ij = 0

#     for q in range(3):
#         for r in range(3):
#             Tp_ij += mat_R[i, q] * mat_R[j, r] * T_pq[q,r]

#     return Tp_ij

# def rotate_3tensor_elt(i, j, k, T_pqr, mat_R):
#     """
#     Calculate the element (i, j, k) of a rotated rank-3 tensor under a given rotation.

#     Args:
#         i (int): First index of the rotated tensor.
#         j (int): Second index of the rotated tensor.
#         k (int): Third index of the rotated tensor.
#         T_pqr (np.ndarray): Original 3x3x3 tensor.
#         mat_R (np.ndarray): 3x3 rotation matrix.

#     Returns:
#         float: Rotated tensor element Tp_ijk.
#     """

#     Tp_ijk = 0

#     for q in range(3):
#         for r in range(3):
#             for s in range(3):
#                     Tp_ijk += mat_R[i, q] * mat_R[j, r] * mat_R[k, s] * T_pqr[q,r,s]

#     return Tp_ijk


# def rotate_Voigt_4tensor_elt(i, j, k, l, T_PQ, mat_R):
#     """
#     Calculate the element (i, j, k, L) of a rotated rank-4 tensor in 6x6 Voigt notation under a given rotation.

#     Args:
#         i, j, k, L (int): Indices of the rotated tensor.
#         T_PQ (np.ndarray): Original 6x6 Voigt tensor.
#         mat_R (np.ndarray): 3x3 rotation matrix.

#     Returns:
#         float: Rotated tensor element Tp_ijkl.
#     """

#     Tp_ijkl = 0

#     for q in range(3):
#         for r in range(3):
#             V_I = to_Voigt_zerobase[q, r]
#             for s in range(3):
#                 for t in range(3):
#                     V_J = to_Voigt_zerobase[s, t]
#                     Tp_ijkl += mat_R[i, q] * mat_R[j, r] * \
#                         mat_R[k, s] * mat_R[l, t] * T_PQ[V_I, V_J]

#     return Tp_ijkl


def parse_rotation_axis(rot_axis_spec):
    """
    Convert a rotation axis specification to a standard unit 3-vector.

    Args:
        rot_axis_spec (str, tuple, list, or np.ndarray): Axis specification (e.g., 'x', 'y', 'z', or a 3-vector).

    Returns:
        np.ndarray: Normalized 3-element unit vector.
    """

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
    Return the SO(3) matrix corresponding to a rotation of theta radians about the specified axis.

    Args:
        rot_axis_spec (str, tuple, list, or np.ndarray): Axis specification.
        theta (float): Angle in radians.

    Returns:
        np.ndarray: 3x3 rotation matrix.
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

def generate_voigt_rotation_matrix(matR):
    """
    Generates the 6x6 transformation matrix M for rotating a rank 4 tensor
    expressed in 6x6 Voigt notation, given a 3x3 rotation matrix R.

    The transformation is C_rotated = M @ C_unrotated @ M.T

    Args:
        matR (numpy.ndarray): A 3x3 NumPy array representing the rotation matrix R.

    Returns:
        numpy.ndarray: A 6x6 NumPy array representing the transformation matrix M.
                       Returns None if matR is not a 3x3 matrix.
    """
    if matR.shape != (3, 3):
        print("Error: Input matrix matR must be a 3x3 NumPy array.")
        return None

    M = np.zeros((6, 6))

    # Extract elements of R for clarity
    R11, R12, R13 = matR[0, 0], matR[0, 1], matR[0, 2]
    R21, R22, R23 = matR[1, 0], matR[1, 1], matR[1, 2]
    R31, R32, R33 = matR[2, 0], matR[2, 1], matR[2, 2]

    # Populate the M matrix according to the derived formula
    # First 3x3 block (normal-normal components)
    M[0, 0] = R11**2; M[0, 1] = R12**2; M[0, 2] = R13**2
    M[1, 0] = R21**2; M[1, 1] = R22**2; M[1, 2] = R23**2
    M[2, 0] = R31**2; M[2, 1] = R32**2; M[2, 2] = R33**2

    # First 3x3 block (normal-shear components, factor of 2)
    M[0, 3] = 2 * R12 * R13
    M[0, 4] = 2 * R11 * R13
    M[0, 5] = 2 * R11 * R12

    M[1, 3] = 2 * R22 * R23
    M[1, 4] = 2 * R21 * R23
    M[1, 5] = 2 * R21 * R22

    M[2, 3] = 2 * R32 * R33
    M[2, 4] = 2 * R31 * R33
    M[2, 5] = 2 * R31 * R32

    # Bottom 3x3 block (shear-normal components)
    M[3, 0] = R21 * R31
    M[3, 1] = R22 * R32
    M[3, 2] = R23 * R33

    M[4, 0] = R11 * R31
    M[4, 1] = R12 * R32
    M[4, 2] = R13 * R33

    M[5, 0] = R11 * R21
    M[5, 1] = R12 * R22
    M[5, 2] = R13 * R23

    # Bottom 3x3 block (shear-shear components)
    M[3, 3] = R22 * R33 + R23 * R32
    M[3, 4] = R21 * R33 + R23 * R31
    M[3, 5] = R21 * R32 + R22 * R31

    M[4, 3] = R12 * R33 + R13 * R32
    M[4, 4] = R11 * R33 + R13 * R31
    M[4, 5] = R11 * R32 + R12 * R31

    M[5, 3] = R12 * R23 + R13 * R22
    M[5, 4] = R11 * R23 + R13 * R21
    M[5, 5] = R11 * R22 + R12 * R21

    return M





def rotate_3vector(vec3, mat_R):
    """
    Rotate a 3-vector using a given rotation matrix.

    Args:
        vec3 (np.ndarray): 3-element vector.
        mat_R (np.ndarray): 3x3 rotation matrix.

    Returns:
        np.ndarray: Rotated 3-element vector.
    """

    return np.matmul(mat_R, vec3)


def _rotate_2tensor(T_ij, mat_R):
    """
    Rotate a rank-2 tensor using a rotation matrix.

    Args:
        T_ij (np.ndarray): 3x3 tensor to be rotated.
        mat_R (np.ndarray): 3x3 rotation matrix.

    Returns:
        np.ndarray: Rotated 3x3 tensor.
    """

    Tp_ij = np.matmul(mat_R, np.matmul(T_ij, mat_R.T))

    return Tp_ij


# def _old_rotate_2tensor(T_ij, mat_R):
#     """
#     Rotate a material tensor by theta radians around a specified rotation_axis.
#     T_ij is a rank-2 tensor expressed in 3x3 zero-indexed standard notation.

#     The complete operation in 3x3 notation is
#     T'_ij  = sum_{pqr} R_ip R_jq T_pq

#     Args:
#         T_ij  (array): Tensor to be rotated.

#         theta  (float): Angle to rotate by in radians.

#         rotation_axis  (str): Axis around which to rotate.
#     """

#     Tp_ij = np.zeros((3, 3), dtype=T_ij.dtype)

#     for i in range(3):
#         for j in range(3):
#                 Tp_ij[i,j] = rotate_2tensor_elt(i, j, T_ij, mat_R)

#     return Tp_ij

# def _rotate_3tensor(T_ijk, mat_R):
#     """
#     Rotate a material tensor by theta radians around a specified rotation_axis.
#     T_ijk is a rank-3 tensor expressed in 3x3x3 zero-indexed standard notation.

#     The complete operation in 3x3x3 notation is
#     T'_ijk  = sum_{pqr} R_ip R_jq R_kr T_pqr.

#     Args:
#         T_ijk (np.ndarray): 3x3x3 tensor to be rotated.
#         mat_R (np.ndarray): 3x3 rotation matrix.

#     Returns:
#         np.ndarray: Rotated 3x3x3 tensor.
#     """

#     Tp_ijk = np.zeros((3, 3, 3), dtype=T_ijk.dtype)

#     for i in range(3):
#         for j in range(3):
#             for k in range(3):
#                 Tp_ijk[i,j,k] = rotate_3tensor_elt(i, j, k, T_ijk, mat_R)

#     return Tp_ijk

# def _rotate_Voigt_4tensor(T_PQ, mat_R):
#     """
#     Rotate a rank-4 acoustic material tensor in 6x6 Voigt notation using a rotation matrix.

#     The complete operation in 3x3x3 notation is
#     T'_ijkl  = sum_{pqr} R_ip R_jq R_kr R_ls T_pqrs.

#     Args:
#         T_PQ (np.ndarray): 6x6 Voigt tensor to be rotated.
#         mat_R (np.ndarray): 3x3 rotation matrix.

#     Returns:
#         np.ndarray: Rotated 6x6 Voigt tensor.
#     """

#     Tp_PQ = np.zeros((6, 6), dtype=T_PQ.dtype)
#     for i in range(3):
#         for j in range(3):
#             V1 = to_Voigt_zerobase[i, j]
#             for k in range(3):
#                 for l in range(3):
#                     V2 = to_Voigt_zerobase[k, l]
#                     Tp_PQ[V1, V2] = rotate_Voigt_4tensor_elt(i, j, k, l, T_PQ, mat_R)

#     return Tp_PQ


def Voigt3_iJ_to_ijk(mat_iJ, fac2mul = False):
    """
    Convert a 3x6 Voigt-indexed tensor to a 3x3x3 tensor, optionally multiplying off-diagonal elements by 2.

    Args:
        mat_iJ (np.ndarray): 3x6 Voigt-indexed tensor.
        fac2mul (bool): Whether to multiply off-diagonal elements by 2.

    Returns:
        np.ndarray: 3x3x3 tensor.
    """

    T_ijk = np.zeros([3,3,3], dtype=mat_iJ.dtype)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                J = to_Voigt3_index[j,k]
                if J>=4 and fac2mul:
                    T_ijk[i,j,k] = 2*mat_iJ[i,J-1]
                else:
                    T_ijk[i,j,k] = mat_iJ[i,J-1]
    return T_ijk

def Voigt3_ijk_to_iJ(mat_ijk, fac2mul = False):
    """
    Convert a 3x3x3 tensor to a 3x6 Voigt-indexed tensor, optionally dividing off-diagonal elements by 2.

    Args:
        mat_ijk (np.ndarray): 3x3x3 tensor.
        fac2mul (bool): Whether to divide off-diagonal elements by 2.

    Returns:
        np.ndarray: 3x6 Voigt-indexed tensor.
    """

    T_iJ = np.zeros([3,6], dtype=mat_ijk.dtype)

    for i in range(3):
        for J in range(1,7):
            (j,k) = from_Voigt[J]
            if J>=4 and fac2mul:
                T_iJ[i,J-1] = mat_ijk[i,j,k] / 2
            else:
                T_iJ[i,J-1] = mat_ijk[i,j,k]
    return T_iJ


class PlainTensor2_ij(object):
    '''A class for representing rank 2 tensors with ij indexing in [0..2, 0..2]
      (regular indexing in the first slot, Voigt in the last two slots)
    '''

    def __init__(self, param_nm, physical_name='', physical_symbol='',
                 unit=None, is_complex=False, is_symmetric=False):
        """
        Initialize a PlainTensor2_ij object.

        Args:
            param_nm (str): Symbol head in parameter dictionary (e.g., 'piezo_epsS').
            physical_name (str): Human-readable name for the tensor.
            physical_symbol (str): Symbol for display.
            unit (tuple or None): Unit information (name, scale).
            is_complex (bool): Whether the tensor is complex-valued.
            is_symmetric (bool): Whether the tensor is symmetric.
        """
        dtype = np.complex128 if is_complex else np.float64

        self._is_complex = is_complex
        self._is_symmetric = is_symmetric

        self.mat = np.zeros([3, 3], dtype=dtype)

        self._param_nm = param_nm  # eg 'piezo_epsS'
        self._physical_name=physical_name
        self._physical_symbol=physical_symbol
        self.unit=unit

    def as_str(self, chop=False, rtol=1e-10, atol=1e-12):
        """
        Return a string representation of the tensor, formatted for display.

        Args:
            chop (bool): Whether to chop small values to zero.
            rtol (float): Relative tolerance for chopping.
            atol (float): Absolute tolerance for chopping.

        Returns:
            str: Formatted string representation of the tensor.
        """

        sh = f'\n {self._physical_name} {self._physical_symbol}'
        if self.unit is not None:
            sh += f', unit: {self.unit[0]}.'
        else:
            sh += '., unit: dimensionless.'

        sh += '\n   Plain 2-tensor:\n'
        with np.printoptions(precision=4):
            if self.unit is not None:
                sd = str(nbtools.chopmat(self.mat[1:, 1:]/self.unit[1], chop=chop, rtol=rtol, atol=atol))
            else:
                sd = str(nbtools.chopmat(self.mat[1:, 1:], chop=chop, rtol=rtol, atol=atol))

        return sh + nbtools.indent_string(sd, indent=4)

    def __str__(self):
        """
        Return a formatted string representation of the tensor, including units and values.
        """

        return self.as_str(chop=False)


    def set_from_params(self, d_params):
        """
        Set tensor elements from a parameter dictionary, handling symmetry and complex values.

        Args:
            d_params (dict): Dictionary of parameter values.
        """

        # must do xy before yx etc
        for i_xyz, is_xyz in enumerate(('x', 'y', 'z')):
            for j_xyz, js_xyz in enumerate(('x', 'y', 'z')):

                elt = f'{self._param_nm}_{is_xyz}{js_xyz}'
                if self._is_symmetric and i_xyz > j_xyz:
                    elt = f'{self._param_nm}_{js_xyz}{is_xyz}'

                if self._is_complex:
                    val = d_params.get(elt+'_r', 0) + 1j*d_params.get(elt+'_i',0.0)
                else:
                    val = d_params.get(elt, 0)

                self.mat[i_xyz, j_xyz] = val

    def set_from_00matrix(self, mat00):
        """
        Set the tensor values from a 3x3 matrix.

        Args:
            mat00 (np.ndarray): 3x3 matrix.
        """

        self.mat = mat00.copy()

    def copy(self):
        """
        Return a deep copy of the tensor object.
        """

        return copy.deepcopy(self)

    def value(self):
        """
        Return a copy of the tensor's matrix data.
        """

        return self.mat.copy()

    def fill_random(self):
        """
        Fill the tensor with random values (uniform in [0,1)).
        """

        self.mat = nprand.rand(3,3)

    def rotate(self, mat_R):
        """
        Rotate the tensor using the SO(3) matrix mat_R.

        Args:
            mat_R (np.ndarray): 3x3 rotation matrix.
        """

        self.mat = _rotate_2tensor(self.mat, mat_R)


class VoigtTensor3_iJ(object):
    '''A class for representing rank 3 tensors with iJ indexing
      (regular indexing in the first slot, Voigt in the last two slots)
    '''

    def __init__(self, param_nm, physical_name='', physical_symbol='',
                 unit=None, transforms_with_factor_2 = False):
        """
        Initialize a VoigtTensor3_iJ object.

        Args:
            param_nm (str): Symbol head in parameter dictionary (e.g., 'piezo_d_IJ').
            physical_name (str): Human-readable name for the tensor.
            physical_symbol (str): Symbol for display.
            unit (tuple or None): Unit information (name, scale).
            transforms_with_factor_2 (bool): Whether to apply factor 2 in Voigt conversion.
        """

        self.mat = np.zeros([3, 6], dtype=float)  # unit indexing in the last two indices

        self._param_nm = param_nm  # eg 'piezo_d_IJ'
        self._physical_name=physical_name
        self._physical_symbol=physical_symbol
        self.unit=unit
        self._transforms_with_factor_2 = transforms_with_factor_2  # converts to _ijk form with stiffness facs

    def _old_rotate(self, mat_R):
        """
        Rotate the tensor using the SO(3) matrix mat_R.

        Args:
            mat_R (np.ndarray): 3x3 rotation matrix.
        """

        T_ijk = Voigt3_iJ_to_ijk(self.mat, self._transforms_with_factor_2)
        Tp_ijk = _rotate_3tensor(T_ijk, mat_R)
        self.mat = Voigt3_ijk_to_iJ(Tp_ijk, self._transforms_with_factor_2)

    def rotate(self, mat_R):
        """
        Rotate the tensor using the SO(3) matrix mat_R, returning a new tensor.

        Args:
            mat_R (np.ndarray): 3x3 rotation matrix.

        """

        mat_M = generate_voigt_rotation_matrix(mat_R)

        self.mat = mat_R @ self.mat @ mat_M.T




    def fill_random(self):
        """
        Fill the tensor with random values (uniform in [0,1)).
        """

        self.mat = nprand.rand(3,6)

    def as_transpose_iJ(self): # data expressed as dT_Ij[0:6][0:3]
        """
        Return the transpose of the tensor (shape [6,3]).
        """

        return self.mat.T


    def as_str(self, chop=False, rtol=1e-10, atol=1e-12):
        """
        Return a string representation of the tensor, formatted for display.

        Args:
            chop (bool): Whether to chop small values to zero.
            rtol (float): Relative tolerance for chopping.
            atol (float): Absolute tolerance for chopping.

        Returns:
            str: Formatted string representation of the tensor.
        """

        sh = f'\n {self._physical_name} {self._physical_symbol}'
        if self.unit is not None:
            sh += f', unit: {self.unit[0]}.'
        else:
            sh += '., unit: dimensionless.'

        sh += '\n   Voigt 3-tensor:\n'
        with np.printoptions(precision=4):
            if self.unit is not None:
                sd = str(nbtools.chopmat(self.mat[1:, 1:]/self.unit[1], chop=chop, rtol=rtol, atol=atol))
            else:
                sd = str(nbtools.chopmat(self.mat[1:, 1:], chop=chop, rtol=rtol, atol=atol))

        return sh + nbtools.indent_string(sd, indent=4)

    def __str__(self):
        """
        Return a formatted string representation of the tensor, including units and values.
        """

        return self.as_str(chop=False)


    def set_from_00matrix(self, mat00):
        """
        Set the tensor values from a 3x6 matrix.

        Args:
            mat00 (np.ndarray): 3x6 matrix.
        """

        self.mat = mat00.copy()

    def set_from_params(self, d_params):
        """
        Set tensor elements from a parameter dictionary.

        Args:
            d_params (dict): Dictionary of parameter values.
        """

        for ci,xyz in enumerate(('x', 'y', 'z')):

            for J in range(1,7):
                elt = f'{self._param_nm}_{xyz}{J}'
                if elt in d_params:
                    self.mat[ci][J-1] = d_params[elt]

    def _siJ_to_rc(self, s_iJ):
        """
        Convert a string index (e.g., 'x1') to row and column indices.

        Args:
            s_iJ (str): String index (e.g., 'x1', 'y3').

        Returns:
            tuple: (row, column) indices.
        """

        e_i = {'x':0, 'y':1, 'z':2}[s_iJ[0]]
        e_J = int(s_iJ[1]) -1  # convert to zero-based index
        return e_i, e_J

    def set_elt_from_param_dict(self, s_iJ, d_params):
        """
        Load a single element (by string index) from d_params.

        Args:
            s_iJ (str): String index (e.g., 'x1').
            d_params (dict): Dictionary of parameter values.
        """

        elt = f'{self._param_nm}_{s_iJ}'
        e_i, e_J = self._siJ_to_rc(s_iJ)

        self.mat[e_i, e_J] =  d_params[elt]

    def elt_iJ(self, s_iJ):
        """
        Return the value of a single element (by string index).

        Args:
            s_iJ (str): String index (e.g., 'x1').

        Returns:
            float: Value of the element.
        """

        e_i, e_J = self._siJ_to_rc(s_iJ)
        return self.mat[e_i, e_J]

    def set_elt_iJ(self, s_iJ, val):
        """
        Set the value of a single element (by string index).

        Args:
            s_iJ (str): String index (e.g., 'x1').
            val (float): Value to set.
        """

        e_i, e_J = self._siJ_to_rc(s_iJ)
        self.mat[e_i, e_J] = val

    def set_from_structure_elts(self, d_params, elts_indep, elts_dep, depwhich):
        """
        Set elements from known sets of independent and dependent values, with optional multipliers.

        Args:
            d_params (dict): Dictionary of parameter values.
            elts_indep (list): List of independent element indices.
            elts_dep (dict): Dictionary of dependent element indices and their relations.
            depwhich (int): Index of the multiplier to use.
        """

        for elt in elts_indep:
            self.set_elt_from_param_dict(elt, d_params)

        for iJ,v in elts_dep.items():
            s_iJ = v[0]
            fac = v[depwhich]
            val =self.elt_iJ(s_iJ)*fac
            self.set_elt_iJ(iJ, val)

    def value(self):
        """
        Return a copy of the tensor's matrix data (shape [3,6]).
        """

        return self.mat.copy()

    def refvalue(self):
        """
        Return a reference to the tensor's matrix data (shape [3,6]).
        """

        return self.mat[:, :]

    def copy(self):
        """
        Return a deep copy of the tensor object.
        """

        return copy.deepcopy(self)

class VoigtTensor4_IJ(object):
    '''A class for representing rank 4 tensors in the compact Voigt representation.

    Internally uses 1-based indexing like the notation. Externally it passes zero-based values.'''

    def __init__(self, param_nm, physical_name='', physical_symbol='', unit=None):
        """
        Initialize a VoigtTensor4_IJ object.

        Args:
            src_dict (dict): Source dictionary for parameter values.
            param_nm (str): Symbol head in parameter dictionary (e.g., 'c', 'p', 'eta').
            physical_name (str): Human-readable name for the tensor.
            physical_symbol (str): Symbol for display.
            unit (tuple or None): Unit information (name, scale).
        """

        self.mat = np.zeros([7, 7], dtype=float)  # unit indexing

        self._physical_name = physical_name
        self._param_nm = param_nm  # eg 'c', 'p', 'eta'

        #self._d_params = src_dict
        self._physical_symbol=physical_symbol
        self.unit=unit

    def as_str(self, chop=False, rtol=1e-10, atol=1e-12):
        """
        Return a string representation of the tensor, formatted for display.

        Args:
            chop (bool): Whether to chop small values to zero.
            rtol (float): Relative tolerance for chopping.
            atol (float): Absolute tolerance for chopping.

        Returns:
            str: Formatted string representation of the tensor.
        """

        sh = f'\n {self._physical_name} {self._physical_symbol}'
        if self.unit is not None:
            sh += f', unit: {self.unit[0]}.'
        else:
            sh += '., unit: dimensionless.'

        sh += '\n   Voigt 4-tensor:\n'
        with np.printoptions(precision=4):
            if self.unit is not None:
                sd = str(nbtools.chopmat(self.mat[1:, 1:]/self.unit[1], chop=chop, rtol=rtol, atol=atol))
            else:
                sd = str(nbtools.chopmat(self.mat[1:, 1:], chop=chop, rtol=rtol, atol=atol))

        return sh + nbtools.indent_string(sd, indent=4)

    def __str__(self):
        """
        Return a formatted string representation of the tensor, including units and values.
        """

        return self.as_str(chop=False)


    # Allow direct indexing of Voigt tensor in [(i,j)] form

    def __getitem__(self, k):
        """
        Get an element of the tensor using (i, j) tuple indexing.

        Args:
            k (tuple): Tuple of indices (i, j).

        Returns:
            float: Value at the specified indices.
        """

        return self.mat[k[0], k[1]]

    def __setitem__(self, k, v):
        """
        Set an element of the tensor using (i, j) tuple indexing.

        Args:
            k (tuple): Tuple of indices (i, j).
            v (float): Value to set.
        """

        self.mat[k[0], k[1]] = v

    def fill_random(self):
        """
        Fill the tensor with random values (uniform in [0,1)).
        """

        self.mat = nprand.rand(7,7)

    def elt_s_IJ(self, s_IJ):
        """
        Get an element by string index (e.g., '12', '23').

        Args:
            s_IJ (str): String index (e.g., '12').

        Returns:
            float: Value at the specified indices.
        """

        eI=int(s_IJ[0])
        eJ=int(s_IJ[1])
        assert(eI>=1 and eI<=6 and eJ>=1 and eJ<=6)
        return self.mat[eI,eJ]

    def set_elt_s_IJ(self, s_IJ, val):
        """
        Set an element by string index (e.g., '12', '23').

        Args:
            s_IJ (str): String index (e.g., '12').
            val (float): Value to set.
        """

        eI=int(s_IJ[0])
        eJ=int(s_IJ[1])
        assert(eI>=1 and eI<=6 and eJ>=1 and eJ<=6)
        self.mat[eI,eJ] = val

    def dump_rawdata(self):
        """
        Print the raw data of the tensor for debugging purposes.
        """

        print(f'\nVoigt tensor {self._physical_name}, tensor {self._param_nm}')
        print(self.mat[1:, 1:])

    def value(self):
        """
        Return a copy of the tensor's matrix data (shape [6,6]).
        """

        return self.mat[1:, 1:].copy()

    def refvalue(self):
        """
        Return a reference to the tensor's matrix data (shape [6,6]).
        """

        return self.mat[1:, 1:]


    def set_from_00matrix(self, mat00):
        """
        Set the tensor values from a 6x6 matrix.

        Args:
            mat00 (np.ndarray): 6x6 matrix.
        """

        if self.mat.dtype != mat00.dtype:
            self.mat = np.zeros([7, 7], dtype=mat00.dtype)
        self.mat[1:7, 1:7] = mat00

    def copy(self):
        """
        Return a deep copy of the tensor object.
        """

        return copy.deepcopy(self)


    def elt_specified(self, d_params, m, n):
        """Check if a tensor element is specified in the parameter dictionary."""

        elt = f'{self._param_nm}_{m}{n}'

        return elt in d_params

    def sIJ_to_rc(self, s_IJ):
        """
        Convert a string index (e.g., '12') to row and column indices.

        Args:
            s_IJ (str): String index (e.g., '12').

        Returns:
            tuple: (row, column) indices.
        """

        e_I = int(s_IJ[0])
        e_J = int(s_IJ[1])
        return e_I, e_J

    def set_elt_from_param_dict(self, I, J, d_params):
        """
        Load a single element (by row and column indices) from d_params.

        Args:
            I (int): Row index (1-based).
            J (int): Column index (1-based).
            d_params (dict): Dictionary of parameter values.
        """

        elt = f'{self._param_nm}_{I}{J}'

        self.mat[I, J] = d_params[elt]

    def set_elt_from_param_dict_as_str(self, s_IJ, d_params):
        """
        Load a single element (by string index) from d_params.

        Args:
            s_IJ (str): String index (e.g., '12').
            d_params (dict): Dictionary of parameter values.
        """

        elt = f'{self._param_nm}_{s_IJ}'
        e_i, e_J = self._siJ_to_rc(s_IJ)

        self.mat[e_i, e_J] =  d_params[elt]



    # def read_from_json(self, m, n, optional=False):
    #     """
    #     Read a tensor element from the parameter dictionary using indices (m, n).

    #     Args:
    #         m (int): Row index (1-based).
    #         n (int): Column index (1-based).
    #         optional (bool): If True, do not raise error if element is missing.

    #     Returns:
    #         bool: True if element was found and set, False otherwise.
    #     """

    #     elt = f'{self._param_nm}_{m}{n}'

    #     if elt not in self._d_params:
    #         if not optional:
    #             reporting.report_and_exit(
    #                 f'Failed to read required tensor element {elt} for material {self._physical_name}')
    #         else:
    #             return False

    #     self.mat[m, n] = self._d_params[elt]
    #     return True

    def load_isotropic_from_json(self, d_params):
        """
        Load isotropic tensor values from the parameter dictionary.
        """

        self.set_elt_from_param_dict(1, 1, d_params)
        self.set_elt_from_param_dict(1, 2, d_params)
        self.set_elt_from_param_dict(4, 4, d_params)
        self.make_isotropic_tensor(self.mat[1, 1], self.mat[1, 2], self.mat[4, 4])

    def make_isotropic_tensor(self, m11, m12, m44):
        """
        Build Voigt matrix from 3 parameters for isotropic geometry.
        (Actually, only two are independent.)

        Args:
            m11 (float): Diagonal value.
            m12 (float): Off-diagonal value.
            m44 (float): Shear value.
        """

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
        """
        Check that the tensor matrix is symmetric and positive definite.
        """

        # Check matrix is symmetric and positive definite

        rtol = 1e-12
        tol = rtol * np.abs(self.mat).max()
        tmat = self.mat - self.mat.T
        mat_is_sym = nbtools.almost_zero(np.linalg.norm(tmat), tol)
        reporting.assertion(
            mat_is_sym, f'Material matrix {self._physical_name}-{self._physical_symbol} is symmetric.\n' + str(self.mat))

    # def _old_rotate(self, mat_R):
    #     """
    #     Rotate the tensor using the SO(3) matrix mat_R.

    #     Args:
    #         mat_R (np.ndarray): 3x3 rotation matrix.
    #     """

    #     rot_tensor = _rotate_Voigt_4tensor(self.mat[1:, 1:], mat_R)
    #     self.mat[1:, 1:] = rot_tensor

    def rotate(self, mat_R):
        """
        Rotate the tensor using the SO(3) matrix mat_R and return a new tensor.

        Args:
            mat_R (np.ndarray): 3x3 rotation matrix.

        """

        mat_M = generate_voigt_rotation_matrix(mat_R)

        self.mat[1:, 1:] = mat_M @ self.mat[1:, 1:] @ mat_M.T


    def set_from_structure_elts(self, d_params, elts_indep, elts_dep):

        for s_IJ in elts_indep:
            eI=int(s_IJ[0])
            eJ=int(s_IJ[1])
            self.set_elt_from_param_dict(eI, eJ, d_params)

        for d_IJ,src in elts_dep.items():
            #if d_IJ == "66":
            #    print('setting elts ', d_IJ, 'from', src)
            dI=int(d_IJ[0])
            dJ=int(d_IJ[1])

            if self.elt_specified(d_params, dI, dJ):
                self.set_elt_from_param_dict(dI, dJ, d_params)
             #   print('found element ', d_IJ, 'in d_params, setting it directly to', self.mat[dI, dJ])
            else:
                for (s_IJ,sf) in src:
                    sI=int(s_IJ[0])
                    sJ=int(s_IJ[1])
                    self.mat[dI, dJ] += self.mat[sI,sJ] * sf

        #print('set Voigt tensor from structure elements, result:\n', self.mat[1:, 1:])