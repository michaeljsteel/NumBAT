
import math
import numpy as np

import reporting
import numbattools
from nbtypes import unit_x, unit_y, unit_z

# Array that converts between 4th rank tensors in terms of x,y,z and Voigt notation
#               [[xx,xy,xz], [yx,yy,yz], [zx,zy,zz]]
to_Voigt = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])


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

    nabla_IJ = np.array([
        [kx,   0.0, 0.0],
        [0.0,  ky,  0.0],
        [0.0,  0.0, kz],
        [0,    kz,  ky],
        [kz,   0.0, kx],
        [ky,   kx, 0.0]
    ])
    return nabla_IJ




def rotate_tensor_elt(i, j, k, l, T_pqrs, mat_R):
    '''
    Calculates the element ijkl of the rotated tensor Tp from the original
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
    if numbattools.almost_zero(nvec):
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

    reporting.assertion(numbattools.almost_unity(
        np.linalg.det(mat_R)), 'Rotation matrix has unit determinant.')

    return mat_R


def rotate_3vector(vec3, mat_R):


    #vrot = 0*vec
    #for i in range(3):
    #    vrot[i] = mat_R[i, 0]*vec[0] + mat_R[i, 1]*vec[1] + mat_R[i, 2]*vec[2]
    #return vrot

    return np.matmul(mat_R, vec3)



def _rotate_Voigt_tensor(T_PQ, mat_R):
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

    # mat_R = _make_rotation_matrix(theta, rotation_axis)

    Tp_PQ = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            V1 = to_Voigt[i, j]
            for k in range(3):
                for l in range(3):
                    V2 = to_Voigt[k, l]
                    Tp_PQ[V1, V2] = rotate_tensor_elt(i, j, k, l, T_PQ, mat_R)

    return Tp_PQ



class VoigtTensor4(object):
    '''A class for representing rank 4 tensors in the compact Voigt representation.

    Internally uses 1-based indexing like the notation. Externally it passes zero-based values.'''

    def __init__(self, material_name, symbol, src_dict=None, physical_name='', unit=None):

        self.mat = np.zeros([7, 7], dtype=float)  # unit indexing

        self.material_name = material_name
        self.symbol = symbol  # eg 'c', 'p', 'eta'

        self._json_data = src_dict
        self._physical_name=physical_name
        self._unit=unit

    # Allow direct indexing of Voigt tensor in [(i,j)] form

    def __getitem__(self, k):
        return self.mat[k[0], k[1]]

    def __setitem__(self, k, v):
        self.mat[k[0], k[1]] = v

    def __str__(self):

        with np.printoptions(precision=4):
            s = f'\nVoigt tensor {self.material_name}, {self._physical_name} {self.symbol}'
            if self._unit is not None:
                s+=f', unit: {self._unit[0]}. \n'
                s += str(self.mat[1:, 1:]/self._unit[1])
            else:
                s+='.\n'
                s += str(self.mat[1:, 1:])

        return s

    def dump_rawdata(self):
        print(f'\nVoigt tensor {self.material_name}, tensor {self.symbol}')
        print(self.mat[1:, 1:])

    def value(self):
        '''Returns copy of Voigt matrix indexed as m[0..5, 0..5].'''
        return self.mat[1:, 1:].copy()

    def refvalue(self):
        '''Returns reference to internal Voigt matrix indexed as m[0..5, 0..5].'''
        return self.mat[1:, 1:]

    def read_from_json(self, m, n, optional=False):
        '''Looks for data in the _json_data dict in form c_12, eta_23, etc.'''

        elt = f'{self.symbol}_{m}{n}'

        if elt not in self._json_data:
            if not optional:
                reporting.report_and_exit(
                    f'Failed to read required tensor element {elt} for material {self.material_name}')
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
        mat_is_sym = numbattools.almost_zero(np.linalg.norm(tmat), tol)
        reporting.assertion(
            mat_is_sym, f'Material matrix {self.material_name}-{self.symbol} is symmetric.\n' + str(self.mat))

    def rotate(self, matR):
        '''Rotates the crystal according to the SO(3) matrix matR.
        '''

        rot_tensor = _rotate_Voigt_tensor(self.mat[1:, 1:], matR)
        self.mat[1:, 1:] = rot_tensor
