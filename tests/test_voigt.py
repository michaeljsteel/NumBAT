
import sys
sys.path.append('backend')

import numpy as np
#import voigt


from voigt import PlainTensor2_ij, VoigtTensor3_iJ, VoigtTensor4_IJ, make_rotation_matrix
#from voigt import _rotate_2tensor_b, _rotate_2tensor

twopi = 2 * np.pi

# def test_plain_mat3x3_rotation():
#     mat = np.random.rand(3, 3)
#     mat_R = make_rotation_matrix((1, 0, 0), twopi / 4)  # Rotate around x-axis by 90 degrees

#     mat_rot1 = _rotate_2tensor_b(mat, mat_R)
#     mat_rot2 = _rotate_2tensor(mat, mat_R)

#     np.testing.assert_allclose(mat_rot1, mat_rot2, atol=1e-12)


def test_plain_tensor2_ij_rotation():

    tens = PlainTensor2_ij('testsym')
    tens.fill_random()

    for nvec in ((1,0,0),(0,1,0),(0,0,1),
        (1,1,0),(1,0,1),(0,1,1), (1,1,1)):

        for nrots in range(1,6):
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()

            for i in range(nrots):
                tens_rot.rotate(mat_R)
            np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)


def test_voigt_tensor3_iJ_rotation():

    tens = VoigtTensor3_iJ('testsym')
    tens.fill_random()

    for nvec in ((1,0,0),(0,1,0),(0,0,1),
        (1,1,0),(1,0,1),(0,1,1), (1,1,1)):
        for nrots in range(1,6):
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()

            for i in range(nrots):
                tens_rot.rotate(mat_R)
            np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)


def test_voigt_tensor4_IJ_rotation():

    tens = VoigtTensor4_IJ('testsym', None)
    tens.fill_random()

    for nvec in ((1,0,0),(0,1,0),(0,0,1),
        (1,1,0),(1,0,1),(0,1,1), (1,1,1)):
        for nrots in range(1,6):
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()

            for i in range(nrots):
                tens_rot.rotate(mat_R)
            np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)

# def test_voigt_tensor4_IJ_rotation_via_M6x6():
#     tens = VoigtTensor4_IJ('testsym', None)
#     tens.fill_random()
#     tens2 = tens.copy()

#     for nvec in ((1,0,0),(0,1,0),(0,0,1),
#         (1,1,0),(1,0,1),(0,1,1), (1,1,1)):l
#         for nrots in range(1,6):
#             mat_R = make_rotation_matrix(nvec, twopi/nrots)
#             tens.rotate(mat_R)
#             tens2.rotate_b(mat_R)

#             np.testing.assert_allclose(tens.mat, tens2.mat, atol=1e-12)

def test_voigt_tensor3_iJ_rotation_via_M6x6():
    tens = VoigtTensor3_iJ('testsym', None)
    tens.fill_random()
    tens2 = tens.copy()

    for nvec in ((1,0,0),(0,1,0),(0,0,1),
        (1,1,0),(1,0,1),(0,1,1), (1,1,1)):
        for nrots in range(1,6):
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens.rotate(mat_R)
            tens2.rotate_b(mat_R)

            np.testing.assert_allclose(tens.mat, tens2.mat, atol=1e-12)