"""Conversions between elastic displacement, strain, stress, and energy flux.

This module provides functions to convert between different representations
of elastic fields in Fourier space, including displacement to strain,
strain to stress, and computation of acoustic Poynting vectors.

All functions work with Voigt notation for tensors where applicable.
"""

from numpy.typing import NDArray
import numpy as np

def qvec3_to_symgrad(vecq: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Construct symmetric gradient operator matrix from wavevector.

    Builds the 6x3 matrix representation of the symmetric gradient operator
    in Voigt notation: symgrad(u) = (1/2)(∇u + ∇uᵀ) expressed in terms of
    the wavevector q for plane wave fields u(r) = û exp(iq·r).

    Parameters
    ----------
    vecq : NDArray[np.complex128]
        Wavevector components [qx, qy, qz] in 1/m.

    Returns
    -------
    NDArray[np.complex128]
        6x3 symmetric gradient operator matrix in Voigt notation.
        Rows correspond to strain components: [xx, yy, zz, yz, xz, xy].
    """
    qx, qy, qz = vecq
    mSymGrad = 1j * np.array(
        [[qx, 0, 0], [0, qy, 0], [0, 0, qz], [0, qz, qy], [qz, 0, qx], [qy, qx, 0]],
        dtype=np.complex128,
    )
    return mSymGrad


def uvec3_to_strain_S6(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Convert displacement amplitude to strain tensor in Voigt notation.

    Computes the strain tensor S from displacement u via S = (1/2)(∇u + ∇uᵀ)
    for a plane wave field u(r) = vecu exp(iq·r).

    Parameters
    ----------
    vecu : NDArray[np.complex128]
        Displacement amplitude vector [ux, uy, uz] in meters.
    vecq : NDArray[np.complex128]
        Wavevector components [qx, qy, qz] in 1/m.

    Returns
    -------
    NDArray[np.complex128]
        Strain components in Voigt notation [Sxx, Syy, Szz, Syz, Sxz, Sxy].
    """
    vS = qvec3_to_symgrad(vecq) @ vecu
    return vS


def strain_S6_to_stress_T6(cstiff: NDArray[np.float64], vecS: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Convert strain to stress using the stiffness tensor.

    Computes stress T from strain S via Hooke's law: T = C:S where C is
    the elastic stiffness tensor in Voigt notation.

    Parameters
    ----------
    cstiff : NDArray[np.float64]
        6x6 elastic stiffness tensor in Voigt notation (Pa).
    vecS : NDArray[np.complex128]
        Strain components in Voigt notation [Sxx, Syy, Szz, Syz, Sxz, Sxy].

    Returns
    -------
    NDArray[np.complex128]
        Stress components in Voigt notation [Txx, Tyy, Tzz, Tyz, Txz, Txy] (Pa).
    """
    vT6 = cstiff @ vecS
    return vT6


def uvec3_to_stress_T6(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.complex128]:
    """Convert displacement amplitude directly to stress tensor in Voigt notation.

    Computes stress from displacement in one step: T = C:S where S = symgrad(u).

    Parameters
    ----------
    vecu : NDArray[np.complex128]
        Displacement amplitude vector [ux, uy, uz] in meters.
    vecq : NDArray[np.complex128]
        Wavevector components [qx, qy, qz] in 1/m.
    cstiff : NDArray[np.float64]
        6x6 elastic stiffness tensor in Voigt notation (Pa).

    Returns
    -------
    NDArray[np.complex128]
        Stress components in Voigt notation [Txx, Tyy, Tzz, Tyz, Txz, Txy] (Pa).
    """
    vT6 = cstiff @ qvec3_to_symgrad(vecq) @ vecu
    return vT6


def uvec3_to_stress_T33(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.complex128]:
    """Convert displacement amplitude to stress tensor in full 3x3 matrix form.

    Computes the stress tensor T and returns it as a 3x3 symmetric matrix
    rather than in Voigt notation.

    Parameters
    ----------
    vecu : NDArray[np.complex128]
        Displacement amplitude vector [ux, uy, uz] in meters.
    vecq : NDArray[np.complex128]
        Wavevector components [qx, qy, qz] in 1/m.
    cstiff : NDArray[np.float64]
        6x6 elastic stiffness tensor in Voigt notation (Pa).

    Returns
    -------
    NDArray[np.complex128]
        3x3 stress tensor matrix (Pa).
    """
    vT6 = uvec3_to_stress_T6(vecu, vecq, cstiff)
    return np.array(
        [[vT6[0], vT6[5], vT6[4]], [vT6[5], vT6[1], vT6[3]], [vT6[4], vT6[3], vT6[2]]],
        dtype=np.complex128,
    )


def T6_to_T33(vecT6: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Convert stress from Voigt notation (6-vector) to 3x3 matrix form.

    Parameters
    ----------
    vecT6 : NDArray[np.complex128]
        Stress components in Voigt notation [Txx, Tyy, Tzz, Tyz, Txz, Txy] (Pa).

    Returns
    -------
    NDArray[np.complex128]
        3x3 symmetric stress tensor matrix (Pa).
    """
    T33 = np.array(
        [[vecT6[0], vecT6[5], vecT6[4]],
         [vecT6[5], vecT6[1], vecT6[3]],
         [vecT6[4], vecT6[3], vecT6[2]]],
        dtype=np.complex128,
    )
    return T33

# This is slow, inline
def uvec3_to_poynting_S3_slow(Om: float, vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.float64]:
    """Compute time-averaged acoustic Poynting vector (reference implementation).

    Computes the acoustic energy flux density S = (1/2)Re(-T·v*) where T is
    stress, v is velocity, and * denotes complex conjugate.

    Note: This is a reference implementation. Use uvec3_to_poynting_S3() for
    better performance.

    Parameters
    ----------
    Om : float
        Angular frequency (rad/s).
    vecu : NDArray[np.complex128]
        Displacement amplitude vector [ux, uy, uz] in meters.
    vecq : NDArray[np.complex128]
        Wavevector components [qx, qy, qz] in 1/m.
    cstiff : NDArray[np.float64]
        6x6 elastic stiffness tensor in Voigt notation (Pa).

    Returns
    -------
    NDArray[np.float64]
        Time-averaged Poynting vector [Sx, Sy, Sz] in W/m².
    """
    T33 = uvec3_to_stress_T33(vecu, vecq, cstiff)
    vecv = -1j * Om * vecu
    # Poynting vector is 0.5 Re(-T . v*)
    Sav = np.real(-T33 @ vecv.conj()) * 0.5
    return Sav

# This is an expensive function for some reason. 17% of total runtime.
# This is inlined version of function immediately above.
def uvec3_to_poynting_S3(Om: float, vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.float64]:
    """Compute time-averaged acoustic Poynting vector (optimized).

    Computes the acoustic energy flux density S = (1/2)Re(-T·v*) where T is
    stress and v is velocity. This optimized version inlines all intermediate
    operations to minimize memory allocations and function call overhead.

    The Poynting vector represents the time-averaged energy flux carried by
    the acoustic wave.

    Parameters
    ----------
    Om : float
        Angular frequency (rad/s).
    vecu : NDArray[np.complex128]
        Displacement amplitude vector [ux, uy, uz] in meters.
    vecq : NDArray[np.complex128]
        Wavevector components [qx, qy, qz] in 1/m.
    cstiff : NDArray[np.float64]
        6x6 elastic stiffness tensor in Voigt notation (Pa).

    Returns
    -------
    NDArray[np.float64]
        Time-averaged Poynting vector [Sx, Sy, Sz] in W/m².
    """

    qx, qy, qz = vecq
    ux, uy, uz = vecu

    symgrad_u = 1j * np.array([ux * qx, uy * qy, uz * qz,
                          uy * qz + uz * qy,
                          uz * qx + ux * qz,
                          ux * qy + uy * qx], dtype=np.complex128)
    vT6 = cstiff @ symgrad_u

    T33 = np.array(
        [[vT6[0], vT6[5], vT6[4]], [vT6[5], vT6[1], vT6[3]], [vT6[4], vT6[3], vT6[2]]],
        dtype=np.complex128,
    )

    # Poynting vector is 0.5 Re(-T . v*) = 0.5 Re(-T . ( -1j Om u)* ) = 0.5j Re( Om T . u* )
    vecv = 0.5j * Om * vecu
    Sav = np.real(T33 @ vecv.conj())
    return Sav




# This is slow, inline
def uvec3_phi_to_poynting_S3_piezo(Om: float, vecu: NDArray[np.complex128], vecq: NDArray[np.complex128],
                                    phi: np.complex128, cstiff: NDArray[np.float64],
                                    eiJ: NDArray[np.float64], permij: NDArray[np.float64]) -> NDArray[np.float64]:
    """Compute time-averaged Poynting vector for piezoelectric materials.

    Computes the acoustic energy flux including piezoelectric coupling:
    S = (1/2)Re(-T·v* + iω φ* D) where T is stress, v is velocity,
    φ is electric potential, D is electric displacement, and ω is frequency.

    The piezoelectric stress is T = C:S - e·E where e is the piezoelectric
    coupling tensor and E is the electric field.

    Parameters
    ----------
    Om : float
        Angular frequency (rad/s).
    vecu : NDArray[np.complex128]
        Displacement amplitude vector [ux, uy, uz] in meters.
    vecq : NDArray[np.complex128]
        Wavevector components [qx, qy, qz] in 1/m.
    phi : np.complex128
        Electric potential amplitude in volts.
    cstiff : NDArray[np.float64]
        6x6 elastic stiffness tensor in Voigt notation (Pa).
    eiJ : NDArray[np.float64]
        3x6 piezoelectric coupling tensor (C/m²).
    permij : NDArray[np.float64]
        3x3 dielectric permittivity tensor (F/m).

    Returns
    -------
    NDArray[np.float64]
        Time-averaged Poynting vector [Sx, Sy, Sz] in W/m².
    """
    S_J = uvec3_to_strain_S6(vecu, vecq)
    vecE = -1j * vecq * phi

    T6 = strain_S6_to_stress_T6(cstiff, S_J) -  vecE @ eiJ
    T33 = T6_to_T33(T6)

    D3 = eiJ @ S_J + permij @ vecE

    vecv = -1j * Om * vecu

    Sav = .5 * np.real(-T33 @ vecv.conj() + 1j * Om * phi* D3.conj())
    return Sav