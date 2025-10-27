from typing import List, Tuple, Optional
from numpy.typing import NDArray
import numpy as np

def qvec3_to_symgrad(vecq: NDArray[np.complex128]) -> NDArray[np.complex128]:
    qx, qy, qz = vecq
    mSymGrad = 1j * np.array(
        [[qx, 0, 0], [0, qy, 0], [0, 0, qz], [0, qz, qy], [qz, 0, qx], [qy, qx, 0]],
        dtype=np.complex128,
    )
    return mSymGrad


def uvec3_to_strain_S6(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128]) -> NDArray[np.complex128]:
    vS = qvec3_to_symgrad(vecq) @ vecu
    return vS


def strain_S6_to_stress_T6(cstiff: NDArray[np.float64], vecS: NDArray[np.complex128]) -> NDArray[np.complex128]:
    vT6 = cstiff @ vecS
    return vT6


def uvec3_to_stress_T6(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.complex128]:
    vT6 = cstiff @ qvec3_to_symgrad(vecq) @ vecu
    return vT6


def uvec3_to_stress_T33(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.complex128]:
    vT6 = uvec3_to_stress_T6(vecu, vecq, cstiff)
    return np.array(
        [[vT6[0], vT6[5], vT6[4]], [vT6[5], vT6[1], vT6[3]], [vT6[4], vT6[3], vT6[2]]],
        dtype=np.complex128,
    )


def T6_to_T33(vecT6: NDArray[np.complex128]) -> NDArray[np.complex128]:
    T33 = np.array(
        [[vecT6[0], vecT6[5], vecT6[4]],
         [vecT6[5], vecT6[1], vecT6[3]],
         [vecT6[4], vecT6[3], vecT6[2]]],
        dtype=np.complex128,
    )
    return T33

# This is slow, inline
def uvec3_to_poynting_S3(Om: float, vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.float64]:
    T33 = uvec3_to_stress_T33(vecu, vecq, cstiff)
    vecv = -1j * Om * vecu
    # Poynting vector is 0.5 Re(-T . v*)
    Sav = np.real(-T33 @ vecv.conj()) * 0.5
    return Sav

# This is slow, inline
def uvec3_phi_to_poynting_S3_piezo(Om: float, vecu: NDArray[np.complex128], vecq: NDArray[np.complex128],
                                    phi: np.complex128, cstiff: NDArray[np.float64],
                                    eiJ: NDArray[np.float64], permij: NDArray[np.float64]) -> NDArray[np.float64]:
    S_J = uvec3_to_strain_S6(vecu, vecq)
    vecE = -1j * vecq * phi

    T6 = strain_S6_to_stress_T6(cstiff, S_J) -  vecE @ eiJ
    T33 = T6_to_T33(T6)

    D3 = eiJ @ S_J + permij @ vecE

    vecv = -1j * Om * vecu

    Sav = .5 * np.real(-T33 @ vecv.conj() + 1j * Om * phi* D3.conj())
    return Sav