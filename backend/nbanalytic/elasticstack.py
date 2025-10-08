import numpy as np
import numpy.linalg as npla
import scipy.linalg as spla
from typing import List, Tuple, Optional
from numpy.typing import NDArray

from enum import Enum

import plotting.plottools as plpl
import matplotlib.pyplot as plt

import nbtypes
import materials
from numtools.optimize import bracket_minima, find_bracketed_minima, BracketError
from nbtypes import SI_GHz, SI_THz, SI_um, SI_kmps
import plotting.plottools as nbpt


class NormalisationConstants:
    def __init__(self):
        # adjustable
        self.t0 = nbtypes.SI_ns
        self.x0 = nbtypes.SI_um
        self.rho0 = 1000  # 1000 kg/m^3
        self.eps0 = nbtypes.SI_permittivity_eps0  # 8.854187817e-12  F/m
        self.V0 = 1.0  # V

        # derived
        self.f0 = 1 / self.t0  # GHz
        self.v0 = self.x0 / self.t0  # km/s

        self.T0 = self.rho0 * self.v0**2  # GPa
        self.c0 = self.T0  # GPa
        self.E0 = self.V0 / self.x0  # V/um
        self.D0 = self.eps0 * self.E0  # C/m^2

        self.e0 = 1e-2  # C/m^2
        # self.p0 = nbtypes.SI_GPa


g_norms = NormalisationConstants()


class ElBCType(Enum):
    VACUUM = "vacuum"
    SEMI_INFINITE = "semiinfinite"
    SHORT = "short"
    CHARGE_FREE = "chargefree"


def qvec_to_symgrad(vecq: NDArray[np.complex128]) -> NDArray[np.complex128]:
    qx, qy, qz = vecq
    mSymGrad = 1j * np.array(
        [[qx, 0, 0], [0, qy, 0], [0, 0, qz], [0, qz, qy], [qz, 0, qx], [qy, qx, 0]],
        dtype=np.complex128,
    )
    return mSymGrad


def uvec3_to_strain_S6(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128]) -> NDArray[np.complex128]:
    vS = qvec_to_symgrad(vecq) @ vecu
    return vS


def strain_S6_to_stress_T6(cstiff: NDArray[np.float64], vecS: NDArray[np.complex128]) -> NDArray[np.complex128]:
    vT6 = cstiff @ vecS
    return vT6


def uvec3_to_stress_T6(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.complex128]:
    vT6 = cstiff @ qvec_to_symgrad(vecq) @ vecu
    return vT6


def uvec3_to_stress_T33(vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.complex128]:
    vT6 = uvec3_to_stress_T6(vecu, vecq, cstiff)
    return np.array(
        [[vT6[0], vT6[5], vT6[4]], [vT6[5], vT6[1], vT6[3]], [vT6[4], vT6[3], vT6[2]]],
        dtype=np.complex128,
    )


def poynting_vector(Om: float, vecu: NDArray[np.complex128], vecq: NDArray[np.complex128], cstiff: NDArray[np.float64]) -> NDArray[np.float64]:
    T33 = uvec3_to_stress_T33(vecu, vecq, cstiff)
    vecv = -1j * Om * vecu
    # Poynting vector is 0.5 Re(-T . v*)
    Sav = np.real(-T33 @ vecv.conj()) * 0.5
    return Sav


# Coefficient matrices defined in the layers3.nb mathematica file
def get_layer_stiffness_submatrices(cs: NDArray[np.float64], use_4d: bool = False) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    # cs is normalised stiffness

    mM1 = np.array(
        [
            [cs[6, 6], cs[6, 2], cs[6, 4]],
            [cs[2, 6], cs[2, 2], cs[2, 4]],
            [cs[4, 6], cs[4, 2], cs[4, 4]],
        ]
    )

    mM2 = np.array(
        [
            [cs[6, 5], cs[6, 4], cs[6, 3]],
            [cs[2, 5], cs[2, 4], cs[2, 3]],
            [cs[4, 5], cs[4, 4], cs[4, 3]],
        ]
    )

    mL1 = np.array(
        [
            [cs[5, 5], cs[5, 4], cs[5, 3]],
            [cs[4, 5], cs[4, 4], cs[4, 3]],
            [cs[3, 5], cs[3, 4], cs[3, 3]],
        ]
    )

    mL2 = np.array(
        [
            [cs[5, 6], cs[5, 2], cs[5, 4]],
            [cs[4, 6], cs[4, 2], cs[4, 4]],
            [cs[3, 6], cs[3, 2], cs[3, 4]],
        ]
    )

    if use_4d:
        # modify for 8-direction notation
        mM1_4d = np.zeros(4, 4, dtype=np.float64)
        mM2_4d = np.zeros(4, 4, dtype=np.float64)
        mL1_4d = np.zeros(4, 4, dtype=np.float64)
        mL2_4d = np.zeros(4, 4, dtype=np.float64)
        mM1_4d[:3, :3] = mM1
        mM2_4d[:3, :3] = mM2
        mL1_4d[:3, :3] = mL1
        mL2_4d[:3, :3] = mL2

        mM1, mM2, mL1, mL2 = mM1_4d, mM2_4d, mL1_4d, mL2_4d

    return mM1, mM2, mL1, mL2


def get_layer_stiffness_submatrices_piezo(cs: NDArray[np.float64], e_iJ: NDArray[np.float64], perm_ij: NDArray[np.float64]) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    # tensors are all unit indexed

    mM1 = np.array(
        [
            [cs[6, 6], cs[6, 2], cs[6, 4], e_iJ[2, 6]],
            [cs[2, 6], cs[2, 2], cs[2, 4], e_iJ[2, 2]],
            [cs[4, 6], cs[4, 2], cs[4, 4], e_iJ[2, 4]],
            [e_iJ[2, 6], e_iJ[2, 2], e_iJ[2, 4], -perm_ij[2, 2]],
        ]
    )

    mM2 = np.array(
        [
            [cs[6, 5], cs[6, 4], cs[6, 3], e_iJ[3, 6]],
            [cs[2, 5], cs[2, 4], cs[2, 3], e_iJ[3, 2]],
            [cs[4, 5], cs[4, 4], cs[4, 3], e_iJ[3, 4]],
            [e_iJ[2, 5], e_iJ[2, 4], e_iJ[2, 3], -perm_ij[2, 3]],
        ]
    )

    mL1 = np.array(
        [
            [cs[5, 5], cs[5, 4], cs[5, 3], e_iJ[3, 5]],
            [cs[4, 5], cs[4, 4], cs[4, 3], e_iJ[3, 4]],
            [cs[3, 5], cs[3, 4], cs[3, 3], e_iJ[3, 3]],
            [e_iJ[3, 5], e_iJ[3, 4], e_iJ[3, 3], -perm_ij[3, 3]],
        ]
    )

    mL2 = np.array(
        [
            [cs[5, 6], cs[5, 2], cs[5, 4], e_iJ[2, 5]],
            [cs[4, 6], cs[4, 2], cs[4, 4], e_iJ[2, 4]],
            [cs[3, 6], cs[3, 2], cs[3, 4], e_iJ[2, 3]],
            [e_iJ[3, 6], e_iJ[3, 2], e_iJ[3, 4], -perm_ij[3, 2]],
        ]
    )

    return mM1, mM2, mL1, mL2


class ElasticLayer:
    def __init__(self, material, #: materials.Material,
                 L: float = 1e-6, use_4d: bool = False) -> None:
        # dimensionless thickness

        self.material = material
        self.L = L / g_norms.x0

        self._build_material_quantities(use_4d)

    def L_phys(self) -> float:
        return self.L * g_norms.x0

    def _build_material_quantities(self, use_4d: bool) -> None:

        self.cstiff_norm = (
            self.material.stiffness_c_IJ.unit_indexed_value() / g_norms.c0
        )
        self.rho_norm = self.material.rho / g_norms.rho0

        if self.material.piezo_active():
            self.eiJ_norm = (
                self.material._piezo._tens_e_iJ.unit_indexed_value() / g_norms.e0
            )
            self.permij_norm = (
                self.material._piezo._tens_releps_S_ij.unit_indexed_value()
                / g_norms.eps0
            )
        else:
            self.eiJ_norm = np.zeros((4, 7), dtype=np.float64)
            self.permij_norm = np.eye(4, dtype=np.float64)

        if self.material.piezo_active():
            mM1, mM2, mL1, mL2 = get_layer_stiffness_submatrices_piezo(
                self.cstiff_norm, self.eiJ_norm, self.permij_norm
            )
        else:
            mM1, mM2, mL1, mL2 = get_layer_stiffness_submatrices(
                self.cstiff_norm, use_4d
            )
        self.mM1, self.mM2, self.mL1, self.mL2 = mM1, mM2, mL1, mL2


    def find_layer_transfer_matrix(self, Om: float, Vp: float) -> Tuple[NDArray[np.complex128], Optional[NDArray[np.complex128]]]:
        # Om, Vp are both normalised

        mM1, mM2, mL1, mL2 = self.mM1, self.mM2, self.mL1, self.mL2

        matT = np.zeros([6, 6], dtype=np.complex128)

        # ddz [T] = [ M_TT  M_TU] [T]
        #     [U] = [ M_UT  M_UU] [U]

        mM1inv = npla.inv(mM1)

        matT[:3, :3] = -1j * Om / Vp * mL2 @ mM1inv  # M_TT
        matT[:3, 3:6] = (Om / Vp) ** 2 * (
            mL1 - self.rho_norm * Vp**2 * np.eye(3) - mL2 @ mM1inv @ mM2
        )  # M_TU

        matT[3:6, 0:3] = mM1inv  # M_UT
        matT[3:6, 3:6] = -1j * Om / Vp * mM1inv @ mM2  # M_UU

        # exponentiate to get layer transfer matrix
        matPdz = spla.expm(matT * self.L) if self.L > 0 else None

        return matT, matPdz

def is_growing_eigenvalue(ev: complex) -> bool:
    return not np.isclose(np.real(ev), 0) and np.real(ev) > 0

def is_decaying_eigenvalue(ev: complex) -> bool:
    return not np.isclose(np.real(ev), 0) and np.real(ev) < 0

class ElasticBoundaryCondition:
    def __init__(self, bc_type: ElBCType,
                 material = None
                 #material: Optional[materials.Material] = None
                 ) -> None:
        self.bc_type = bc_type
        self.material = material
        self.layer = ElasticLayer(material, 0.0) if material is not None else None
        if self.is_semi_infinite() and material is None:
            raise ValueError("Semi-infinite boundary condition requires a material.")

        self._semiinf_evals = None
        self._semiinf_evecs = None


    def is_vacuum(self) -> bool:
        return self.bc_type == ElBCType.VACUUM

    def is_semi_infinite(self) -> bool:
        return self.bc_type == ElBCType.SEMI_INFINITE

    def is_short(self) -> bool:
        return self.bc_type == ElBCType.SHORT

    def is_charge_free(self) -> bool:
        return self.bc_type == ElBCType.CHARGE_FREE

    def __str__(self) -> str:
        return f"ElasticBoundaryCondition({self.bc_type})"

    def analyse_semiinfinite_eigenspace(self, Om: float, Vp: float,
                                        growing_is_good: bool ) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        # find eigenvalues and eigenvectors of semi-infinite layer transfer matrix
        # identify growing and "forwards" modes and put them at front of the list

        mat_T, mat_expTL = self.layer.find_layer_transfer_matrix(Om, Vp)

        eigvals, eigvecs = npla.eig(mat_T)

        ev_good = []
        ev_bad = []
        ngrow = 0
        ndecay = 0
        nfwd = 0

        v_is_growing = np.zeros(len(eigvals), dtype=bool)
        v_is_decaying = np.zeros(len(eigvals), dtype=bool)
        v_is_fwds = np.zeros(len(eigvals), dtype=bool)
        v_is_pureimag = np.zeros(len(eigvals), dtype=bool)
        v_poynting_Sz = np.zeros((len(eigvals), 3), dtype=np.float64)

        for i in range(len(eigvals)):
            ev = eigvals[i]  # this is lambda = i q_y

            is_growing = is_growing_eigenvalue(ev)
            is_decaying = is_decaying_eigenvalue(ev)
            is_fwds = np.isclose(np.real(ev), 0) and np.imag(ev) > 0

            v_is_growing[i] = is_growing
            v_is_decaying[i] = is_decaying
            v_is_pureimag[i] = np.isclose(np.real(ev), 0)
            v_is_fwds[i] = is_fwds

            vecq = np.array([0, ev / 1j, Om / Vp], dtype=np.complex128)
            vecu = eigvecs[3:6, i]
            Svec = poynting_vector(Om, vecu, vecq, self.layer.cstiff_norm[1:,1:]) # cstiff is unit-indexed
            v_poynting_Sz[i, :] = Svec
            #is_fwds = Sz > 0


            if growing_is_good:
                if is_growing or (not is_decaying and is_fwds):
                    ev_good.append(i)
                else:
                    ev_bad.append(i)
            else:
                if is_decaying or (not is_growing and is_fwds):
                    ev_good.append(i)
                else:
                    ev_bad.append(i)

        ngrow = np.sum(v_is_growing)
        ndecay = np.sum(v_is_decaying)
        nfwd = np.sum(v_is_fwds)

        # sometimes after a branch cut (above shear velocity?) we can get 2 pure imaginary forward modes
        # for rayleigh, this may not matter, but for leaky modes it does
        # so perhaps we look for ones with Sy away from the interface?
        if len(ev_good) != 3:
            print('Unexpected number of good eigenvalues in semi-infinite layer:')
            print(f' Om:{Om:.4f}, Vp:{Vp:.4f}, (ngrow={ngrow}, ndecay={ndecay}, nfwd={nfwd}, ngood={len(ev_good)}).')
            for iev in v_is_growing.nonzero()[0]:
                Sx,Sy,Sz = v_poynting_Sz[iev,:] * 100
                print(f'  ev[{iev}] = {eigvals[iev]: 12.5f} grow,  Svec=({Sx:7.3f},{Sy:7.3f},{Sz:7.3f})')

            for iev in v_is_decaying.nonzero()[0]:
                Sx,Sy,Sz = v_poynting_Sz[iev,:]* 100
                print(f'  ev[{iev}] = {eigvals[iev]: 12.5f} decay, Svec=({Sx:7.3f},{Sy:7.3f},{Sz:7.3f})')

            for iev in v_is_pureimag.nonzero()[0]:
                Sx,Sy,Sz = v_poynting_Sz[iev,:]* 100
                print(f'  ev[{iev}] = {eigvals[iev]: 12.5f} imag,  Svec=({Sx:7.3f},{Sy:7.3f},{Sz:7.3f})')


            #raise ValueError(
            #    f"Semi-infinite layer: Expected 3 good eigenvalues, found {len(ev_good)} "
            #    )

        mVeigs = np.zeros((6, 6), dtype=np.complex128)
        ev_order = ev_good + ev_bad
        vEvals = eigvals[ev_order]
        mVeigs = eigvecs[:, ev_order]

        # print(f"Semi-infinite  layer: Kept evs", ev_good, eigvals[ev_good])
        self._semiinf_evals = vEvals
        self._semiinf_evecs = mVeigs

        return vEvals, mVeigs


    def suggested_L_left(self, Om: float, Vp: float) -> float:
        if not self.is_semi_infinite():
            return 0.0

        vEvals, mVeigs = self.analyse_semiinfinite_eigenspace(Om, Vp, growing_is_good=True)
        # find smallest decay length
        lamrs = np.real(vEvals[0:3])
        lammin =np.min(lamrs)
        if not np.isclose(lammin,0):
            return 6.0 / abs(lammin)  # 3 decay lengths
        else:
            return 5.0 # arbitrary


    def suggested_L_right(self, Om: float, Vp: float) -> float:
        return self.suggested_L_left(Om, Vp)


    def null_vector_to_field(self, v_null):
        v_psi0 = np.zeros(6, dtype=np.complex128)

        if self.is_vacuum():
            v_psi0[0:3] = 0.0
            v_psi0[3:6] = v_null  # null is literally initial displacement field
        elif self.is_semi_infinite():
            v_psi0 = self._semiinf_evecs[:, 0:3] @ v_null  # initial field from evectors
        else:
            raise NotImplementedError(
                "Boundary condition type not implemented for null vector to field conversion.")



        return v_psi0



class ElasticStack:
    def __init__(
        self, l_layers: List[ElasticLayer], bcs: List[ElasticBoundaryCondition]
    ):
        self._def_Vmin = 0.01  # per micron
        self._def_Vmax = 10.0  # per micron
        self._n_layers = len(l_layers)

        self._l_layers = l_layers
        self._bcs = bcs

        self.check_layers_and_bcs()

    def check_layers_and_bcs(self) -> None:
        if len(self._bcs) != 2:
            raise ValueError(
                "Boundary conditions list must have exactly two elements (top and bottom)."
            )


        if (self._bcs[0].is_vacuum() and self._bcs[1].is_vacuum() and self._n_layers == 0 ):
            raise ValueError(
                "At least one layer is required with vacuum boundary conditions on both sides."
            )

    def _extract_quadrant(self, mat_G: NDArray[np.complex128]) -> NDArray[np.complex128]:
        bc0_vaclike = not self._bcs[0].is_semi_infinite()
        bc1_vaclike = not self._bcs[1].is_semi_infinite()

        match (bc0_vaclike, bc1_vaclike):
            case (True, True):
                return mat_G[0:3, 3:6] # G12 = B
            case (True, False):
                return mat_G[3:6, 3:6] # G22 = D
            case (False, True):
                return mat_G[0:3, 0:3] # G11 = A
            case (False, False):
                return mat_G[3:6, 0:3] # G21 = C



    def make_layer_Pdz_mat(self):
        raise NotImplementedError(
            "Base class method called - should use subclass method."
        )

    def make_whole_stack_Pmat(self):
        raise NotImplementedError(
            "Base class method called - should use subclass method."
        )

    def dispersion_relation_find_Vs_at_Omegas(self, v_Omega, Vmin=0, Vmax=0):
        raise NotImplementedError(
            "Base class method called - should use subclass method."
        )

    def dispersion_relation_find_Vs_at_Omega(self, Omega, Vmin=0, Vmax=0):
        raise NotImplementedError(
            "Base class method called - should use subclass method."
        )

    def build_profile_domain_y(self, ny_per_layer: int = 50,
                               L_left: float = 0, L_right: float = 0) -> NDArray[np.float64]:
        # this is all in x0 units

        # get semi-inf borders
        if self._bcs[0].is_semi_infinite(): # add a bit of space for decay
            L_left = L_left if L_left > 0 else 5.0
        else:
            L_left = 0.0
        if self._bcs[1].is_semi_infinite(): # add a bit of space for decay
            L_right = L_right if L_right > 0 else 5.0
        else:
            L_right = 0.0


        #L_stack = sum(layer.L for layer in self._l_layers)
        #L_y = L_stack + L_left + L_right

        n_y = self._n_layers * ny_per_layer

        if L_left > 0:
            n_y += ny_per_layer
        if L_right > 0:
            n_y += ny_per_layer

        n_y +=1 # to include endpoint

        v_y = np.zeros(n_y)  # physical location
        v_ilayer = np.zeros(n_y, dtype=np.int32)  # index of layer (-1 for left semi-inf, -2 for right semi-inf)

        iyoff = 0
        ylo = 0
        yhi = 0
        if L_left > 0:
            ylo -= L_left
            dy = L_left / ny_per_layer
            v_y[iyoff:ny_per_layer] = np.linspace(ylo, yhi - dy, ny_per_layer)
            v_ilayer[iyoff:ny_per_layer] = -1
            iyoff += ny_per_layer
            ylo += L_left

        for ilay,layer in enumerate(self._l_layers):
            L = layer.L
            yhi += L
            dy = L / ny_per_layer
            v_y[iyoff : iyoff + ny_per_layer] = np.linspace(ylo, yhi - dy, ny_per_layer)
            v_ilayer[iyoff : iyoff + ny_per_layer] = ilay
            iyoff += ny_per_layer
            ylo += L

        if L_right > 0:
            yhi += L_right
            dy = L_right / ny_per_layer
            v_y[iyoff : iyoff + ny_per_layer] = np.linspace(ylo, yhi - dy, ny_per_layer)
            v_ilayer[iyoff : iyoff + ny_per_layer] = -2
            iyoff += ny_per_layer
            ylo += L_right

        v_y[-1] = v_y[-2] + dy
        v_ilayer[-1] = v_ilayer[-2]

        return v_y, v_ilayer

    def find_mode_profile_1d(self, Omega: float, Vs: float,
                             ny_per_layer: int = 50, L_left: float = -1, L_right: float = -1) -> NDArray[np.float64]:
        pass






class NonPiezoStack(ElasticStack):
    def __init__(self, l_layers: List[ElasticLayer], bcs: List[ElasticBoundaryCondition]) -> None:
        super().__init__(l_layers, bcs)

    # def eig_fwds(ev):
    #    val = np.isclose(np.real(ev), 0) and np.imag(ev) > 0
    #    return val

    def _find_full_transfer_matrix(self, Om: float, Vp: float) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        # build layer parts
        mat_G = np.eye(6, dtype=np.complex128)

        for layer in self._l_layers:
            mat_T, mat_expTL = layer.find_layer_transfer_matrix(Om, Vp)
            mat_G = mat_expTL @ mat_G


        # add boundary layer parts
        for ibcs, bc in enumerate(self._bcs):
            if bc.is_semi_infinite():

                growing_is_good = ibcs == 0
                # sort eigenvectors into growing and decaying/forwards and backwards modes
                v_evals, m_eigvecs = bc.analyse_semiinfinite_eigenspace(Om, Vp, growing_is_good)
                #print(f"BC {ibcs} Vp {Vp} eigvals:", npla.det(m_eigvecs))
                #print('BC eigmag:', ibcs,
                #      np.max(np.abs(m_eigvecs)),
                #      np.max(np.abs(npla.inv(m_eigvecs)))
                #      )

                if ibcs == 0:
                    mat_G = mat_G @ m_eigvecs
                else:
                    # Inverted eigenvectors can get large, but attempting to scale
                    # leads to unsmooth determinant scans
                    #m_inveigs = npla.inv(m_eigvecs)
                    #m_inveigs /= np.max(np.abs(m_inveigs))
                    #mat_G = m_inveigs @ mat_G

                    mat_G = npla.inv(m_eigvecs) @ mat_G

        mat_G_ij = self._extract_quadrant(mat_G)
        #print(f"detG Om {Om:.3f} Vp {Vp:.3f}:", npla.det(mat_G), npla.det(mat_G_ij))
        #print(mat_G)
        #print(mat_G_ij)
        #np.set_printoptions(precision=3)
        #print(mat_G_ij)
        mat_G_ij /= np.max(np.abs(mat_G_ij))  # scale to avoid overflow in det

        return mat_G, mat_G_ij

    def dispersion_relation_find_Vs_at_Omega(
        self,
        Omega_SI: float,
        Vmin: float = 10,
        Vmax: float = 10000.0,
        show_scan_plot: bool = False,
        show_progress: bool = False,
        fmin_thresh: float = 1e-5,
        find_mode_color: bool = False,
        prefix: str = "tmpdrscan",
    ) -> List[Tuple[float, float]]:


        # move to normalised units here
        Om = Omega_SI / g_norms.f0

        n_Vp = 500
        v_Vp = np.linspace(Vmin, Vmax, n_Vp) / g_norms.v0

        v_det = np.zeros(n_Vp, dtype=np.complex128)

        for ivP, vP in enumerate(v_Vp):
            mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
            v_det[ivP] = npla.det(mat_G_ij)
            #print(f"detG Om {Om:.3f} Vp {vP:.3f}:", npla.det(mat_G_ij))
            #print('detmag', ivP, np.max(np.abs(mat_G_ij)))

        if show_scan_plot:
            plot_det_scan(Om, v_Vp, v_det, prefix=prefix)


        l_braks = bracket_minima(v_Vp, np.abs(v_det))
        if show_progress:
            print(f"\nAt Nu={Om/(2*np.pi):.10f}, found {len(l_braks)} brackets:")

            for br in l_braks:
                (xm1, fm1), (x0, f0), (xp1, fp1) = np.array(br).tolist() # awkward but converts to plain floats
                print(f"    [{xm1:.5f}, {x0:.5f}, {xp1:.5f}] with f=[{fm1:.3e}, {f0:.3e}, {fp1:.3e}]")

        if len(l_braks) == 0:
            raise BracketError("No minima found in determinant scan - no modes found.")

        def minfunc(vP):
            mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
            return np.abs(npla.det(mat_G_ij))

        l_mins = find_bracketed_minima(minfunc, l_braks)
        if show_progress:
            print("  Found minima:")
            for v, f in l_mins:
                print(f"    Vp={v:.9f}, Det(f)={f:.4e}")

        # remove any minima with f above threshold or with multiple singular values as probably not real modes
        fmin_norm_thresh = fmin_thresh * np.max(np.abs(v_det))
        l_mins_SI_keep = []
        l_mins_bad_threshold = []
        l_mins_bad_singval = []

        for v, f in l_mins:
            if f > fmin_norm_thresh:
                l_mins_bad_threshold.append((v,f))
            elif self.count_singular_values(Om, v) != 1:
                l_mins_bad_singval.append((v,f))  # could be zero or two singvals
            else:
                l_mins_SI_keep.append((v * g_norms.v0,f))

        #l_mins_SI_keep = [(v * g_norms.v0, f) for v, f in l_mins if f < fmin_norm_thresh]
        #l_mins_SI_drop = [(v * g_norms.v0, f) for v, f in l_mins if f >= fmin_norm_thresh]

        n_dropped = len(l_mins_bad_threshold) + len(l_mins_bad_singval)
        if show_progress and n_dropped > 0:
            print(f"  Dropping {n_dropped} modes:")
            for v, f in l_mins_bad_singval:
                print(f"    Vp={v:.3f} m/s has multiple small singular values.")
            for v, f in l_mins_bad_threshold:
                print(f"    Vp={v:.3f} has relative det minimum f={f:.3e} above fmin_thresh={fmin_norm_thresh:.3e}:")

        # color modes by polarization at the surface.
        # better to do by energy weighted average over all domain but that is more expensive
        v_mode_cols = np.zeros((len(l_mins_SI_keep), 3), dtype=np.float64 )

        if find_mode_color:
            for iv, (vSI,f) in enumerate(l_mins_SI_keep):
                v_psi0 = self.find_mode_psi0(Om, vSI/g_norms.v0)
                #ux,uy,uz = np.abs(v_psi0[3:6])
                v_mode_cols[iv,:] = nbpt.get_rgb_for_poln(*v_psi0[3:6])

        return l_mins_SI_keep, v_mode_cols

    def count_singular_values(self, Om: float, vP: float,
                              thresh: float = 1e-6) -> int:
        mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
        u, s, vh = npla.svd(mat_G_ij)
        nsmall = np.sum(s < thresh*np.max(s)) # Limited contrast precision in SVD
        return nsmall

    def find_mode_psi0(self, Om: float, vP: float,
                       show_progress: bool = False) -> NDArray[np.complex128]:

        mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)

        u, s, vh = npla.svd(mat_G_ij)

        if show_progress:
            print('Det check:', npla.det(mat_G_ij))
            print('SVD check:', s)

        nsmall = np.sum(s < 1e-6*np.max(s)) # Limited contrast precision in SVD

        if show_progress:
            print('Num small singular values:', nsmall)

        if nsmall == 0:
            raise ValueError('No nullspace found - not a mode?')
        if nsmall > 1:
            raise ValueError('Multiple small singular values found - may be non-genuine mode.')


        v_null = vh.conj().T[:,-1]  # take the last column of Vh (smallest singular value)

        v_psi0 = self._bcs[0].null_vector_to_field(v_null)
        v_psi0 /= np.max(np.abs(v_psi0[3:6]))  # normalise max displacement to 1

        if show_progress:
            with np.printoptions(precision=4, suppress=True, linewidth=200):
                print(f'Initial field at y=0:\n null vec: {v_null}\n psi0:  {v_psi0}')

        return v_psi0

    def find_mode_profile_1d(self, Omega: float, Vs: float,
                             ny_per_layer: int = 100, L_left: float = -1, L_right: float = -1) -> NDArray[np.float64]:
        # find nullspace vector and initial field at y=0
        Om = Omega / g_norms.f0
        vP = Vs / g_norms.v0

        try:
            v_psi0 = self.find_mode_psi0(Om, vP)
        except ValueError as e:
            raise ValueError(f"Could not find mode at Omega={Omega:.3e} ({Om:.3f} norm), Vs={Vs:.3f} ({vP:.3f} norm): {e}") from e

        # find profile through stack

        # user can override default L_left and L_right for semi-infinite layers
        L_left = L_left if L_left > 0 else self._bcs[0].suggested_L_left(Om, vP)
        L_right = L_right if L_right > 0 else self._bcs[1].suggested_L_right(Om, vP)

        v_y, v_ilayer = self.build_profile_domain_y(ny_per_layer, L_left, L_right)

        # fill profile
        n_y = len(v_y)
        m_psi_y = np.zeros((6, n_y), dtype=np.complex128)

        # precompute layer transfer matrices
        has_left = v_ilayer[0] == -1
        has_right = v_ilayer[-1] == -2
        if has_left:
            mat_Tleft, mat_expTL = self._bcs[0].layer.find_layer_transfer_matrix(Om, vP)

        v_mat_Tlayer = []
        for ilay,layer in enumerate(self._l_layers):
            mat_Tlayer, mat_expTL = layer.find_layer_transfer_matrix(Om, vP)
            v_mat_Tlayer.append(mat_Tlayer)


        if has_right:
            mat_Tright, mat_expTR = self._bcs[1].layer.find_layer_transfer_matrix(Om, vP)
            # we just put the right hand matrix onto the end of the list. Doesn't need to be handled specially
            v_mat_Tlayer.append(mat_Tright)
        #mat_Glayers_done = np.eye(6, dtype=np.complex128) # hold the full layers we've already been through

        y_layer_start = 0
        ilayer_last = 0
        v_psi_layer_0 = v_psi0.copy() # initial field at start of current layer

        for iy,y in enumerate(v_y):

            ilayer = v_ilayer[iy]
            #print('lays', iy, ilayer, ilayer_last)
            if ilayer>0 and ilayer != ilayer_last: #entered a new layer, don't do it for Left layer to zero layer transfer
                L_last_layer = y-y_layer_start
                G_last_layer = spla.expm(v_mat_Tlayer[ilayer_last] * L_last_layer)
                v_psi_layer_0 = G_last_layer @ v_psi_layer_0 # push forward the start of the layer by one layer
                #mat_Glayers_done = G_last_layer @ mat_Glayers_done
                y_layer_start = y
                ilayer_last = ilayer



            if ilayer == -1:
                m_psi_y[:,iy] = spla.expm(mat_Tleft * y) @ v_psi0
            else:
                Dy = y - y_layer_start  # distance into current layer
                mat_T = v_mat_Tlayer[ilayer] #if ilayer >0 else mat_Tright
                m_psi_y[:,iy] = spla.expm(mat_T * Dy) @ v_psi_layer_0

        #with np.printoptions(precision=3, suppress=True, linewidth=200,
        #                     formatter={'float': '{: 0.3f}'.format}):
            #print(v_y, m_psi_y)
            #print(v_y, np.abs(m_psi_y[5,:]))
            #print(v_y[-5:],'\n', m_psi_y[:,-5:])

        #TODO: should probably unscale the position v_y
        return v_y, m_psi_y





class PiezoStack(ElasticStack):
    def __init__(self, l_layers: List[ElasticLayer], bcs: List[ElasticBoundaryCondition]) -> None:
        super().__init__(l_layers, bcs)











def plot_det_scan(Omega: float, v_Vp: NDArray[np.float64], v_det: NDArray[np.complex128], prefix: str = "tmp") -> None:

    fig, ax = plt.subplots()

    ax.plot(
        v_Vp, np.abs(v_det), "-+", label=r"$|\Delta|$", markersize=1.0, linewidth=0.5
    )

    ax.set_yscale("log")
    ax.set_xlabel(r"$V$ [km/s]", fontsize=10)
    ax.set_ylabel(r"$|T|$ ", fontsize=10)
    ax.set_title(rf"$\Omega/(2\pi)={Omega/(2*np.pi):.3f}\,$GHz", fontsize=10)
    ax.tick_params(axis='both', labelsize=10)
    ax.legend(fontsize=10)

    plpl.save_and_close_figure(fig, prefix + "_det.png")


def plot_rayleigh_det_scan(v_Vp: NDArray[np.float64], v_det: NDArray[np.complex128], Vs: float, Vl: float,
                           material#: materials.Material
                           ) -> None:

    fig, ax = plt.subplots()
    ax.plot(v_Vp * 1e-3, np.real(v_det), "o", label=r"$\Delta_r$", markersize=0.5)
    ax.plot(v_Vp * 1e-3, np.imag(v_det), "x", label=r"$\Delta_i$", markersize=0.5)
    ax.plot(v_Vp * 1e-3, np.abs(v_det), "+", label=r"$|\Delta|$", markersize=0.5)
    ax.set_xlabel(r"$V$ [km/s]", fontsize=10)
    ax.set_ylabel(r"$|T|$ ", fontsize=10)
    ax.tick_params(axis='both', labelsize=10)
    # ax.axvline(Vs*.001, ls=':', c='gray')
    # ax.axvline(Vacs[0], ls=':', c='gray')
    # ax.axvline(Vacs[1], ls=':', c='gray')
    # ax.axvline(Vacs[2], ls=':', c='gray')
    ax.legend(fontsize=10)

    vR = 0
    if material.is_isotropic():
        vR = materials.find_Rayleigh_velocity_isotropic(Vs, Vl)
        # ax.axvline(vR*.001, ls=':', c='gray')

    plpl.save_and_close_figure(fig, "ray_det.png")
    # print('Shear vel', Vs, vR)


def find_Rayleigh_velocity_anisotropic(material#: materials.Material
                                       ) -> float:
    """Find Rayleigh velocity for arbitrary material along the zhat=(0,0,1) direction with surface normal along yhat=(0,1,0)."""
    # use transfer matrix approach and search for vanishing determinant
    # Solution(s) should lie below the bulk shear velocity

    bcs = (
        ElasticBoundaryCondition(ElBCType.SEMI_INFINITE, material),
        ElasticBoundaryCondition(ElBCType.VACUUM),
    )
    l_layers = []


    Omega = 2*np.pi * 1.0 * SI_GHz
    estack = NonPiezoStack(l_layers, bcs)
    Vs = estack.dispersion_relation_find_Vs_at_Omega(Omega,
                                                         Vmin=1*SI_kmps, Vmax=8.0*SI_kmps)


    if Vs:
        return Vs[0][0]  # return first found mode
    else:
        return 0
