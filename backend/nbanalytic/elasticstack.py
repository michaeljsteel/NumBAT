import numpy as np
import numpy.linalg as npla
import scipy.linalg as spla
from typing import List, Tuple, Optional
from numpy.typing import NDArray

from enum import Enum

import matplotlib.pyplot as plt

import numbat
import nbtypes
import materials
from numtools.optimize import bracket_minima, find_bracketed_minima, BracketError
from nbtypes import SI_GHz, SI_um, SI_kmps
import plotting.plottools as nbpt
from nbanalytic.elasticfieldconversions import uvec3_to_poynting_S3

from nbanalytic.elasticmodeplots import ModeFunction1D, QuantLabel

g_norms = nbtypes.NormalisationConstants()


class ElBCType(Enum):
    VACUUM = "vacuum"
    SEMI_INFINITE = "semiinfinite"
    SHORT = "short"
    CHARGE_FREE = "chargefree"


def check_magnitude(x: float, name: str, minval: float, maxval: float) -> None:
    if not (minval <= abs(x) <= maxval):
        raise ValueError(f"{name} value {x} out of range [{minval}, {maxval}]")

# Coefficient matrices defined in the layers3.nb mathematica file
def get_layer_stiffness_submatrices(cs: NDArray[np.float64], use_4d: bool = False) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    # cs is normalised stiffness

    assert cs.shape == (7, 7)

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
        # modify for 8-component case
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

    assert cs.shape == (7, 7)
    assert e_iJ.shape == (4, 7)
    assert perm_ij.shape == (4, 4)


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
                 thickness_SI: float = 1e-6, use_4d: bool = False) -> None:
        # dimensionless thickness

        self.material = material
        self.L = thickness_SI / g_norms.x0

        self._build_material_quantities(use_4d)

    def thickness_SI(self) -> float:
        return self.L * g_norms.x0

    def _build_material_quantities(self, use_4d: bool) -> None:

        self.cstiff_norm = (
            self.material.stiffness_c_IJ.unit_indexed_value() / g_norms.c0
        )
        self.rho_norm = self.material.rho / g_norms.rho0

        if self.material.piezo_active():
            eiJ_norm = (
                self.material._piezo._tens_e_iJ.unit_indexed_value() / g_norms.e0
            )
            permij_norm = (
                self.material._piezo._tens_releps_S_ij.unit_indexed_value()
                / g_norms.eps0
            )
        else: # trivial values for non-piezo materials in a piezo-aware stack
            eiJ_norm = np.zeros((4, 7), dtype=np.float64)
            permij_norm = np.eye(4, dtype=np.float64)

        if self.material.piezo_active():
            mM1, mM2, mL1, mL2 = get_layer_stiffness_submatrices_piezo(
                self.cstiff_norm, eiJ_norm, permij_norm
            )
        else:
            mM1, mM2, mL1, mL2 = get_layer_stiffness_submatrices(
                self.cstiff_norm, use_4d
            )
        self.mM1, self.mM2, self.mL1, self.mL2 = mM1, mM2, mL1, mL2


    def find_layer_transfer_matrix(self, Om: float, Vp: float) -> Tuple[NDArray[np.complex128], Optional[NDArray[np.complex128]]]:
        # Om, Vp are both normalised
        check_magnitude(Om, "Om", 0.05, 1000)
        check_magnitude(Vp, "Vp", 0.01, 50)

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

        self.layer = None
        if material is not None:
            self.layer = ElasticLayer(material, 0.0)

        if self.is_semi_infinite() and material is None:
            raise ValueError("Semi-infinite boundary condition requires a material.")

        self._semiinf_evals = None
        self._semiinf_evecs = None

        self._is_front_side = True # is bc at front (left side) of stack


    def is_vacuum(self) -> bool:
        return self.bc_type.value == ElBCType.VACUUM.value  # using .value since we seem to be getting multiple visibly identical but not identical objects

    def is_semi_infinite(self) -> bool:
        return self.bc_type.value == ElBCType.SEMI_INFINITE.value

    def is_short(self) -> bool:
        return self.bc_type.value == ElBCType.SHORT.value

    def is_charge_free(self) -> bool:
        return self.bc_type.value == ElBCType.CHARGE_FREE.value

    def set_is_front(self):
        self._is_front_side = True

    def set_is_back(self):
        self._is_front_side = False

    def __str__(self) -> str:
        return f"ElasticBoundaryCondition({self.bc_type})"

    def analyse_semiinfinite_eigenspace(self, Om: float, Vp: float) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        # find eigenvalues and eigenvectors of semi-infinite layer transfer matrix
        # identify growing and "forwards" modes and put them at front of the list

       # Om, Vp are both normalised
        check_magnitude(Om, "Om", .05, 1000)
        check_magnitude(Vp, "Vp", 0.01, 50)

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
        v_poynting_Sav = np.zeros((len(eigvals), 3), dtype=np.float64)

        growing_is_good = self._is_front_side

        qz = Om/Vp

        for iev, ev in enumerate(eigvals): # ev is lambda = i q_y

            qy = ev / 1j
            vecq = np.array([0, qy, qz], dtype=np.complex128)
            vecu = eigvecs[3:6, iev]
            Svec = uvec3_to_poynting_S3(Om, vecu, vecq, self.layer.cstiff_norm[1:,1:]) # cstiff is unit-indexed


            v_is_growing[iev] = is_growing_eigenvalue(ev)
            v_is_decaying[iev] = is_decaying_eigenvalue(ev)
            v_is_fwds[iev] = np.isclose(np.real(ev), 0) and np.imag(ev) > 0
            v_is_pureimag[iev] = np.isclose(np.real(ev), 0)
            v_poynting_Sav[iev, :] = Svec


        ngrow = np.sum(v_is_growing)
        ndecay = np.sum(v_is_decaying)
        nfwd = np.sum(v_is_fwds)

        #take all the growing/decaying, and then as many outwards S modes as we need to fill it up
        if growing_is_good: # front bc
            sign_Sy = -1
            ev_good.extend(np.nonzero(v_is_growing)[0]) # take all the growing
        else: # back bc
            sign_Sy = 1
            ev_good.extend(np.nonzero(v_is_decaying)[0])

        # fill up the rest of the good slots with outward (leftward) going propagating waves
        for i in range(6):
            if v_is_pureimag[i] and v_poynting_Sav[i,1] * sign_Sy > 0: #  Sy is outwards
                ev_good.append(i)
            if len(ev_good)==3: break


        if len(ev_good) != 3:
            print(f'Unexpected number of good eigenvalues {len(ev_good)} in semi-infinite layer:')
            print(f' Om:{Om:.4f}, Vp:{Vp:.4f}, (ngrow={ngrow}, ndecay={ndecay}, nfwd={nfwd}, ngood={len(ev_good)}).')
            for (vprop, nm) in (
                (v_is_growing, "growing"),
                (v_is_decaying, "decaying"),
                (v_is_pureimag, "pureimag"),):

                for iev in vprop.nonzero()[0]:
                    Sx,Sy,Sz = v_poynting_Sav[iev,:] * 100
                    print(f'  ev[{iev}] = {eigvals[iev]: 12.5f} {nm},  Svec=({Sx:7.3f},{Sy:7.3f},{Sz:7.3f})')


        # order the eigensolutions

        ev_bad = list(set(range(6)) - set(ev_good)) # the bads are the values {0,1,2,3,4,5} minus the known good indices
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

        vEvals, mVeigs = self.analyse_semiinfinite_eigenspace(Om, Vp)
        # find smallest decay length
        lamrs = np.real(vEvals[0:3])
        lammin =np.min(lamrs)

        if not np.isclose(lammin,0):
            L= 6.0 / abs(lammin)  # 3 decay lengths
        else:
            L= 5.0 # arbitrary

        if L>100: # something has gone wrong, probably a nondecaying eigenvalue
            L=5.0
        return L


    def suggested_L_right(self, Om: float, Vp: float) -> float:
        return self.suggested_L_left(Om, Vp)


    def null_vector_to_field_psi0(self, v_null):
        hdim = len(v_null)

        if self.is_vacuum():
            v_psi0 = np.zeros(2*hdim, dtype=np.complex128)
            v_psi0[hdim:] = v_null  # null is literally initial displacement field
        elif self.is_semi_infinite():
            v_psi0 = self._semiinf_evecs[:,0:hdim] @ v_null  # initial field from evectors
        else:
            raise NotImplementedError(
                "Boundary condition type not implemented for null vector to field conversion.")

        return v_psi0


class ElasticStack:
    def __init__(
        self, layers: List[ElasticLayer], bcs: List[ElasticBoundaryCondition]
    ):
        self._def_Vmin = 0.01  # km/s
        self._def_Vmax = 10.0  # km/s
        self._n_layers = len(layers)

        self._layers = layers
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

        self._bcs[0].set_is_front()
        self._bcs[1].set_is_back()

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

    # def make_whole_stack_Pmat(self):
    #     raise NotImplementedError(
    #         "Base class method called - should use subclass method."
    #     )

    # def dispersion_relation_find_Vs_at_Omegas(self, v_Omega, Vmin=0, Vmax=0):
    #     raise NotImplementedError(
    #         "Base class method called - should use subclass method."
    #     )

    def dispersion_relation_find_Vs_at_Omega(self, Omega, Vmin=0, Vmax=0):
        raise NotImplementedError(
            "Base class method called - should use subclass method."
        )

    def build_profile_domain_y(self, ny_per_layer: int = 50,
                               L_left: float = 0, L_right: float = 0) -> NDArray[np.float64]:
        """Make vectors v_y and v_ilayer for the y-domain of the mode profile."""

        # get semi-inf borders
        if self._bcs[0].is_semi_infinite(): # add a bit of space for decay
            L_left = L_left if L_left > 0 else 5.0
        else:
            L_left = 0.0
        if self._bcs[1].is_semi_infinite(): # add a bit of space for decay
            L_right = L_right if L_right > 0 else 5.0
        else:
            L_right = 0.0


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

        layer_lengths = [layer.L for layer in self._layers]
        if L_right>0:
            layer_lengths.append(L_right)

        for ilay,L in enumerate(layer_lengths):
            yhi += L
            dy = L / ny_per_layer
            v_y[iyoff : iyoff + ny_per_layer] = np.linspace(ylo, yhi - dy, ny_per_layer)
            v_ilayer[iyoff : iyoff + ny_per_layer] = ilay
            iyoff += ny_per_layer
            ylo += L

        # if L_right > 0:
        #     yhi += L_right
        #     dy = L_right / ny_per_layer
        #     v_y[iyoff : iyoff + ny_per_layer] = np.linspace(ylo, yhi - dy, ny_per_layer)
        #     v_ilayer[iyoff : iyoff + ny_per_layer] = -2
        #     iyoff += ny_per_layer
        #     ylo += L_right

        v_y[-1] = v_y[-2] + dy
        v_ilayer[-1] = v_ilayer[-2]

        if L_right>0: # overwrite code for
            v_ilayer[-(ny_per_layer+1):] = -2


        return v_y, v_ilayer

    def find_mode_profile_1d(self, Omega: float, Vs: float,
                             ny_per_layer: int = 50, L_left: float = -1, L_right: float = -1) -> NDArray[np.float64]:
        pass






class NonPiezoStack(ElasticStack):
    def __init__(self, l_layers: List[ElasticLayer], bcs: List[ElasticBoundaryCondition]) -> None:
        super().__init__(l_layers, bcs)

    def _find_full_transfer_matrix(self, Om: float, Vp: float) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        # build layer parts
        mat_G = np.eye(6, dtype=np.complex128)

        for layer in self._layers:
            mat_T, mat_expTL = layer.find_layer_transfer_matrix(Om, Vp)
            mat_G = mat_expTL @ mat_G


        # add boundary layer parts
        for ibcs, bc in enumerate(self._bcs):

            if bc.is_semi_infinite():

                # sort eigenvectors into growing and decaying/forwards and backwards modes
                v_evals, m_eigvecs = bc.analyse_semiinfinite_eigenspace(Om, Vp)

                if ibcs._is_front_side:
                    mat_G = mat_G @ m_eigvecs
                else:
                    mat_G = npla.inv(m_eigvecs) @ mat_G

        mat_G_ij = self._extract_quadrant(mat_G)
        mat_G_ij /= np.max(np.abs(mat_G_ij))  # scale to avoid overflow in det

        return mat_G, mat_G_ij

    def dispersion_relation_find_Vs_at_Omega(
        self,
        Omega_SI: float,
        Vmin_SI: float = 10,
        Vmax_SI: float = 10000.0,
        nVpts: int = 1000,
        find_mode_color: bool = False,
        prefix: str = "tmpdrscan",
        label: str = '',
        show_scan_plot: bool = False,
        show_progress: bool = False,
        fmin_thresh: float = 1e-5,

    ) -> List[Tuple[float, float]]:

        check_magnitude(Omega_SI, "Omega", 1e8,100e9*2*np.pi)
        check_magnitude(Vmin_SI, "Vmin", 10, 20000)
        check_magnitude(Vmax_SI, "Vmax", 100, 20000)


        # move to normalised units here
        Om = Omega_SI / g_norms.f0

        if label:
            pllabel=label
        else:
            label = f'Nu={Om/(2*np.pi):.10f}'
            pllabel = rf'\Omega/(2\pi)={Om/(2*np.pi)} [GHz]'


        n_Vp = nVpts

        v_Vp = np.linspace(Vmin_SI, Vmax_SI, n_Vp) / g_norms.v0
        # choose to space the velocities according to equal qz
        #v_Vp = 1.0/np.linspace(1/Vmin, 1/Vmax, n_Vp) / g_norms.v0


        v_det = np.zeros(n_Vp, dtype=np.complex128)

        for ivP, vP in enumerate(v_Vp):
            mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
            v_det[ivP] = npla.det(mat_G_ij)
            #print(f"detG Om {Om:.3f} Vp {vP:.3f}:", npla.det(mat_G_ij))
            #print('detmag', ivP, np.max(np.abs(mat_G_ij)))

        if show_scan_plot:
            plot_det_scan(Om, v_Vp, v_det, prefix=prefix, pllabel=pllabel)


        l_braks = bracket_minima(v_Vp, np.abs(v_det))
        if show_progress:
            print(f"\nAt {label}, found {len(l_braks)} brackets:")

            for br in l_braks:
                (xm1, fm1), (x0, f0), (xp1, fp1) = np.array(br).tolist() # awkward but converts to plain floats
                print(f"    [{xm1:.5f}, {x0:.5f}, {xp1:.5f}] with f=[{fm1:.3e}, {f0:.3e}, {fp1:.3e}]")

        if len(l_braks) == 0:
            plot_det_scan(Om, v_Vp, v_det, prefix=prefix, pllabel=pllabel)
            raise BracketError("No minima found in determinant scan - no modes found. Determinant scan plot dumped.")

        def minfunc(vP):
            mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
            return np.abs(npla.det(mat_G_ij))

        l_mins = find_bracketed_minima(minfunc, l_braks)
        if show_progress:
            print(f"  Found {len(l_mins)} minima:")
            for v, f in l_mins:
                print(f"    Vp={v:.9f}, Det(f)={f:.4e}")



        # remove any minima with f above threshold or with multiple singular values as probably not real modes
        fmin_norm_thresh = fmin_thresh * np.max(np.abs(v_det))
        l_mins_V_SI_keep = []
        l_mins_bad_threshold = []
        l_mins_bad_singval = []

        for v, f in l_mins:
            if f > fmin_norm_thresh:
                l_mins_bad_threshold.append((v,f))
            #elif self.count_singular_values(Om, v) != 1:
            #    l_mins_bad_singval.append((v,f))  # could be zero or two singvals
            else:
                l_mins_V_SI_keep.append((v * g_norms.v0,f))

        n_dropped = len(l_mins_bad_threshold) + len(l_mins_bad_singval)
        if show_progress and n_dropped > 0:
            print(f"  Dropping {n_dropped} modes:")
            for v, f in l_mins_bad_singval:
                print(f"    Vp={v:.3f} m/s has multiple small singular values.")
            for v, f in l_mins_bad_threshold:
                print(f"    Vp={v:.3f} has relative det minimum f={f:.3e} above fmin_thresh={fmin_norm_thresh:.3e}:")

        # color modes by polarization at the surface.
        # better to do by energy weighted average over all domain but that is more expensive
        v_mode_cols = None

        full_mode_col = False
        if find_mode_color:
            v_mode_cols = np.zeros((len(l_mins_V_SI_keep), 3), dtype=np.float64 )
            for iv, (vSI,f) in enumerate(l_mins_V_SI_keep):
                if full_mode_col:
                    pol_vec = self.find_mode_polarization_over_profile_1d(Omega_SI, vSI)
                else:
                    pol_vec = self.find_mode_polarization_at_y0(Omega_SI, vSI)

                v_mode_cols[iv,:] = nbpt.get_rgb_for_poln(*pol_vec)
        res = { 'Vphase': l_mins_V_SI_keep,
                'poln_cols': v_mode_cols,
                'det_scan': (v_Vp, v_det)
                }

        return res

    def find_mode_polarization_at_y0(self, Omega_SI: float, vP_SI: float) -> NDArray[np.float64]:
        v_psi0 = self.find_mode_psi0(Omega_SI/g_norms.f0, vP_SI/g_norms.v0)
        u_x, u_y, u_z = np.abs(v_psi0[3:6])
        poln_vec = np.array([u_x, u_y, u_z], dtype=np.float64)
        poln_vec /= np.linalg.norm(poln_vec)
        return poln_vec

    def find_mode_polarization_over_profile_1d(self, Omega_SI: float, vP_SI: float,
                             ny_per_layer: int = 100, L_left: float = -1, L_right: float = -1) -> NDArray[np.float64]:


        v_y, m_psi_y = self.find_mode_profile_1d(Omega_SI, vP_SI)
        m_u_y = m_psi_y[3:6, :]  # displacement part
        poln_fracs = np.sum(np.abs(m_u_y)**2, axis=1)
        poln_fracs /= np.linalg.norm(poln_fracs)
        return poln_fracs



    def count_singular_values(self, Om: float, vP: float,
                              thresh: float = 1e-6) -> int:
        mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
        u, s, vh = npla.svd(mat_G_ij)
        nsmall = np.sum(s < thresh*np.max(s)) # Limited contrast precision in SVD
        return nsmall

    def find_mode_psi0(self, Om: float, vP: float, thresh_singavl: float = 1e-6,
                       show_progress: bool = False) -> NDArray[np.complex128]:

        mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)

        u, s, vh = npla.svd(mat_G_ij)

        if show_progress:
            print('Det check:', npla.det(mat_G_ij))
            print('SVD check:', s)

        nsmall = np.sum(s < thresh_singavl*np.max(s)) # Limited contrast precision in SVD

        if show_progress:
            print('Num small singular values:', nsmall)

        if nsmall == 0:
            raise ValueError('No nullspace found - not a mode?')
        if nsmall > 1:
            print('Warning: Multiple small singular values found - may be non-genuine mode.')
            #raise ValueError('Multiple small singular values found - may be non-genuine mode.')


        for isv in range(len(s)):
            #if s[isv] < thresh_singavl*np.max(s):
            #    print(f'  Small singval {isv}: {s[isv]:.3e}')

            v_null = vh.conj().T[:,-1]  # take the last column of Vh (smallest singular value)

            v_psi0 = self._bcs[0].null_vector_to_field_psi0(v_null)
            v_psi0 /= np.max(np.abs(v_psi0[3:6]))  # normalise max displacement to 1

            if show_progress:
                with np.printoptions(precision=4, suppress=True, linewidth=200):
                    print(f'Initial field at y=0:\n null vec: {v_null}\n psi0:  {v_psi0}')

        return v_psi0

    def find_mode_profile_1d(self, Omega_SI: float, Vp_SI: float,
                             ny_per_layer: int = 100, L_left: float = -1, L_right: float = -1) -> NDArray[np.float64]:
        # find nullspace vector and initial field at y=0
        Om = Omega_SI / g_norms.f0
        vP = Vp_SI / g_norms.v0

        try:
            v_psi0 = self.find_mode_psi0(Om, vP)
        except ValueError as e:
            raise ValueError(f"Could not find mode at Omega={Om:.3e} ({Om:.3f} norm), Vs={vP:.3f} ({vP:.3f} norm): {e}") from e

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
        for ilay,layer in enumerate(self._layers):
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

            # We've entered a new layer, update the starting vector for the layer by propagating over the full previous layer.
            # We don't do it for the Left layer to zero layer transfer since nothing has happened
            if ilayer>0 and ilayer != ilayer_last:
                L_last_layer = y-y_layer_start
                G_last_layer = spla.expm(v_mat_Tlayer[ilayer_last] * L_last_layer)
                v_psi_layer_0 = G_last_layer @ v_psi_layer_0 # push forward the start of the layer by one layer
                #mat_Glayers_done = G_last_layer @ mat_Glayers_done
                y_layer_start = y
                ilayer_last = ilayer



            if ilayer == -1: # in the first layer. Back-propagate from the y=0 point (y is negative)
                m_psi_y[:,iy] = spla.expm(mat_Tleft * y) @ v_psi0
            else:
                Dy = y - y_layer_start  # distance into current layer
                mat_T = v_mat_Tlayer[ilayer] if ilayer >0 else v_mat_Tlayer[-1]  # if we are in right layer (ilayer=-2), then use the last layer!
                m_psi_y[:,iy] = spla.expm(mat_T * Dy) @ v_psi_layer_0

        #with np.printoptions(precision=3, suppress=True, linewidth=200,
        #                     formatter={'float': '{: 0.3f}'.format}):
            #print(v_y, m_psi_y)
            #print(v_y, np.abs(m_psi_y[5,:]))
            #print(v_y[-5:],'\n', m_psi_y[:,-5:])

        #TODO: should probably unscale the position v_y

        v_y_SI = v_y * g_norms.x0


        m_uxyz = m_psi_y[3:,:].T.copy()
        m_T6 = m_psi_y[0:3,:].T.copy()

        # to make physically realistic values  will treat displacements as being in nm, so stresses come down by 1000
        m_T6 /=1000

        mode_prof_u = ModeFunction1D(Omega_SI, Vp_SI, v_y_SI, m_uxyz)
        mode_prof_T6 = ModeFunction1D(Omega_SI, Vp_SI, v_y_SI, m_T6)

        mode_prof_T6.set_labels(dvar=QuantLabel('Normal stress', r'$\vec T_i$', r'$|\vec T_i|$', '[GPa]'))


        #return v_y, m_psi_y, mode_prof_u, mode_prof_T6

        return  mode_prof_u, mode_prof_T6




class PiezoStack(ElasticStack):
    def __init__(self, l_layers: List[ElasticLayer], bcs: List[ElasticBoundaryCondition]) -> None:
        super().__init__(l_layers, bcs)











def plot_det_scan(Omega: float, v_Vp: NDArray[np.float64], v_det: NDArray[np.complex128],
                  prefix: str = "tmp", pllabel: str = "") -> None:

    fig, ax = plt.subplots()

    ax.plot(
        v_Vp, np.abs(v_det), "-+", label=r"$|\Delta|$", markersize=1.0, linewidth=0.5
    )

    ax.set_yscale("log")
    ax.set_xlabel(r"$V$ [km/s]", fontsize=10)
    ax.set_ylabel(r"$|T|$ ", fontsize=10)
    #ax.set_title(rf"$\Omega/(2\pi)={Omega/(2*np.pi):.3f}\,$GHz", fontsize=10)
    ax.set_title(pllabel, fontsize=10)
    ax.tick_params(axis='both', labelsize=10)
    ax.legend(fontsize=10)

    nbpt.save_and_close_figure(fig, prefix + "_det.png")


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

    nbpt.save_and_close_figure(fig, "ray_det.png")
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
    res = estack.dispersion_relation_find_Vs_at_Omega(Omega, Vmin_SI=1*SI_kmps, Vmax_SI=8.0*SI_kmps)
    Vp = res['Vphase']

    if Vp:
        return Vp[0][0]  # return first found mode
    else:
        return 0


# def plot_profile_2d(v_y, m_psi_y, Omega, Vs, imd, prefix, Vbulks=None):
#     fig,ax=plt.subplots(3,2, figsize=(6,4))
#     fig.set_tight_layout(True)

#     nutil=Omega/(2*np.pi*SI_GHz)
#     Vstil=Vs/SI_kmps
#     title = fr'Mode {imd} for $\Omega/2\pi$={nutil:.3f} GHz, V={Vstil:.3f} km/s'
#     if Vbulks is not None:
#         title += f'\nBulk Velocities: {Vbulks[0]:.3f}, {Vbulks[1]:.3f}, {Vbulks[2]:.3f} km/s'
#     fig.suptitle(title, fontsize=8)

#     q= Omega / Vs * SI_um
#     v_z = np.linspace(0, 2*np.pi/q, 100)
#     v_Y, v_Z = np.meshgrid(v_y, v_z)

#     v_T6 = np.abs(m_psi_y[0,:])
#     v_T2 = np.abs(m_psi_y[1,:])
#     v_T4 = np.abs(m_psi_y[2,:])
#     v_ux = np.abs(m_psi_y[3,:])
#     v_uy = np.abs(m_psi_y[4,:])
#     v_uz = np.abs(m_psi_y[5,:])
#     m_eiphiz = np.exp(1j * q * v_Z).T
#     m_T6 = np.real(np.outer(v_T6, np.ones(len(v_z))) * m_eiphiz)
#     m_T2 = np.real(np.outer(v_T2, np.ones(len(v_z))) * m_eiphiz)
#     m_T4 = np.real(np.outer(v_T4, np.ones(len(v_z))) * m_eiphiz)
#     m_ux = np.real(np.outer(v_ux, np.ones(len(v_z))) * m_eiphiz)
#     m_uy = np.real(np.outer(v_uy, np.ones(len(v_z))) * m_eiphiz)
#     m_uz = np.real(np.outer(v_uz, np.ones(len(v_z))) * m_eiphiz)

#     exts = (v_y[0], v_y[-1], v_z[0], v_z[-1])
#     fs = 6
#     plotprefs = numbat.NumBATPlotPrefs()
#     cm = plotprefs.cmap_ac_field_signed
#     kwimgs = dict(extent=exts, aspect='auto', origin='lower', cmap=cm)

#     for iax, (m_T,mlab) in enumerate(
#         ([m_T6, r'$\Re\{T_6\}$ [GPa]'],
#          [m_T2, r'$\Re\{T_2\}$ [GPa]'],
#          [m_T4, r'$\Re\{T_4\}$ [GPa]'],
#          [m_ux, r'$\Re\{u_x\}$ [$\mu$m]'],
#          [m_uy, r'$\Re\{u_y\}$ [$\mu$m]'],
#          [m_uz, r'$\Re\{u_z\}$ [$\mu$m]'])):
#         irow = iax % 3
#         icol = iax // 3
#         im = ax[irow, icol].imshow(m_T.T, **kwimgs)
#         cbar = fig.colorbar(im, ax=ax[irow, icol])
#         cbar.ax.tick_params(labelsize=fs)
#         ax[irow, icol].set_title(mlab, fontsize=fs)


#     for tax in ax.flatten():
#         tax.tick_params(axis='both', labelsize=fs)
#         tax.set_xlabel(r'$y$ [$\mu$m]', fontsize=fs)
#         tax.set_ylabel(r'$z$ [$\mu$m]', fontsize=fs)

#     fname = prefix+'_profile_2d.png'
#     nbpt.save_and_close_figure(fig, fname)
#     return fname


# def plot_profile_1d(v_y, m_psi_y, Omega, Vs, imd, prefix, Vbulks=None):
#     fig,ax=plt.subplots(2,2, figsize=(6,4))
#     fig.set_tight_layout(True)

#     nutil=Omega/(2*np.pi*SI_GHz)
#     Vstil=Vs/SI_kmps
#     title = fr'Mode {imd} for $\Omega/2\pi$={nutil:.3f} GHz, V={Vstil:.3f} km/s'
#     if Vbulks is not None:
#         title += f'\nBulk Velocities: {Vbulks[0]:.3f}, {Vbulks[1]:.3f}, {Vbulks[2]:.3f} km/s'
#     fig.suptitle(title, fontsize=8)

#     ms = 1
#     v_T6 = m_psi_y[0,:]
#     v_T2 = m_psi_y[1,:]
#     v_T4 = m_psi_y[2,:]
#     v_ux = m_psi_y[3,:]
#     v_uy = m_psi_y[4,:]
#     v_uz = m_psi_y[5,:]

#     lw=.5
#     lsr='-'
#     lsi='--'
#     marker='o'
#     c1='red'
#     c2='green'
#     c3='blue'
#     fs=6
#     kwlns = dict(lw=lw, marker=marker, markersize=ms)

#     ax[0,0].plot(v_y, np.real(v_T6), linestyle=lsr, **kwlns, color=c1, label=r'$T_6$')
#     ax[0,0].plot(v_y, np.real(v_T2), linestyle=lsr, **kwlns, color=c2, label=r'$T_2$')
#     ax[0,0].plot(v_y, np.real(v_T4), linestyle=lsr, **kwlns, color=c3, label=r'$T_4$')
#     ax[0,0].plot(v_y, np.imag(v_T6), linestyle=lsi, **kwlns, color=c1)
#     ax[0,0].plot(v_y, np.imag(v_T2), linestyle=lsi, **kwlns, color=c2)
#     ax[0,0].plot(v_y, np.imag(v_T4), linestyle=lsi, **kwlns, color=c3)

#     ax[1,0].plot(v_y, np.abs(v_T6), linestyle=lsr, **kwlns, color=c1, label=r'$T_6$')
#     ax[1,0].plot(v_y, np.abs(v_T2), linestyle=lsr, **kwlns, color=c2, label=r'$T_2$')
#     ax[1,0].plot(v_y, np.abs(v_T4), linestyle=lsr, **kwlns, color=c3, label=r'$T_4$')


#     ax[0,1].plot(v_y, np.real(v_ux), linestyle=lsr, **kwlns, color=c1, label=r'$u_x$')
#     ax[0,1].plot(v_y, np.real(v_uy), linestyle=lsr, **kwlns, color=c2, label=r'$u_y$')
#     ax[0,1].plot(v_y, np.real(v_uz), linestyle=lsr, **kwlns, color=c3, label=r'$u_z$')
#     ax[0,1].plot(v_y, np.imag(v_ux), linestyle=lsi, **kwlns, color=c1)
#     ax[0,1].plot(v_y, np.imag(v_uy), linestyle=lsi, **kwlns, color=c2)
#     ax[0,1].plot(v_y, np.imag(v_uz), linestyle=lsi, **kwlns, color=c3)


#     ax[1,1].plot(v_y, np.abs(v_ux), linestyle=lsr, **kwlns, color=c1, label=r'$u_x$')
#     ax[1,1].plot(v_y, np.abs(v_uy), linestyle=lsr, **kwlns, color=c2, label=r'$u_y$')
#     ax[1,1].plot(v_y, np.abs(v_uz), linestyle=lsr, **kwlns, color=c3, label=r'$u_z$')


#     for tax in ax.flatten():
#         tax.set_xlabel(r'Position [$\mu$m]', fontsize=fs)
#         tax.legend(fontsize=fs)
#         tax.tick_params(axis='both', labelsize=fs)

#     ax[0,0].set_ylabel(r'Normal stress $\vec T_{r,i}$ [GPa]', fontsize=fs)
#     ax[1,0].set_ylabel(r'Normal stress $|\vec T|$ [GPa]', fontsize=fs)

#     ax[0,1].set_ylabel(r'Displacement $\vec u_{r,i}$ [$\mu$m]', fontsize=fs)
#     ax[1,1].set_ylabel(r'Displacement $|\vec u|$ [$\mu$m]', fontsize=fs)

#     fname = prefix+'_profile_1d.png'
#     nbpt.save_and_close_figure(fig, fname)
#     return fname


