import numpy as np
import numpy.linalg as npla
import scipy.linalg as spla
from typing import List, Tuple, Optional, Dict, Any, Sequence
from numpy.typing import NDArray

from enum import Enum

import matplotlib.pyplot as plt

import nbtypes
from numtools.optimize import bracket_minima, find_bracketed_minima, BracketError
import numbattools as nbt
from nbtypes import SI_GHz, SI_kmps
import plotting.plottools as nbpt
from nbanalytic.elasticfieldconversions import uvec3_to_poynting_S3, uvec3_phi_to_poynting_S3_piezo

from nbanalytic.elasticmodeplots import ModeFunction1D, QuantLabel

g_norms = nbtypes.NormalisationConstants()


"""Solves elastic wave propagation in layered media using the transfer matrix method.
   Assumes propagation along z with layers along y.
"""


class ElBCType(Enum):
    """Boundary condition types for the elastic stack edges.

    Members
    -------
    VACUUM
        Free surface (zero tractions at the boundary).
    SEMI_INFINITE
        Semi-infinite half-space made of a material, used to absorb/emit waves.
    SHORT
        Electrically shorted boundary for piezoelectric problems (swaps rows to
        enforce phi=0 at the surface).
    CHARGE_FREE
        Charge-free boundary for piezoelectric problems (tau_phi coupling used).
    """
    VACUUM = "vacuum"
    SEMI_INFINITE = "semiinfinite"
    SHORT = "short"
    CHARGE_FREE = "chargefree"







class ElasticLayer:
    """Single homogeneous layer in the elastic stack.

    Parameters
    ----------
    material
        Material object providing stiffness, density, and (optionally) piezo tensors.
    thickness_SI
        Physical thickness in metres.
    use_4d
        Force an 8-component state dimension (piezo-aware) even for non-piezo materials.

    Notes
    -----
    All internal calculations are performed in normalised units defined by
    ``nbtypes.NormalisationConstants``. Access physical thickness via
    ``thickness_SI()``.
    """
    def __init__(
        self,
        material,  #: materials.Material,
        thickness_SI: float = 0.0,
        use_4d: bool = False,
    ) -> None:

        nbt.check_magnitude(thickness_SI, "thickness_SI", 0.0, 1e-4) # max length 100 microns

        self.material = material
        self.L = thickness_SI / g_norms.x0

        self.dim = 8 if use_4d or material.piezo_active() else 6
        self.hdim = self.dim // 2

        self.cstiff_norm: NDArray[np.float64] = None
        self.rho_norm: NDArray[np.float64] = None
        self.eiJ_norm: NDArray[np.float64] = None
        self.permij_norm: NDArray[np.float64] = None

        self._build_material_quantities(use_4d)

    def thickness_SI(self) -> float:
        """Return the physical thickness of the layer in metres."""
        return self.L * g_norms.x0

    def is_piezo(self) -> bool:
        """Whether the layer is piezoelectrically active."""
        return self.material.piezo_active()

    # TODO: this belongs somewhere else but involves a lot of arg passing.
    def poynting_vector_Sav(self, Om: float, Vp: float, vecq: NDArray[np.complex128],
                            evec: NDArray[np.complex128]) -> NDArray[np.float64]:
        """Compute time-averaged elastic Poynting vector S for a mode.

        Parameters
        ----------
        Om, Vp
            Normalised angular frequency and phase velocity.
        vecq
            Wavevector components as a 3-vector [qx, qy, qz] in normalised units.
        evec
            Eigenvector (state vector) corresponding to the local eigenvalue.

        Returns
        -------
        ndarray
            S = [Sx, Sy, Sz] in normalised power-flow units.
        """

        if self.dim == 6: # plain case
            vecu = evec[self.hdim:]
            Svec = uvec3_to_poynting_S3(Om, vecu, vecq, self.cstiff_norm[1:, 1:])
        else:  # piezo case
            vecu = evec[self.hdim:self.dim-1]
            phi = evec[-1]
            Svec = uvec3_phi_to_poynting_S3_piezo(Om, vecu, vecq, phi,
                                        self.cstiff_norm[1:, 1:],
                                        self.eiJ_norm[1:, 1:],
                                        self.permij_norm[1:, 1:])
        return Svec


    def find_layer_transfer_matrix(
        self, Om: float, Vp: float
    ) -> Tuple[NDArray[np.complex128], Optional[NDArray[np.complex128]]]:
        """Construct the layer transfer generator and its exponential.

        Parameters
        ----------
        Om, Vp
            Normalised angular frequency and phase velocity.

        Returns
        -------
        matT : ndarray (dim x dim)
            First-order system matrix d/dy psi = matT psi.
        matPdz : ndarray or None
            exp(matT * L) if the layer thickness L>0, else None for zero-thickness
            layers used to represent boundary conditions.
        """
        # Om, Vp are both normalised
        #nbt.check_magnitude(Om, "Om", 0.05, 1000)
        #nbt.check_magnitude(Vp, "Vp", 0.01, 50)

        dim, hdim = self.dim, self.hdim

        matT = np.zeros([dim, dim], dtype=np.complex128)

        qz = Om / Vp
        mjqz = -1j * qz
        # ddz [T] = [ M_TT  M_TU] [T]
        #     [U] = [ M_UT  M_UU] [U]

        matT[:hdim, :hdim] = mjqz * self.mL2_mM1inv  # M_TT

        matT[:hdim, hdim:dim] = (qz) ** 2 * (
            self.mL1 - self.rho_norm_I * Vp**2 - self.mL2_mM1inv_mM2
        )  # M_TU

        matT[hdim:dim, 0:hdim] = self.mM1inv  # M_UT

        matT[hdim:dim, hdim:dim] = mjqz * self.mM1inv_mM2  # M_UU

        # exponentiate to get layer transfer matrix
        matPdz = spla.expm(matT * self.L) if self.L > 0 else None

        return matT, matPdz # type: ignore[arg-type]  # pylance worries that matPdz might be a real matrix


    def _build_material_quantities(self, use_4d: bool) -> None:
        """Populate normalised material tensors used by the transfer matrix.

        Sets stiffness, density, and (if present) piezo coupling and permittivity,
        and precomputes helper matrices used in the layer transfer operator.
        """

        self.cstiff_norm = (
            self.material.stiffness_c_IJ.unit_indexed_value() / g_norms.c0
        )
        self.rho_norm = self.material.rho / g_norms.rho0

        if self.material.piezo_active():
            self.eiJ_norm = self.material._piezo._tens_e_iJ.unit_indexed_value() / g_norms.e0
            self.permij_norm = (
                self.material._piezo._tens_relepsS_ij.unit_indexed_value()
            )
        else:  # trivial values for non-piezo materials in a piezo-aware stack
            self.eiJ_norm = np.zeros((4, 7), dtype=np.float64)
            self.permij_norm = np.eye(4, dtype=np.float64)

        if self.material.piezo_active():
            mM1, mM2, mL1, mL2 = _get_layer_stiffness_submatrices_piezo(
                self.cstiff_norm, self.eiJ_norm, self.permij_norm
            )
        else:
            mM1, mM2, mL1, mL2 = _get_layer_stiffness_submatrices(
                self.cstiff_norm, use_4d
            )

        self.rho_norm_I = self.rho_norm * np.eye(self.hdim)
        self.mL1 = mL1
        self.mM1inv = npla.inv(mM1)
        self.mL2_mM1inv = mL2 @ self.mM1inv
        self.mL2_mM1inv_mM2 = mL2 @ self.mM1inv @ mM2
        self.mM1inv_mM2 = self.mM1inv @ mM2


class ElasticBoundaryCondition:
    """Boundary condition wrapper, optionally with a semi-infinite layer.

    Depending on ``bc_type``, this may host an internal zero-thickness
    ``ElasticLayer`` (for SEMI_INFINITE), or just encode algebraic constraints
    (VACUUM/SHORT/CHARGE_FREE).
    """
    def __init__(
        self,
        bc_type: ElBCType,
        material=None,
        # material: Optional[materials.Material] = None
    ) -> None:
        self.bc_type = bc_type
        self.material = material

        self._layer = None
        if material is not None:
            self._layer = ElasticLayer(material)

        if self.is_semi_infinite() and material is None:
            raise ValueError("Semi-infinite boundary condition requires a material.")

        self._semiinf_evals = None
        self._semiinf_evecs = None

        self._is_front_side = True  # is bc at front (left side) of stack



    def is_piezo(self) -> bool:
        """Whether this boundary condition engages piezo behaviour."""
        return self._layer is not None and self._layer.is_piezo()

    def is_vacuum(self) -> bool:
        """True for a free-surface (vacuum) boundary."""
        # using .value since we seem to be getting multiple visibly identical but not identical objects
        return self.bc_type.value == ElBCType.VACUUM.value

    def is_semi_infinite(self) -> bool:
        """True if this side is a semi-infinite half-space of a material."""
        return self.bc_type.value == ElBCType.SEMI_INFINITE.value

    def is_short(self) -> bool:
        """True for electrically shorted boundary (piezo)."""
        return self.bc_type.value == ElBCType.SHORT.value

    def is_charge_free(self) -> bool:
        """True for charge-free boundary (piezo)."""
        return self.bc_type.value == ElBCType.CHARGE_FREE.value

    def set_is_front(self) -> None:
        """Mark this BC as the front (left) side of the stack."""
        self._is_front_side = True

    def set_is_back(self) -> None:
        """Mark this BC as the back (right) side of the stack."""
        self._is_front_side = False

    def __str__(self) -> str:
        """Human-readable description of the boundary condition."""
        return f"ElasticBoundaryCondition({self.bc_type})"

    def find_layer_transfer_matrix(
        self, Om: float, Vp: float
    ) -> Tuple[NDArray[np.complex128], Optional[NDArray[np.complex128]]]:
        """Get the transfer matrix of the semi-infinite layer (if any).

        Parameters
        ----------
        Om, Vp
            Normalised angular frequency and phase velocity.

        Returns
        -------
        matT : ndarray (dim x dim)
            First-order system matrix d/dy psi = matT psi.
        matPdz : ndarray or None
            exp(matT * L) if the layer thickness L>0, else None for zero-thickness
            layers used to represent boundary conditions.
        """

        assert self._layer is not None, "Layer must be initialized for semi-infinite BC" # helps with typing

        return self._layer.find_layer_transfer_matrix(Om, Vp)

    def analyse_semiinfinite_eigenspace(
        self, Om: float, Vp: float
    ) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        """Analyse the semi-infinite medium and sort eigenmodes.

        Finds eigenvalues/vectors of the transfer operator for the half-space and
        orders the modes so that the first half correspond to physically
        admissible solutions (growing away from the interface on the front side
        or decaying on the back side, plus outward-propagating real modes).

        Parameters
        ----------
        Om, Vp
            Normalised angular frequency and phase velocity.

        Returns
        -------
        eigvals, eigvecs : ndarray
            Eigenvalues (ordered) and the corresponding eigenvectors as columns.
        """



        # Om, Vp are both normalised
        nbt.check_magnitude(Om, "Om", 0.05, 1000)
        nbt.check_magnitude(Vp, "Vp", 0.01, 50)

        mat_T, mat_expTL = self.find_layer_transfer_matrix(Om, Vp)

        evals: NDArray[np.complex128]
        evecs: NDArray[np.complex128]
        evals, evecs = spla.eig(mat_T)  # type: ignore[assignment]

        evals_sorted, evecs_sorted = self._classify_seminifite_eigenspace(
            Om, Vp, evals, evecs)

        self._semiinf_evals = evals_sorted
        self._semiinf_evecs = evecs_sorted

        return evals_sorted, evecs_sorted

    def _classify_seminifite_eigenspace(self, Om:float, Vp: float,
                                    eigvals: NDArray[np.complex128],
                                    eigvecs: NDArray[np.complex128]):

        assert self._layer is not None, "Layer must be initialized for semi-infinite BC" # helps with typing

        dim = len(eigvals)
        hdim = dim//2

        ev_good = []
        ev_bad = []
        ngrow = 0
        ndecay = 0
        nfwd = 0

        v_is_growing = np.zeros(dim, dtype=bool)
        v_is_decaying = np.zeros(dim, dtype=bool)
        v_is_fwds = np.zeros(dim, dtype=bool)
        v_is_pureimag = np.zeros(dim, dtype=bool)
        v_poynting_Sav = np.zeros((dim, 3), dtype=np.float64)



        tol = 1e-11
        qz = Om / Vp

        # this loop is time critical. don't factor into smaller functions or use np.isclose.
        # everything starts false, so only need to update if possibly true
        for iev, ev in enumerate(eigvals):  # ev is lambda = i q_y

            qy = ev / 1j
            evec = eigvecs[:, iev]
            vecq = np.array([0, qy, qz], dtype=np.complex128)

            evr, evi = np.real(ev), np.imag(ev)
            is_pureimag = np.abs(evr) < tol

            if is_pureimag:
                v_is_fwds[iev] = evi > 0
                v_is_pureimag[iev] = True
                v_poynting_Sav[iev, :] =self._layer.poynting_vector_Sav(Om, Vp, vecq, evec) # This is one of the time critical parts
            else:
                v_is_growing[iev] = evr > 0
                v_is_decaying[iev] = evr < 0



        ngrow = np.sum(v_is_growing)
        ndecay = np.sum(v_is_decaying)
        nfwd = np.sum(v_is_fwds)

        # take all the growing/decaying, and then as many outwards S modes as we need to fill it up
        if self._is_front_side:
            sign_Sy = -1
            ev_good.extend(np.nonzero(v_is_growing)[0])  # take all the growing
        else:  # back bc
            sign_Sy = 1
            ev_good.extend(np.nonzero(v_is_decaying)[0])

        missing_ev = hdim - len(ev_good)

        # fill up the rest of the good slots with outward (leftward) going propagating waves
        for i in np.nonzero(v_is_pureimag)[0]:
            if (
                #v_is_pureimag[i] and
                v_poynting_Sav[i, 1] * sign_Sy > 0
            ):  #  Sy is outwards
                ev_good.append(i)
                missing_ev -= 1
                if missing_ev == 0:
                    break

        if len(ev_good) != hdim:
            print(
                f"Unexpected number of good eigenvalues {len(ev_good)} in semi-infinite layer:"
            )
            print(
                f" Om:{Om:.4f}, Vp:{Vp:.4f}, (ngrow={ngrow}, ndecay={ndecay}, nfwd={nfwd}, ngood={len(ev_good)})."
            )
            for vprop, nm in (
                (v_is_growing, "growing"),
                (v_is_decaying, "decaying"),
                (v_is_pureimag, "pureimag"),
            ):

                for iev in vprop.nonzero()[0]:
                    Sx, Sy, Sz = v_poynting_Sav[iev, :] * 100
                    print(
                        f"  ev[{iev}] = {eigvals[iev]: 12.5f} {nm},  Svec=({Sx:7.3f},{Sy:7.3f},{Sz:7.3f})"
                    )

        # order the eigensolutions

        ev_bad = list(
            set(range(dim)) - set(ev_good)
        )  # the bads are the values {0,1,2,3,4,5} minus the known good indices
        ev_order = ev_good + ev_bad

        vEvals = eigvals[ev_order]
        mVeigs = eigvecs[:, ev_order]
        return vEvals, mVeigs

    def suggested_L_left(self, Om: float, Vp: float) -> float:
        """Suggest a finite truncation length for a left semi-infinite side.

        Uses the smallest decay length among admissible modes; returns a few
        decay lengths for robust profile plotting. Value is in normalised units.
        """
        if not self.is_semi_infinite():
            return 0.0

        vEvals, mVeigs = self.analyse_semiinfinite_eigenspace(Om, Vp)
        hdim = len(vEvals) // 2
        # find smallest decay length
        lamrs = np.real(vEvals[:hdim])  # only the kept evs
        lammin = np.min(lamrs)

        if not np.isclose(lammin, 0):
            L = 4.0 / abs(lammin)  # 4 decay lengths
        else:
            L = 5.0  # arbitrary

        if L > 100:  # something has gone wrong, probably a nondecaying eigenvalue
            L = 5.0
        return L

    def suggested_L_right(self, Om: float, Vp: float) -> float:
        """Suggest a finite truncation length for a right semi-infinite side."""
        return self.suggested_L_left(Om, Vp)

    def null_vector_to_field_psi0(self, v_null: NDArray[np.complex128]) -> NDArray[np.complex128]:
        """Map a nullspace vector into a physical state at y=0.

        For SEMI_INFINITE, the null vector are coefficients of selected eigenvectors.
        For VACUUM/SHORT/CHARGE_FREE, the null vector corresponds to displacements
        directly and algebraic constraints are applied.
        """
        hdim = len(v_null)
        dim = 2 * hdim

        if self.is_semi_infinite():                     # null vector is coeffs of eigenvectors
            assert self._semiinf_evecs is not None, "Semi-infinite eigenspace not analysed yet." # helps with typing
            v_psi0 = self._semiinf_evecs[:, :hdim] @ v_null
        else:
            v_psi0 = np.zeros(dim, dtype=np.complex128)
            v_psi0[hdim:] = v_null                      # null vector is literally initial displacement field

            if self.is_short():
                v_psi0 = self.Sswap @ v_psi0
            elif self.is_charge_free():
                v_psi0 = self.tau_phi_plus @ v_psi0

        return v_psi0


class ElasticStack:
    """Transfer-matrix solver for elastic waveguides composed of layers.

    Parameters
    ----------
    layers
        Sequence of ``ElasticLayer`` objects (front to back).
    bcs
        Two boundary conditions, for the front (left) and back (right) sides.

    Notes
    -----
    All computations are performed in normalised units. Public APIs that accept
    SI values are explicitly named with ``*_SI`` parameters.
    """
    def __init__(self, layers: Sequence[ElasticLayer], bcs: Sequence[ElasticBoundaryCondition]):
        self._def_Vmin = 0.01  # km/s
        self._def_Vmax = 10.0  # km/s

        self._n_layers = len(layers)
        # Store as lists internally to support list-specific operations elsewhere
        self._layers = list(layers)
        self._bcs = list(bcs)

        self.check_bcs()
        self.analyse_layers()

        # row swap matrix for shorted bc
        self.Sswap = np.eye(8, dtype=np.float64)
        row8 = self.Sswap[7, :].copy()
        self.Sswap[7, :] = self.Sswap[3, :]
        self.Sswap[3, :] = row8

        # charge-free BC matrices
        self.tau_phi_plus = np.eye(8, dtype=np.float64)
        self.tau_phi_minus = np.eye(8, dtype=np.float64)


    def check_bcs(self) -> None:
        """Validate boundary-condition list and basic configuration invariants."""
        if len(self._bcs) != 2:
            raise ValueError(
                "Boundary conditions list must have exactly two elements (top and bottom)."
            )

        if (
            self._bcs[0].is_vacuum()
            and self._bcs[1].is_vacuum()
            and self._n_layers == 0
        ):
            raise ValueError(
                "At least one layer is required with vacuum boundary conditions on both sides."
            )

    def analyse_layers(self) -> None:
        """Determine stack dimensionality and prepare slices based on piezo presence."""
        self._bcs[0].set_is_front()
        self._bcs[1].set_is_back()

        # check if any piezo layers or bcs
        self.dim = 6
        all_layers = self._layers + [self._bcs[0], self._bcs[1]]
        for layer in all_layers:
            if layer.is_piezo():
                self.dim = 8
                break

        self.hdim = self.dim // 2

        self.slice_Tn = slice(0, 3)
        self.slice_uvec = slice(self.hdim, self.hdim+3)

        self._prepare_quadrant_slice()

    def _prepare_quadrant_slice(self) -> None:
        """Choose which quadrant of the full G matrix is used for the determinant.

        The chosen submatrix depends on whether each side is vacuum-like or
        semi-infinite. This avoids ill-conditioning from including constrained
        rows/columns that do not affect mode existence.
        """

        bc0_vaclike = not self._bcs[0].is_semi_infinite()
        bc1_vaclike = not self._bcs[1].is_semi_infinite()

        # build slices for 4 quadrants of the full G matrix
        sltop = slice(0, self.hdim)
        slbot = slice(self.hdim, self.dim)
        slice_B = (sltop, slbot)
        slice_D = (slbot, slbot)
        slice_A = (sltop, sltop)
        slice_C = (slbot, sltop)

        match (bc0_vaclike, bc1_vaclike):
            case (True, True):
                the_slice = slice_B # G12 = B
            case (True, False):
                the_slice = slice_D  # G22 = D
            case (False, True):
                the_slice = slice_A  # G11 = A
            case (False, False):
                the_slice = slice_C  # G21 = C
        self.the_quadrant_slice = the_slice


    def _build_profile_domain_y(
        self, ny_per_layer: int = 50, L_left: float = 0, L_right: float = 0
    ) -> Tuple[NDArray[np.float64], NDArray[np.int32]]:
        """Make vectors ``v_y`` and ``v_ilayer`` for the profile domain.

        Parameters
        ----------
        ny_per_layer
            Number of samples per physical layer (uniform within layers).
        L_left, L_right
            Optional semi-infinite truncation lengths (normalised units).

        Returns
        -------
        v_y : ndarray
            Monotonic y-positions (normalised) spanning left BC through layers to right BC.
        v_ilayer : ndarray[int]
            Layer index per sample; -1 for left semi-infinite region, -2 for right.
        """

        # get semi-inf borders
        if self._bcs[0].is_semi_infinite():  # add a bit of space for decay
            L_left = L_left if L_left > 0 else 5.0
        else:
            L_left = 0.0
        if self._bcs[1].is_semi_infinite():  # add a bit of space for decay
            L_right = L_right if L_right > 0 else 5.0
        else:
            L_right = 0.0

        n_y = self._n_layers * ny_per_layer

        if L_left > 0:
            n_y += ny_per_layer
        if L_right > 0:
            n_y += ny_per_layer

        n_y += 1  # to include endpoint

        v_y = np.zeros(n_y)  # physical location
        # index of layer (-1 for left semi-inf, -2 for right semi-inf)
        v_ilayer = np.zeros(n_y, dtype=np.int32)

        iyoff = 0
        ylo = 0.0
        yhi = 0.0
        dy = 0.0
        if L_left > 0:
            ylo -= L_left
            dy = L_left / ny_per_layer
            v_y[iyoff:ny_per_layer] = np.linspace(ylo, yhi - dy, ny_per_layer)
            v_ilayer[iyoff:ny_per_layer] = -1
            iyoff += ny_per_layer
            ylo += L_left

        layer_lengths = [layer.L for layer in self._layers]
        if L_right > 0:
            layer_lengths.append(L_right)

        for ilay, L in enumerate(layer_lengths):
            yhi += L
            dy = L / ny_per_layer
            v_y[iyoff : iyoff + ny_per_layer] = np.linspace(ylo, yhi - dy, ny_per_layer)
            v_ilayer[iyoff : iyoff + ny_per_layer] = ilay
            iyoff += ny_per_layer
            ylo += L

        v_y[-1] = v_y[-2] + dy
        v_ilayer[-1] = v_ilayer[-2]

        if L_right > 0:  # overwrite code for
            v_ilayer[-(ny_per_layer + 1) :] = -2

        return v_y, v_ilayer

    def update_tau_phi_matrices(self, Om: float, Vp: float) -> None:
        """Build the charge-free boundary-condition coupling matrices.

        For piezoelectric problems the charge-free condition introduces a
        coupling between stress and potential; the coefficient is q=Omega/V.
        """
        q = Om / Vp
        #self.tau_phi_plus = np.eye(8, dtype=np.float64)
        #self.tau_phi_minus = np.eye(8, dtype=np.float64)
        self.tau_phi_plus[3,7] = q
        self.tau_phi_minus[3,7] = -q




    def _find_full_transfer_matrix(
        self, Om: float, Vp: float
    ) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        """Assemble the full transfer matrix across the stack and a quadrant.

        Returns the full G matrix and the selected quadrant used for det-based
        mode finding, including contributions from layers and boundary
        conditions.
        """
        # build layer parts
        mat_G = np.eye(self.dim, dtype=np.complex128)
        self.update_tau_phi_matrices(Om, Vp)

        for layer in self._layers:
            mat_T, mat_expTL = layer.find_layer_transfer_matrix(Om, Vp)
            mat_G = mat_expTL @ mat_G  # type: ignore[arg-type]  # pylance worries that mat_expTL might be None

        # add boundary layer parts
        for ibcs, bc in enumerate(self._bcs):

            if bc.is_semi_infinite():

                # sort eigenvectors into growing and decaying/forwards and backwards modes
                v_evals, m_eigvecs = bc.analyse_semiinfinite_eigenspace(Om, Vp)

                if bc._is_front_side:
                    mat_G = mat_G @ m_eigvecs
                else:
                    mat_G = npla.inv(m_eigvecs) @ mat_G
            elif bc.is_charge_free():
                if bc._is_front_side:
                    mat_G = mat_G @ self.tau_phi_plus
                else:
                    mat_G = self.tau_phi_minus @ mat_G
            elif bc.is_short():
                if bc._is_front_side:
                    mat_G = mat_G @ self.Sswap
                else:
                    mat_G = self.Sswap @ mat_G

        mat_G_ij = mat_G[self.the_quadrant_slice]


        mat_G_ij /= np.max(np.abs(mat_G_ij))  # scale to avoid overflow in det


        # scale to avoid overflow in det and minimise meaningless phase changes
        # maxG_flatindex = np.argmax(np.abs(mat_G_ij))
        # maxG_loc = np.unravel_index(maxG_flatindex, mat_G_ij.shape)
        # maxG = mat_G_ij[maxG_loc]
        # mat_G_ij /= maxG


        return mat_G, mat_G_ij


    def _find_determinant_minima(self,  Om: float, v_Vp: NDArray[np.float64],
                                 show_scan_plot: bool, show_progress: bool, prefix: str,
                                 label_scanstep: str, pllabel_scanstep: str,
                                 vBulks=None) -> Tuple[NDArray[np.complex128],
                                                       List[Tuple[float, float]],
                                                       float]:
        """Scan det vs phase velocity and bracket local minima.

        Returns the complex det samples, the list of bracketed minima (value,
        det magnitude), and the maximum det magnitude (for normalisation).
        """


        v_det: NDArray[np.complex128] = np.zeros(len(v_Vp), dtype=np.complex128)

        for ivP, vP in enumerate(v_Vp):
            mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
            t_det = npla.det(mat_G_ij) # 1j # can be useful to pull out factor i for plots sometimes

            # # tune the phase
            # if np.abs(np.imag(t_det)) > 1e-6 * np.abs(t_det):
            #     ang = np.angle(t_det)
            #     t_det *= np.exp(-1j * ang)

            v_det[ivP] = t_det

        if show_scan_plot:
            plot_det_scan(Om, v_Vp, v_det, prefix=prefix, pllabel=pllabel_scanstep, Vbulks=vBulks)

        l_braks = bracket_minima(v_Vp, np.abs(v_det))

        print(f"\nAt {label_scanstep}:")
        if vBulks is not None:
            print(f"  Vbulks: {vBulks[2]:.6f},  {vBulks[1]:.6f},  {vBulks[0]:.6f}")

        show_brack = False and show_progress
        if show_brack:
            print(f"\n  Found {len(l_braks)} brackets:")

            for br in l_braks:
                (xm1, fm1), (x0, f0), (xp1, fp1) = np.array(
                    br
                ).tolist()  # awkward but converts to plain floats
                print(
                    f"    [{xm1:.5f}, {x0:.5f}, {xp1:.5f}] with f=[{fm1:.3e}, {f0:.3e}, {fp1:.3e}]"
                )

        if len(l_braks) == 0:
            plot_det_scan(Om, v_Vp, v_det, prefix=prefix, pllabel=pllabel_scanstep, Vbulks=vBulks)
            raise BracketError(
                "No minima found in determinant scan - no modes found. Determinant scan plot dumped."
            )

        def minfunc(vP: float) -> float:
            mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
            return np.abs(npla.det(mat_G_ij))

        l_mins = find_bracketed_minima(minfunc, l_braks)
        if show_progress:
            print(f"  Found {len(l_mins)} minima:")
            for v, f in l_mins:
                bulkdist = ""
                if vBulks is not None:
                    dbulk = np.min(np.abs(vBulks-v))
                    bulkdist= f", dist to bulk mode = {dbulk:.3e}"
                print(f"    Vp={v:.9f} km/s, Det(f)={f:.4e}{bulkdist}")

        max_det = np.max(np.abs(v_det))

        return v_det, l_mins, max_det

    def _filter_determinant_minima(self, l_mins: List[Tuple[float, float]],
                                 detmin_thresh: float, max_det: float, show_progress: bool) -> List[Tuple[float, float]]:
        """Filter spurious minima using a relative threshold and (optionally) SVD checks."""
        # remove any minima with f above threshold or with multiple singular values as probably not real modes
        fmin_norm_thresh = detmin_thresh * max_det # np.max(np.abs(v_det))
        l_mins_V_keep = []
        l_mins_bad_threshold = []
        l_mins_bad_singval = []

        for v, f in l_mins:
            if f > fmin_norm_thresh:
                l_mins_bad_threshold.append((v, f))
            # elif self.count_singular_values(Om, v) != 1:
            #    l_mins_bad_singval.append((v,f))  # could be zero or two singvals
            else:
                l_mins_V_keep.append((v , f))

        n_dropped = len(l_mins_bad_threshold) + len(l_mins_bad_singval)
        if show_progress and n_dropped > 0:
            print(f"  Dropping {n_dropped} modes:")
            for v, f in l_mins_bad_singval:
                print(f"    Vp={v:.3f} m/s has multiple small singular values.")
            for v, f in l_mins_bad_threshold:
                print(
                    f"    Vp={v:.3f} has relative det minimum f={f:.3e} above detmin_thresh={fmin_norm_thresh:.3e}:"
                )

        return l_mins_V_keep

    def _color_modes_by_poln(self, Omega_SI: float, l_V_SI: List[Tuple[float, float]]) -> NDArray[np.float64]:
        """Assign RGB colors to modes based on displacement polarization fractions."""
        full_mode_col = False
        v_mode_cols = np.zeros((len(l_V_SI), 3), dtype=np.float64)
        for iv, (vSI, f) in enumerate(l_V_SI):
            if full_mode_col:
                pol_vec = self.find_mode_polarization_over_profile_1d(Omega_SI, vSI)
            else:
                pol_vec = self.find_mode_polarization_at_y0(Omega_SI, vSI)

            v_mode_cols[iv, :] = nbpt.get_rgb_for_poln(*pol_vec)
        return v_mode_cols


    def _get_characteristic_bulk_velocities(self) -> Optional[List[float]]:
        """Return bulk wave velocities for a single half-space Rayleigh problem if applicable."""
        vBulks = None
        if not len(self._layers):
            if self._bcs[0].is_semi_infinite():
                vBulks = self._bcs[0]._layer.material.Vac_phase()  # type: ignore[arg-type]
            elif self._bcs[1].is_semi_infinite():
                vBulks = self._bcs[1]._layer.material.Vac_phase()  # type: ignore[arg-type]
        return vBulks

    def dispersion_relation_find_Vs_at_Omega(
        self,
        Omega_SI: float,
        Vmin_SI: float = 10,
        Vmax_SI: float = 10000.0,
        nVpts: int = 1000,
        find_mode_color: bool = False,
        prefix: str = "tmpdrscan",
        label: str = "",
        show_scan_plot: bool = False,
        show_progress: bool = False,
        detmin_thresh: float = 1e-5,
    ) -> Dict[str, Any]:
        """Find guided-mode phase velocities at a fixed frequency.

        Parameters
        ----------
        Omega_SI
            Angular frequency in rad/s.
        Vmin_SI, Vmax_SI
            Scan bounds for phase velocity in m/s.
        nVpts
            Number of scan points.
        find_mode_color
            If True, also compute a color per mode indicating polarization.
        prefix, label, show_scan_plot, show_progress, detmin_thresh
            Plotting/diagnostics controls and minimum-det threshold (relative).

        Returns
        -------
        dict
            Keys: "Vphase" (List[float]), "detmins" (List[float]),
            "poln_cols" (Optional ndarray Nx3), and "det_scan" (tuple of arrays).
        """

        nbt.check_magnitude(Omega_SI, "Omega", 1e8, 100e9 * 2 * np.pi)
        nbt.check_magnitude(Vmin_SI, "Vmin", 10, 20000)
        nbt.check_magnitude(Vmax_SI, "Vmax", 100, 20000)

        # move to normalised units here
        Om = Omega_SI / g_norms.f0
        Vmin = Vmin_SI / g_norms.v0
        Vmax = Vmax_SI / g_norms.v0

        if label:
            pllabel = label
        else:
            label = f"Nu={Om/(2*np.pi):.6f}"
            pllabel = rf"$\Omega/(2\pi)={Om/(2*np.pi):.4f}$ GHz"

        # if doing a rayleigh problem, add in the bulk velocities on the det scan for reference
        vBulks = self._get_characteristic_bulk_velocities()

        v_Vp = np.linspace(Vmin, Vmax, nVpts)
        # choose to space the velocities according to equal qz
        # v_Vp = 1.0/np.linspace(1/Vmin, 1/Vmax, n_Vp)


        v_det, l_mins, max_det = self._find_determinant_minima(Om, v_Vp,
                                                     show_scan_plot, show_progress, prefix,
                                                     label, pllabel, vBulks=vBulks)

        l_mins_V_keep = self._filter_determinant_minima(l_mins, detmin_thresh, max_det, show_progress)

        l_mins_V_SI_keep = [(v * g_norms.v0, f) for v, f in l_mins_V_keep]

        l_Vp_SI = [v * g_norms.v0 for v, f in l_mins_V_keep]
        l_detmins = [f for v, f in l_mins_V_keep]


        v_mode_cols = self._color_modes_by_poln(Omega_SI, l_mins_V_SI_keep) if find_mode_color else None


        res = {
            "Vphase": l_Vp_SI,
            "detmins": l_detmins,
            "poln_cols": v_mode_cols,
            "det_scan": (v_Vp, v_det),
        }

        return res

    def find_mode_polarization_at_y0(
        self, Omega_SI: float, vP_SI: float
    ) -> NDArray[np.complex128]:
        """Return the unit-norm displacement polarization at the interface (y=0)."""
        (v_psi0, rel_sing_val) = self.find_mode_psi0(Omega_SI / g_norms.f0, vP_SI / g_norms.v0)[0]

        uvec = v_psi0[self.slice_uvec]
        uvec /= np.linalg.norm(uvec)

        return uvec

    def find_mode_polarization_over_profile_1d(
        self,
        Omega_SI: float,
        vP_SI: float,
        ny_per_layer: int = 100,
        L_left: float = -1,
        L_right: float = -1,
    ) -> NDArray[np.float64]:
        """Compute polarization fractions integrated over the 1D profile.

        Returns an array of three non-negative values that sum to 1, representing
        relative energy in x/y/z displacement components.
        """

        modefunc_u, _= self.find_mode_profile_1d(Omega_SI, vP_SI)
        m_u_y = modefunc_u.m_uxyz[self.slice_uvec, :]  # displacement part
        poln_fracs = np.sum(np.abs(m_u_y) ** 2, axis=1)
        poln_fracs /= np.linalg.norm(poln_fracs)
        return poln_fracs

    def count_singular_values(self, Om: float, vP: float, thresh: float = 1e-6) -> int:
        """Return how many singular values of the determinant quadrant fall below a threshold."""
        mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)
        u, s, vh = npla.svd(mat_G_ij)
        nsmall = np.sum(s < thresh * np.max(s))  # Limited contrast precision in SVD
        return nsmall

    def find_mode_psi0(
        self,
        Om: float,
        vP: float,
        thresh_singval: float = 1e-4,
        show_progress: bool = False,
    ) -> List[Tuple[NDArray[np.complex128], float]]:
        """Find nullspace vectors at y=0 corresponding to candidate modes.

        Uses SVD of the determinant quadrant. Returns a list of (psi0, relative
        singular value) pairs; multiple entries may appear near apparent mode crossings.
        """

        mat_G, mat_G_ij = self._find_full_transfer_matrix(Om, vP)

        u, s, vh = npla.svd(mat_G_ij)

        if show_progress:
            print("Det check:", npla.det(mat_G_ij))
            print("SVD check:", s)

        nsmall = np.sum(
            s < thresh_singval * np.max(s)
        )  # Limited contrast precision in SVD

        if show_progress:
            print("Num small singular values:", nsmall)

        if nsmall == 0:
            print(f"Warning: largish singular values found for mode at Omega={Om:.3e}, Vs={vP:.3e}: not a mode?")
            print('S values:', s)
            v_psi0 = None

        if nsmall > 1:
            print(
                "Warning: Multiple small singular values found - may be non-genuine mode."
            )

        v_psi0s = []
        for isv in range(max(1, nsmall)):
            t_rel_singval = s[-1-isv] / np.max(s)
            v_null = vh.conj().T[
                :, -1-isv
            ]  # take the last column of Vh (smallest singular value)

            v_psi0 = self._bcs[0].null_vector_to_field_psi0(v_null)
            vecu = v_psi0[self.slice_uvec]
            v_psi0 /= np.max(np.abs(vecu))  # normalise max displacement to 1
            v_psi0s.append((v_psi0, t_rel_singval))
            if show_progress:
                with np.printoptions(precision=4, suppress=True, linewidth=200):
                    print(
                        f"Initial field at y=0:\n null vec: {v_null}\n psi0:  {v_psi0}"
                    )

        return v_psi0s

    def find_mode_profile_1d(
            self,
            Omega_SI: float,
            Vp_SI: float,
            ny_per_layer: int = 100,
            L_left: float = -1,
            L_right: float = -1,
            singval: int = 0,  # access multiple modes if they show up
            singval_thresh: float = 1e-4,
        ) -> Tuple["ModeFunction1D", "ModeFunction1D"]:
        """Compute 1D mode profiles for displacement and (normal) stress.

        Returns two ModeFunction1D objects: displacement components and normal
        stresses (Tyy, Txy, Tzy). The profiles are normalised for convenient plotting.
        """
        # find nullspace vector and initial field at y=0
        Om = Omega_SI / g_norms.f0
        vP = Vp_SI / g_norms.v0

        try:
            (v_psi0, rel_singval) = self.find_mode_psi0(Om, vP,
                                                        thresh_singval=singval_thresh)[singval]
        except ValueError as e:
            raise ValueError(
                f"Could not find mode at Omega={Om:.3e} ({Om:.3f} norm), Vs={vP:.3f} ({vP:.3f} norm): {e}"
            ) from e

        # find profile through stack

        # user can override default L_left and L_right for semi-infinite layers
        L_left = L_left if L_left > 0 else self._bcs[0].suggested_L_left(Om, vP)
        L_right = L_right if L_right > 0 else self._bcs[1].suggested_L_right(Om, vP)

        v_y, v_ilayer = self._build_profile_domain_y(ny_per_layer, L_left, L_right)

        # fill profile
        n_y = len(v_y)
        m_psi_y = np.zeros((self.dim, n_y), dtype=np.complex128)

        # precompute layer transfer matrices
        has_left = v_ilayer[0] == -1
        has_right = v_ilayer[-1] == -2
        if has_left:
            mat_Tleft, _ = self._bcs[0].find_layer_transfer_matrix(Om, vP) # type: ignore[attr-defined]

        v_mat_Tlayer = []
        for ilay, layer in enumerate(self._layers):
            mat_Tlayer, _ = layer.find_layer_transfer_matrix(Om, vP)
            v_mat_Tlayer.append(mat_Tlayer)

        if has_right:
            mat_Tright, _ = self._bcs[1].find_layer_transfer_matrix(Om, vP) # type: ignore[attr-defined]

            # we just put the right hand matrix onto the end of the list. Doesn't need to be handled specially
            v_mat_Tlayer.append(mat_Tright)

        y_layer_start = 0
        ilayer_last = 0
        v_psi_layer_0 = v_psi0.copy()  # initial field at start of current layer

        for iy, y in enumerate(v_y):

            ilayer = v_ilayer[iy]

            # We've entered a new layer, update the starting vector for the layer by propagating over the full previous layer.
            # We don't do it for the Left layer to zero layer transfer since nothing has happened
            if ilayer > 0 and ilayer != ilayer_last:
                L_last_layer = y - y_layer_start
                G_last_layer = spla.expm(v_mat_Tlayer[ilayer_last] * L_last_layer)
                v_psi_layer_0 = (
                    G_last_layer @ v_psi_layer_0
                )  # push forward the start of the layer by one layer
                # mat_Glayers_done = G_last_layer @ mat_Glayers_done
                y_layer_start = y
                ilayer_last = ilayer

            if (
                ilayer == -1
            ):  # in the first layer. Back-propagate from the y=0 point (y is negative)
                m_psi_y[:, iy] = spla.expm(mat_Tleft * y) @ v_psi0 # type: ignore[operator]
            else:
                Dy = y - y_layer_start  # distance into current layer
                mat_T = (
                    v_mat_Tlayer[ilayer] if ilayer > 0 else v_mat_Tlayer[-1]
                )  # if we are in right layer (ilayer=-2), then use the last layer!
                m_psi_y[:, iy] = spla.expm(mat_T * Dy) @ v_psi_layer_0


        # TODO: should probably unscale the position v_y

        v_y_SI = v_y * g_norms.x0

        m_uxyz = m_psi_y[self.slice_uvec, :].T.copy()
        m_T6 = m_psi_y[self.slice_Tn, :].T.copy()  # only normal stresses Tyy, Txy, Tzy

        # to make physically realistic values  will treat displacements as being in nm, so stresses come down by 1000
        m_T6 /= 1000

        mode_prof_u = ModeFunction1D(Omega_SI, Vp_SI, v_y_SI, m_uxyz, relative_sing_val=rel_singval)
        mode_prof_T6 = ModeFunction1D(Omega_SI, Vp_SI, v_y_SI, m_T6, relative_sing_val=rel_singval)

        mode_prof_T6.set_labels(
            dvar=QuantLabel("Normal stress", r"$\vec T_i$", r"$|\vec T_i|$", "GPa")
        )

        return mode_prof_u, mode_prof_T6


def plot_det_scan(
    Omega: float,
    v_Vp: NDArray[np.float64],
    v_det: NDArray[np.complex128],
    prefix: str = "tmp",
    pllabel: str = "",
    Vbulks = None,
) -> None:
    """Plot the determinant scan used to bracket guided modes and save to disk."""
    fig, ax = plt.subplots()

    ax.plot(
        v_Vp, np.abs(v_det), "-+", label=r"$|\Delta|$", markersize=1.0, linewidth=0.5
    )


    ax2=ax.twinx()
    vdr, vdi = np.real(v_det), np.imag(v_det)
    ax2.plot(v_Vp, nbt.signed_log10(vdr,power=4), '-o', linewidth=0.5, color='r',
             label=r"Re($\Delta$)", markersize=1.0)
    ax2.plot(v_Vp, nbt.signed_log10(vdi,power=4), '-o', linewidth=0.5, color='g',
             label=r"Im($\Delta$)", markersize=1.0)

    # ax2.plot(v_Vp, np.angle(v_det), '-x', markersize=1.0, linewidth=0.5, color='C3', label=r"Arg($\Delta$)")

    ax.set_yscale("log")
    ax.set_xlabel(r"$V$ [km/s]", fontsize=10)
    ax.set_ylabel(r"$|T|$ ", fontsize=10)
    ax.set_title(pllabel, fontsize=10)
    ax.tick_params(axis="both", labelsize=10)
    ax.legend(fontsize=10)
    ax.set_xlim(np.min(v_Vp), np.max(v_Vp))
    if Vbulks is not None:
        for Vb in Vbulks:
            ax.axvline(Vb, color='gray', linestyle='--', linewidth=0.8)

    nbpt.save_and_close_figure(fig, prefix + "_det.png")


def plot_rayleigh_det_scan(
    v_Vp: NDArray[np.float64],
    v_det: NDArray[np.complex128],
    Vs: float,
    Vl: float,
    material: Any,  #: materials.Material
) -> None:
    """Plot determinant components for the Rayleigh half-space problem and save."""
    fig, ax = plt.subplots()
    ax.plot(v_Vp * 1e-3, np.real(v_det), "o", label=r"$\Delta_r$", markersize=0.5)
    ax.plot(v_Vp * 1e-3, np.imag(v_det), "x", label=r"$\Delta_i$", markersize=0.5)
    ax.plot(v_Vp * 1e-3, np.abs(v_det), "+", label=r"$|\Delta|$", markersize=0.5)
    ax.set_xlabel(r"$V$ [km/s]", fontsize=10)
    ax.set_ylabel(r"$|T|$ ", fontsize=10)
    ax.tick_params(axis="both", labelsize=10)
    ax.legend(fontsize=10)

    nbpt.save_and_close_figure(fig, "ray_det.png")


def find_Rayleigh_velocity_anisotropic(material: Any) -> float:  #: materials.Material
    """Find Rayleigh velocity for arbitrary material along the zhat=(0,0,1) direction with surface normal along yhat=(0,1,0)."""

    bcs = [
        ElasticBoundaryCondition(ElBCType.SEMI_INFINITE, material),
        ElasticBoundaryCondition(ElBCType.VACUUM),
    ]
    l_layers = []

    Omega = 2 * np.pi * 1.0 * SI_GHz
    estack = ElasticStack(l_layers, bcs)
    res = estack.dispersion_relation_find_Vs_at_Omega(
        Omega, Vmin_SI=1 * SI_kmps, Vmax_SI=8.0 * SI_kmps
    )
    Vp = res["Vphase"]

    if Vp:
        return Vp[0] #[0]  # return first found mode
    else:
        return 0






# Coefficient matrices defined in the layers3.nb mathematica file
def _get_layer_stiffness_submatrices(
    cs: NDArray[np.float64], use_4d: bool = False
) -> Tuple[
    NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]
]:
    """Build 3x3 (or 4x4) submatrices of the constitutive tensor for a layer.

    Parameters
    ----------
    cs
        Normalised 7x7 stiffness matrix in Voigt-like indexing (unit-indexed).
    use_4d
        When True, embed the 3x3 blocks into 4x4 blocks for an 8-component
        (piezo-aware) state vector; the extra row/column are padded with zeros.

    Returns
    -------
    mM1, mM2, mL1, mL2 : ndarray
        The submatrices used to construct the transfer matrix blocks for the
        non-piezoelectric case. Shapes are 3x3 for use_4d=False, otherwise 4x4.
    """

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
        mM1_4d = np.zeros((4, 4), dtype=np.float64)
        mM2_4d = np.zeros((4, 4), dtype=np.float64)
        mL1_4d = np.zeros((4, 4), dtype=np.float64)
        mL2_4d = np.zeros((4, 4), dtype=np.float64)
        mM1_4d[:3, :3] = mM1
        mM2_4d[:3, :3] = mM2
        mL1_4d[:3, :3] = mL1
        mL2_4d[:3, :3] = mL2

        mM1, mM2, mL1, mL2 = mM1_4d, mM2_4d, mL1_4d, mL2_4d

    return mM1, mM2, mL1, mL2


def _get_layer_stiffness_submatrices_piezo(
    cs: NDArray[np.float64], e_iJ: NDArray[np.float64], perm_ij: NDArray[np.float64]
) -> Tuple[
    NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]
]:
    """Build 4x4 submatrices for a piezoelectric layer.

    Parameters
    ----------
    cs
        Normalised 7x7 stiffness matrix (unit-indexed).
    e_iJ
        Normalised 4x7 piezoelectric coupling tensor (unit-indexed; i in 0..3, J in 0..6).
    perm_ij
        Relative permittivity tensor (4x4, unit-indexed).

    Returns
    -------
    mM1, mM2, mL1, mL2 : ndarray
        The 4x4 piezo-aware submatrices to construct transfer matrix blocks for
        piezoelectric media.
    """

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
