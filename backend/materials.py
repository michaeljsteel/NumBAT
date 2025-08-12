# materials.py is a subroutine of NumBAT that defines Material structure,
# these represent dispersive lossy refractive indices and possess
# methods to interpolate n from tabulated data.

# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

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


# - build eps from eps_xx, eps_yy, eps_zz, _xy, ...
# - load e_ijk rather than d_ijk, and construct d_ijk from it
# - print out e_ijk
# - print optical props


import copy
import json
import math
import pathlib
import re
import traceback

import numpy as np
import numpy.linalg as npla
import scipy.optimize as sciopt

from nbtypes import CrystalGroup, unit_x, unit_y, unit_z, SI_permittivity_eps0
import reporting
import voigt
from bulkprops import solve_christoffel, solve_dielectric
#import numbattools as nbtools
from crystalsyms import construct_crystal_for_symmetry_class


import plotting.materials as plms


class BadMaterialFileError(Exception):
    """Exception raised for errors in material file parsing or validation."""
    pass


g_material_library = None


def make_material(s):
    global g_material_library

    if g_material_library is None:
        g_material_library = MaterialLibrary()

    return g_material_library.get_material(s)


class MaterialLibrary:
    """
    Manages the loading and retrieval of material definitions from JSON files in the material_data directory.
    """
    def __init__(self):
        self._material_data_path = ""
        self._materials = {}

        # identify mat data directory:  backend/material_data
        this_dir = pathlib.Path(__file__).resolve().parent
        self._material_data_path = this_dir / "material_data"

        self._load_materials()

    def get_material(self, matname):
        """
        Retrieve a Material object by name.

        Args:
            matname (str): Name of the material to retrieve.
        Returns:
            Material: The requested material object.
        Raises:
            SystemExit: If the material is not found.
        """
        try:
            mat = self._materials[matname]
        except KeyError:
            reporting.report_and_exit(
                f"Material {matname} not found in material_data folder.\nEither the material file is missing or the name field in the material file has been incorrectly specified."
            )

        return mat

    def _load_materials(self):
        for fname in pathlib.Path(self._material_data_path).glob("*.json"):
            json_data = None
            with open(fname, "r") as fin:
                # read whole file and remove comments
                s_in = "".join(fin.readlines())
                s_in = re.sub(r"//.*\n", "\n", s_in)

                try:
                    json_data = json.loads(s_in)
                except Exception as err:
                    traceback.print_exc()  # TODO: why this traceback?
                    reporting.report_and_exit(
                        f"JSON parsing error: {err} for file {fname}"
                    )

            try:
                mat_name = json_data["material_name"]
            except KeyError:
                reporting.report_and_exit(
                    f"Material file {fname} has no 'material_name' field."
                )

            if mat_name in self._materials:
                reporting.report_and_exit(
                    f"Material file {fname} has the same name as an existing material {mat_name}."
                )

            try:
                new_mat = Material(json_data, fname)
            except BadMaterialFileError as err:
                reporting.report_and_exit(str(err))

            self._materials[mat_name] = new_mat


class PiezoElectricProperties:
    """
    Encapsulates piezoelectric tensor properties and their manipulation for a material.
    Handles loading, rotation, and conversion between different piezoelectric tensor forms.
    """
    def __init__(self, d_piezo, stiffness_cE_IJ):
        """
        Initialize piezoelectric properties from a parameter dictionary and a stiffness tensor.

        Args:
            d_piezo (dict): Piezoelectric parameters.
            stiffness_cE_IJ (VoigtTensor4_IJ): Stiffness tensor at fixed electric field.
        """
        # See Auld vol 1, section 8.B for these definitions.

        # cE denotes stiffness tensor at fixed electric field E
        # epsT denotes dielecric tensor at fixed stress T

        self.active = d_piezo.get("piezo_active", 0)

        # incoming base tensors
        self._tens_cE_IJ = stiffness_cE_IJ.copy()

        # Units of [Coulomb/Newton], order 1e-12
        self._tens_d_iJ = voigt.VoigtTensor3_iJ(
            "piezo_d",
            "Strain piezo coefficient",
            "d_iJ",
            ("pC/N", 1e-12),
            transforms_like_piezo_d=True,
        )

        self._tens_e_iJ = voigt.VoigtTensor3_iJ(
            "piezo_e",
            "Stress piezo coefficient",
            "e_iJ",
            ("C/m^2", 1.0),
            transforms_like_piezo_d=False,
        )

        self._tens_relepsS_ij = voigt.PlainTensor2_ij(
            "piezo_epsS",  # epsT_ij - d_iJ cE_JI dbar_Ij
            physical_name="Constant strain dielectric tensor",
            physical_symbol="eps^S_ij",
            unit=None,
            is_complex=False,
            is_symmetric=True,
        )

        self._tens_relepsT_ij = voigt.PlainTensor2_ij(
            "piezo_epsT",
            physical_name="Constant stress dielectric tensor",
            physical_symbol="eps^T_ij",
            unit=None,
            is_complex=False,
            is_symmetric=True,
        )

        self._load_piezo_d_IJ_from_json(d_piezo)

        self._load_diels_eps_ij(d_piezo)  # must come after loading of d_IJ

        self._store_original_tensors()

    def _load_diels_eps_ij(self, d_piezo):
        d_iJ = self._tens_d_iJ.value()
        cE_IJ = self._tens_cE_IJ.value()
        d_cE_d_ij = d_iJ @ cE_IJ @ d_iJ.T / SI_permittivity_eps0  # Auld 8.41

        if d_piezo.get("piezo_use_epsT", 0):
            self._tens_relepsT_ij.set_from_params(d_piezo)
            epsS_ij = self._tens_relepsT_ij.value() - d_cE_d_ij
            self._tens_relepsS_ij.set_from_matrix(epsS_ij)

        else:
            self._tens_relepsS_ij.set_from_params(d_piezo)
            epsT_ij = self._tens_relepsS_ij.value() + d_cE_d_ij
            self._tens_relepsT_ij.set_from_matrix(epsT_ij)

    def _load_piezo_d_IJ_from_json(self, d_piezo):
        if d_piezo["crystal_sym"] == "no sym":
            self._tens_d_iJ.set_from_params(d_piezo)
        else:
            self._fill_piezo_elt_mappings(d_piezo)

    def _fill_piezo_elt_mappings(self, d_piezo):
        sym_group = d_piezo["crystal_sym"]

        elts_indep = []
        elts_dep = []

        # first index of elts_dep is d_ijk multipliers, second is e_ijk multipliers
        if sym_group == "3m":
            elts_indep = "x5", "y2", "z1", "z3"
            elts_dep = {
                "x6": ("y2", -2.0, -1.0),
                "y1": ("y2", -1.0, -1.0),
                "y4": ("x5", 1.0, 1.0),
                "z2": ("z1", 1.0, 1.0),
            }

        elif sym_group == "4'3m":
            elts_indep = ("x4",)
            elts_dep = {"y5": ("x4", 1.0, 1.0), "z6": ("x4", 1.0, 1.0)}

        elif sym_group == "6mm":
            elts_indep = "x5", "z1", "z3"
            elts_dep = {"y4": ("x5", 1.0, 1.0), "z2": ("z1", 1.0, 1.0)}

        if d_piezo.get(
            "piezo_use_e", 0
        ):  # load e_iJ and convert to d_iJ using stiffness

            self._tens_e_iJ.set_from_structure_elts(
                d_piezo, elts_indep, elts_dep, depwhich=2
            )

            # d_iK = e_iJ . s_JK
            tens_s_compliance = npla.inv(self._tens_cE_IJ.value())
            self._tens_d_iJ.set_from_matrix(
                self._tens_e_iJ.value() @ tens_s_compliance
            )

        else:  # default: use the d_iJ parameters

            self._tens_d_iJ.set_from_structure_elts(
                d_piezo, elts_indep, elts_dep, depwhich=1
            )

            tens_cE = self._tens_cE_IJ.value()

            self._tens_e_iJ.set_from_matrix(
                self._tens_d_iJ.value() @ tens_cE
            )  # Auld 8.40

    def make_piezo_stiffened_stiffness_for_kappa(self, v_kap):
        """
        Returns the piezo-stiffened stiffness tensor for propagation along a given direction.
        See Auld 8.147.

        Args:
            v_kap (np.ndarray): Propagation direction (3-vector).
        Returns:
            VoigtTensor4_IJ: Piezo-stiffened stiffness tensor.
        """

        e_iL = self._tens_e_iJ.value()  # [0..2, 0..5]
        e_Kj = self._tens_e_iJ.value().T  # [0..5, 0..2]

        e_Kj_lj = e_Kj @ v_kap
        li_e_iL = v_kap @ e_iL

        l_eps_l = v_kap @ self._tens_relepsS_ij.value() @ v_kap * SI_permittivity_eps0
        corr_stiff = np.outer(e_Kj_lj, li_e_iL) / l_eps_l
        cE_KL_stiff = self._tens_cE_IJ.value() + corr_stiff
        # with np.printoptions(precision=4, suppress=True):
        #     print('\n\nkap:', v_kap,
        #           '\n eiL', e_iL, '\n eKj', e_Kj,
        #           '\n li_e_iL', li_e_iL, '\n e_Kj_lj', e_Kj_lj,
        #           '\n l_eps_l', l_eps_l, '\n corr_stiff\n', corr_stiff /1e10,
        #           '\n cE_KL_stiff\n', cE_KL_stiff /1e10,
        #           )
        # convert to a voigt4
        vt_cE_KL_stiff = self._tens_cE_IJ.copy()
        vt_cE_KL_stiff.set_from_matrix(cE_KL_stiff)

        return vt_cE_KL_stiff

    def rotate(self, matR):
        """
        Rotate all piezoelectric tensors by a given SO(3) rotation matrix.

        Args:
            matR (np.ndarray): 3x3 rotation matrix.
        """

        self._tens_cE_IJ.rotate(matR)

        self._tens_d_iJ.rotate(matR)
        self._tens_e_iJ.rotate(matR)

        self._tens_relepsS_ij.rotate(matR)
        self._tens_relepsT_ij.rotate(matR)

    def reset_orientation(self):
        """
        Restore all piezoelectric tensors to their original (unrotated) state.
        """

        self._tens_d_iJ = self._tens_d_iJ_orig.copy()
        self._tens_cE_IJ = self._tens_cE_IJ_orig.copy()
        self._tens_e_iJ = self._tens_e_iJ_orig.copy()
        self._tens_relepsS_ij = self._tens_relepsS_ij_orig.copy()
        self._tens_relepsT_ij = self._tens_relepsT_ij_orig.copy()

    def _store_original_tensors(self):
        """Save tensors to optionally restore after rotations."""

        self._tens_d_iJ_orig = copy.deepcopy(self._tens_d_iJ)
        self._tens_cE_IJ_orig = copy.deepcopy(self._tens_cE_IJ)
        self._tens_e_iJ_orig = copy.deepcopy(self._tens_e_iJ)
        self._tens_relepsS_ij_orig = copy.deepcopy(self._tens_relepsS_ij)
        self._tens_relepsT_ij_orig = copy.deepcopy(self._tens_relepsT_ij)

    def d_iJ_to_d_ijk(self):
        """
        Convert the d_iJ tensor to full d_ijk notation (not implemented).
        """
        pass

    def d_ijk_to_d_iJ(self):
        """
        Convert the d_ijk tensor to Voigt d_iJ notation (not implemented).
        """
        pass

    def __str__(self):
        """
        Return a formatted string representation of the piezoelectric properties.
        """
        # print('and now', self._tens_d_iJ)
        s = "\n Piezoelectric properties:"
        if not self.active:
            s += "\n  Piezo effects supported but disabled"
        else:
            s += "\n  Piezo effects enabled:"

            with np.printoptions(
                precision=4, floatmode="fixed", sign=" ", suppress=True
            ):
                s += "\n " + str(self._tens_d_iJ)
                s += "\n " + str(self._tens_e_iJ)
                s += "\n " + str(self._tens_relepsS_ij)
                s += "\n " + str(self._tens_relepsT_ij)

        return s


class Material(object):
    """
    Represents a waveguide material, including optical, elastic, and piezoelectric properties.
    Should be constructed via materials.get_material().
    """
    def __init__(self, json_data, filename):
        """
        Initialize a Material object from JSON data and filename.

        Args:
            json_data (dict): Material data loaded from JSON.
            filename (str): Path to the JSON file.
        """

        # a,b,c crystal axes according to standard conventions
        self._crystal_axes = []

        self._anisotropic = False
        self.nuPoisson = 0.0
        self.EYoung = 0.0
        self.rho = None

        self.refindex_n = 0

        # current rotated tensors
        self.total_matR = np.eye(3)
        self.stiffness_c_IJ = None
        self.viscosity_eta_IJ = None
        self.photoel_p_IJ = None
        self.optdiel_eps_ij = None

        # home crystal orientation tensors
        self._stiffness_c_IJ_orig = None
        self._photoel_p_IJ_orig = None
        self._viscosity_eta_IJ_orig = None
        self._optdiel_eps_ij_orig = None

        self._piezo = None

        self._parse_json_data(json_data, filename)

    def __str__(self):
        """
        Return a summary string for the material, including name, file, author, and date.
        """
        s = (
            f"Material: {self.chemical}\n"
            f"  File: {self.material_name}\n"
            f"  Source: {self.author}\n"
            f"  Date: {self.date}"
        )
        if len(self.comment):
            s += f"\nComment: {self.comment}"
        return s

    def piezo_supported(self):
        """
        Returns True if the material supports piezoelectric effects.
        """
        return self._piezo is not None

    def enable_piezoelectric_effects(self):
        """
        Enable piezoelectric effects for this material (if supported).
        """

        if self.piezo_supported():
            self._piezo.active = True

    def disable_piezoelectric_effects(self):
        """
        Disable piezoelectric effects for this material (if supported).
        """

        if self.piezo_supported():
            self._piezo.active = False

    def copy(self):
        """
        Return a deep copy of the material object.
        """

        return copy.deepcopy(self)

    def full_str(self, chop=False, chop_rtol=1e-10, chop_atol=1e-12):
        """
        Return a detailed string representation of the material, including all tensors and piezo properties.
        """
        s = str(self)

        s += f"\n  Crystal class:  {self.crystal_class.name}"
        s += f"\n  Crystal group:  {self.crystal_group}"

        s += "\n" + self.optdiel_eps_ij.as_str()

        s += "\n" + self.stiffness_c_IJ.as_str()
        s += "\n" + self.viscosity_eta_IJ.as_str()
        s += "\n" + self.photoel_p_IJ.as_str()

        if self.piezo_supported():
            s += "\n\n" + str(self._piezo)

        s += "\n"

        return s

    def optical_properties(self):
        """
        Return a string describing the optical properties of the material.
        """

        """Returns a string containing key elastic properties of the material."""
        #dent = "\n  "

        try:
            s = f"Optical properties of material {self.material_name}"

        except Exception as err:
            s = (
                f"Unknown/undefined optical parameters in material {self.material_name}"
                + str(err)
            )
        return s

    def elastic_properties(self, chop=False, chop_rtol=1e-10, chop_atol=1e-12):
        """
        Return a string describing the elastic properties of the material.
        """

        """Returns a string containing key elastic properties of the material."""

        dent = "\n  "
        try:
            s = f"Elastic properties of material {self.material_name}"
            s += dent + f"Density:        {self.rho:.3f} kg/m^3"
            s += dent + f"Ref. index:     {self.refindex_n:.4f} "

            s += dent + f"Crystal class:  {self.crystal_class.name}"
            s += dent + f"Crystal group:  {self.crystal_group}"

            if self.is_isotropic():
                s += (
                    dent
                    + f"c11:            {self.stiffness_c_IJ.mat[1, 1]*1e-9:.3f} GPa"
                )
                s += (
                    dent
                    + f"c12:            {self.stiffness_c_IJ.mat[1, 2]*1e-9:.3f} GPa"
                )
                s += (
                    dent
                    + f"c44:            {self.stiffness_c_IJ.mat[4, 4]*1e-9:.3f} GPa"
                )
                s += dent + f"Young's mod E:  {self.EYoung*1e-9:.3f} GPa"
                s += dent + f"Poisson ratio:  {self.nuPoisson:.3f}"
                s += dent + f"Velocity long.: {self.Vac_longitudinal():.3f} m/s"
                s += dent + f"Velocity shear: {self.Vac_shear():.3f} m/s"
                s += dent + f"Velocity Rayleigh: {self.Vac_Rayleigh():.3f} m/s"
            else:

                # find wave properties for z propagation
                v_phase, v_evecs, v_vgroup = solve_christoffel(
                    unit_z, self.stiffness_c_IJ, self.rho
                )

                with np.printoptions(
                    precision=4, floatmode="fixed", sign=" ", suppress=True
                ):
                    s += dent + " " + self.stiffness_c_IJ.as_str(chop=chop, rtol=chop_rtol, atol=chop_atol) + "\n"

                    for m in range(3):
                        vgabs = np.linalg.norm(v_vgroup[m])
                        s += (
                            dent
                            + f"Wave mode {m+1}: v_p={v_phase[m]:.4f} km/s,  |v_g|={vgabs:.4f} km/s,  "
                            + "u_j="
                            + str(v_evecs[:, m])
                            + ",  v_g="
                            + str(v_vgroup[m])
                            + " km/s"
                        )

            if self.piezo_supported():
                s += "\n"
                s += str(self._piezo)

        except Exception as err:
            s = (
                f"Unknown/undefined elastic parameters in material {self.material_name}:\n"
                + str(err)
            )
        return s

    def Vac_Rayleigh(self):
        """
        For an isotropic material, return the Rayleigh wave elastic phase velocity.
        """
        assert self.is_isotropic()

        Vs = self.Vac_shear()
        Vl = self.Vac_longitudinal()

        if not self.rho or self.rho == 0:  # Catch vacuum cases
            return 0.0
        else:
            return find_Rayleigh_velocity(Vs, Vl)

    def Vac_longitudinal(self):
        """
        For an isotropic material, return the longitudinal (P-wave) elastic phase velocity.
        """
        assert self.is_isotropic()

        if not self.rho or self.rho == 0:  # Catch vacuum cases
            return 0.0
        else:
            return math.sqrt(self.stiffness_c_IJ[1, 1] / self.rho)

    def Vac_shear(self):
        """
        For an isotropic material, return the shear (S-wave) elastic phase velocity.
        """
        assert self.is_isotropic()

        if not self.rho or self.rho == 0:  # Catch vacuum cases
            return 0.0
        else:
            return math.sqrt(self.stiffness_c_IJ[4, 4] / self.rho)

    def Vac_phase(self):
        """
        Return the three acoustic phase velocities for propagation along z for the current orientation.
        """
        vkap = np.array([0, 0, 1], dtype=np.float64)
        v_vphase, vecs, v_vgroup = solve_christoffel(
            vkap, self.stiffness_c_IJ, self.rho
        )
        return v_vphase

    def em_phase_index(self):
        """
        Return the three optical phase indices for propagation along z for the current orientation.
        """
        vkap = np.array([0, 0, 1], dtype=np.float64)
        m_diel = self.optdiel_eps_ij
        v_nphase, vecs, v_ngroup = solve_dielectric(vkap, m_diel)
        return v_nphase

    def has_elastic_properties(self):
        """
        Returns True if the material has elastic properties defined.
        """
        return self.rho is not None

    def _parse_meta_data(self, json_data, fname):
        # Name of this file, will be used as identifier and must be present
        self.material_name = json_data.get("material_name", "NOFILENAME")
        if self.material_name == "NOFILENAME":
            raise BadMaterialFileError(
                f"Material file {fname} has no 'material_name' field."
            )

        self.format = json_data.get("format", "NOFORMAT")
        if self.format == "NOFORMAT":
            raise BadMaterialFileError(f"Material file {fname} has no 'format' field.")

        if self.format != "NumBATMaterial-fmt-2.0":
            raise BadMaterialFileError(
                f"Material file {fname} must be in format 'NumBATMaterial-Fmt-2.0'."
            )

        self.chemical = json_data["chemical"]  # Chemical composition
        self.author = json_data["author"]  # Author of data
        # Year of data publication/measurement
        self.date = json_data["date"]
        # Source institution
        self.institution = json_data["institution"]
        # doi or, failing that, the http address
        self.doi = json_data["doi"]

        # general comment for any purpose
        self.comment = json_data.get("comment", "")

    def _parse_optical_data(self, json_data):
        self.optdiel_eps_ij = voigt.PlainTensor2_ij(
            "eps",
            "Optical dielectric constant",
            "eps_ij",
            is_complex=True,
            is_symmetric=True,
        )

        if "eps_xx_r" in json_data:  # use strain dielectric constants
            self.optdiel_eps_ij.set_from_params(json_data)
        else:
            Re_n = json_data["Re_n"]  # Real part of refractive index []
            Im_n = json_data["Im_n"]  # Imaginary part of refractive index []
            self.refindex_n = Re_n + 1j * Im_n  # Complex refractive index []
            self.optdiel_eps_ij.set_from_matrix(np.eye(3) * self.refindex_n**2)

    def _parse_crystal_types(self, json_data, fname):
        if "crystal_class" not in json_data:
            raise BadMaterialFileError(
                f"Material file {fname} has no 'crystal_class' field."
            )
        try:
            self.crystal_class = CrystalGroup[json_data["crystal_class"]]
        except ValueError as exc:
            print("Unknown crystal class in material data file")
            raise BadMaterialFileError(
                f"Unknown crystal class in material data file {fname}"
            ) from exc

        self.crystal_group = json_data.get("crystal_sym", "no sym")

    def _parse_json_data(self, json_data, fname):
        """
        Load material data from a JSON file.

        Args:
            json_data (dict): Material data.
            fname (str): Filename of the JSON file.
        """
        self._d_params = json_data  # Do without this?
        self.json_file = fname

        self._parse_meta_data(json_data, fname)
        self._parse_optical_data(json_data)

        if self.is_vacuum():  # no mechanical properties available
            return

        self.EYoung = None
        self.nuPoisson = None

        self.rho = json_data["s"]  # Density [kg/m3]

        self._parse_crystal_types(json_data, fname)


        self.stiffness_c_IJ = voigt.VoigtTensor4_IJ(
            "c", "Stiffness", "c_IJ", ("GPa", 1.0e9), transforms_like_stiffness=True
        )
        self.viscosity_eta_IJ = voigt.VoigtTensor4_IJ(
            "eta", "Viscosity", "eta_IJ", transforms_like_stiffness=True
        )
        self.photoel_p_IJ = voigt.VoigtTensor4_IJ(
            "p", "Photoelasticity", "p_IJ", transforms_like_stiffness=True  # TODO: what is correct form?
        )

        self.construct_crystal()


        dpiezo = {k: v for k, v in json_data.items() if k.startswith("piezo")}
        if dpiezo:
            dpiezo["crystal_sym"] = self.crystal_group
            self._piezo = PiezoElectricProperties(dpiezo, self.stiffness_c_IJ)

        self._store_original_tensors()

    def _store_original_tensors(self):
        self._stiffness_c_IJ_orig = self.stiffness_c_IJ.copy()
        self._photoel_p_IJ_orig = self.photoel_p_IJ.copy()
        self._viscosity_eta_IJ_orig = self.viscosity_eta_IJ.copy()
        self._optdiel_eps_ij_orig = self.optdiel_eps_ij.copy()

    def is_vacuum(self):
        """
        Returns True if the material is vacuum.
        """
        return self.chemical == "Vacuum"

    # (don't really need this as isotropic materials are the same)


    def construct_crystal_general(self):
        """
        Fill all tensor elements for a general anisotropic crystal from the parameter dictionary.
        """

        try:  # full anisotropic tensor components
            for I in range(1, 7):
                for J in range(1, 7):
                    #self.stiffness_c_IJ.read_from_json(i, j)
                    #self.photoel_p_IJ.read_from_json(i, j)
                    #self.viscosity_eta_IJ.read_from_json(i, j)

                    self.stiffness_c_IJ.set_elt_from_param_dict(I, J, self._d_params)
                    self.photoel_p_IJ.set_elt_from_param_dict(I, J, self._d_params)
                    self.viscosity_eta_IJ.set_elt_from_param_dict(I, J, self._d_params)



        except KeyError:
            reporting.report_and_exit(
                "Failed to load anisotropic crystal class in material data file {self.json_file}"
            )

    # TODO: make a property
    def set_refractive_index(self, nr, ni=0.0):
        """
        Set the complex refractive index for the material.

        Args:
            nr (float): Real part.
            ni (float): Imaginary part (default 0.0).
        """

        self.refindex_n = nr + 1j * ni

    def is_isotropic(self):
        """
        Returns True if the material is isotropic.
        """

        return not self._anisotropic

    def rotate(self, rot_axis_spec, theta, save_rotated_tensors=False):
        """
        Rotate the crystal axes and all tensors by theta radians about the specified axis.

        Args:
            rot_axis_spec (str, tuple, list, or np.ndarray): Axis specification.
            theta (float): Angle in radians.
            save_rotated_tensors (bool): If True, save rotated tensors to CSV (default False).
        Returns:
            Material: The rotated material object.
        """

        rotation_axis = voigt.parse_rotation_axis(rot_axis_spec)
        matR = voigt.make_rotation_matrix(rotation_axis, theta)

        self.total_matR = matR @ self.total_matR

        self.stiffness_c_IJ.rotate(matR)
        self.photoel_p_IJ.rotate(matR)
        self.viscosity_eta_IJ.rotate(matR)
        self.optdiel_eps_ij.rotate(matR)

        if self._piezo is not None:
            self._piezo.rotate(matR)

        self.stiffness_c_IJ.check_symmetries()

        caxes = self._crystal_axes.copy()
        self.set_crystal_axes(
            voigt.rotate_3vector(caxes[0], matR),
            voigt.rotate_3vector(caxes[1], matR),
            voigt.rotate_3vector(caxes[2], matR),
        )

    # restore orientation to original axes in spec file.
    def reset_orientation(self):
        """
        Restore the material and all tensors to their original orientation.
        """

        self.total_matR = np.eye(3)

        self.stiffness_c_IJ = copy.deepcopy(self._stiffness_c_IJ_orig)
        self.photoel_p_IJ = copy.deepcopy(self._photoel_p_IJ_orig)
        self.viscosity_eta_IJ = copy.deepcopy(self._viscosity_eta_IJ_orig)
        self.optdiel_eps_ij = copy.deepcopy(self._optdiel_eps_ij_orig)

        if self._piezo:
            self._piezo.reset_orientation()

        self.set_crystal_axes(unit_x, unit_y, unit_z)

    # rotate original crystal to specific named-orientation, eg x-cut, y-cut. '111' etc.
    def set_orientation(self, label):
        """
        Set the orientation of the crystal to a specific named orientation (e.g., 'x-cut', '111').

        Args:
            label (str): Orientation label defined in the material file.
        """

        self.reset_orientation()

        try:
            ocode = self._d_params[f"orientation_{label.lower()}"]
        except KeyError:
            reporting.report_and_exit(
                f'Orientation "{label}" is not defined for material {self.material_name}.'
            )

        if ocode == "ident":  # native orientation is the desired one
            return

        try:
            ux, uy, uz, rot = map(float, ocode.split(","))
        except Exception:
            reporting.report_and_exit(
                f"Can't parse crystal orientation code {ocode} for material {self.material_name}."
            )

        rot_axis = np.array((ux, uy, uz))
        theta = rot * np.pi / 180

        self.rotate(rot_axis, theta)

    def set_crystal_axes(self, va, vb, vc):
        """
        Set the crystal axes vectors.

        Args:
            va, vb, vc (np.ndarray): Crystal axis vectors.
        """

        self._crystal_axes = [va, vb, vc]

    def construct_crystal_isotropic(self):
        """
        Construct the tensors for an isotropic crystal from the parameter dictionary.
        """

        # ordinary Cartesian axes for the crystal axes
        self.set_crystal_axes(unit_x, unit_y, unit_z)

        self._anisotropic = False



        # Try to read isotropic from stiffness and then from Young's modulus and Poisson ratio
        if "c_11" in self._d_params and "c_12" in self._d_params and "c_44" in self._d_params:

            self.stiffness_c_IJ.load_isotropic_from_json(self._d_params)
            mu = self.stiffness_c_IJ.mat[4, 4]
            lam = self.stiffness_c_IJ.mat[1, 2]
            r = lam / mu
            self.nuPoisson = 0.5 * r / (1 + r)
            self.EYoung = 2 * mu * (1 + self.nuPoisson)
            self.Lame_mu = mu
            self.Lame_lambda = lam

        elif "EYoung" in self._d_params and "nuPoisson" in self._d_params:
            self.EYoung = self._d_params["EYoung"]
            self.nuPoisson = self._d_params["nuPoisson"]
            c44 = 0.5 * self.EYoung / (1 + self.nuPoisson)
            c12 = (
                self.EYoung
                * self.nuPoisson
                / ((1 + self.nuPoisson) * (1 - 2 * self.nuPoisson))
            )
            c11 = c12 + 2 * c44
            self.stiffness_c_IJ = voigt.VoigtTensor4_IJ(
                "c", "Stiffness", "c_IJ", unit=("GPa", 1.0e9)
            )
            self.stiffness_c_IJ.make_isotropic_tensor(c11, c12, c44)

            self.Lame_mu = self.stiffness_c_IJ.mat[4, 4]
            self.Lame_lambda = self.stiffness_c_IJ.mat[1, 2]
        else:
            reporting.report_and_exit(
                "Broken isotropic material file:" + self.json_file
            )

        self.photoel_p_IJ.load_isotropic_from_json(self._d_params)
        self.viscosity_eta_IJ.load_isotropic_from_json(self._d_params)

        self.stiffness_c_IJ.check_symmetries()

    def construct_crystal(self):
        """
        Construct the tensors for the crystal based on its symmetry class.
        """

        if self.crystal_class == CrystalGroup.Isotropic:
            self.construct_crystal_isotropic()
            return

        self._anisotropic = True

        # TODO: change to match/case
        if self.crystal_class == CrystalGroup.GeneralAnisotropic:
            self.construct_crystal_general()
        else:
            crystal_axes, C_indep, C_dep, pIJ_indep, pIJ_dep = construct_crystal_for_symmetry_class(
                self.crystal_class)

            self.set_crystal_axes(*crystal_axes)

            self.stiffness_c_IJ.set_from_structure_elts(
                self._d_params, C_indep, C_dep
            )
            self.viscosity_eta_IJ.set_from_structure_elts(
                self._d_params, C_indep, C_dep
            )
            self.photoel_p_IJ.set_from_structure_elts(
                self._d_params, pIJ_indep, pIJ_dep
            )

        self.stiffness_c_IJ.check_symmetries()

    def plot_bulk_dispersion_ivp(
        self, pref, label=None, show_poln=True, flip_x=False, flip_y=False,
        cut_plane="xz", mark_velocities=() ):
        """
        Plot the bulk dispersion as inverse phase velocity. Returns the filename of the generated image.
        """
        return plms.plot_bulk_dispersion_ivp(
            self, pref, label, show_poln, flip_x, flip_y, cut_plane, mark_velocities
        )


    def plot_bulk_dispersion_vg(
        self, pref, label=None, show_poln=True, flip_x=False, flip_y=False,
        cut_plane="xz", mark_velocities=()
    ):
        """
        Plot the bulk dispersion as the group velocity (ray surface). Returns the filename of the generated image.
        """

        return plms.plot_bulk_dispersion_vg(
            self, pref, label, show_poln, flip_x, flip_y, cut_plane, mark_velocities
        )



    def plot_bulk_dispersion_all(
        self, pref, label=None, show_poln=True, flip_x=False, flip_y=False,
        cut_plane="xz", mark_velocities=()
    ):
        """Plot the bulk dispersion and ray surfaces in the x-z plane for the current orientation (all modes).
        Returns the filename of the generated image.

        Solving the Christoffel equation: D C D^T u = -\rho v_p^2 u, for eigenvalue v_p
        and eigengector u.

        C is the Voigt form stiffness.
        D = [
        [kapx  0   0   0  kapz  kapy  ]
        [0   kapy  0   kapz 0   kapx  ]
        [0   0   kapz  kapy kapx  0]] where kap=(cos phi, 0, sin phi).

        Returns filename of the generated image.
        """

        return plms.plot_bulk_dispersion_2D_all(
            self, pref, label, show_poln, flip_x, flip_y, cut_plane, mark_velocities
        )

    def plot_bulk_dispersion_3D(self, pref):
        """
        Generate isocontour surfaces of the bulk dispersion in 3D k-space.
        Returns the filename of the generated image.
        """

        return plms.plot_bulk_dispersion_3D(self, pref)

    def make_crystal_axes_plot(self, pref):
        """
        Build a crystal coordinates diagram using an external application. Returns the output PNG filename.
        """

        return plms.make_crystal_axes_plot(pref, self._crystal_axes)

    def plot_photoelastic_IJ(self, prefix, v_comps):
        """
        Plot photoelastic tensor components as a function of rotation angle about the y-axis.

        Args:
            prefix (str): Prefix for output file names.
            v_comps (list): List of element strings (e.g., '11', '12').
        """

        return plms.plot_material_photoelastic_IJ(prefix, v_comps, self)

    def get_stiffness_for_kappa(self, vkap):
        """
        Get the stiffness tensor, piezo-hardened if piezoelectric effects are active.

        Args:
            vkap (np.ndarray): Propagation direction (3-vector).
        Returns:
            VoigtTensor4_IJ: Stiffened stiffness tensor.
        """

        if self._piezo is None or not self._piezo.active:
            return self.stiffness_c_IJ
        else:
            return self._piezo.make_piezo_stiffened_stiffness_for_kappa(vkap)

    # deprecated
    def rotate_axis(self, rotation_axis, theta, save_rotated_tensors=False):
        """
        Deprecated. Use rotate().
        """
        reporting.register_warning("rotate_axis function is deprecated. Use rotate()")
        self.rotate(rotation_axis, theta, save_rotated_tensors)


def compare_bulk_dispersion(mat1, mat2, pref):
    return plms.compare_bulk_dispersion(mat1, mat2, pref)


def isotropic_stiffness(E, v):
    """
    Calculate the stiffness matrix components of isotropic
    materials, given the two free parameters.

    Ref: www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm

    Args:
        E  (float): Youngs modulus

        v  (float): Poisson ratio
    """
    c_11 = E * (1 - v) / ((1 + v) * (1 - 2 * v))
    c_12 = E * (v) / ((1 + v) * (1 - 2 * v))
    c_44 = (E * (1 - 2 * v) / ((1 + v) * (1 - 2 * v))) / 2

    return c_11, c_12, c_44


def disprel_rayleigh(vR, vl):

    return vR**6 - 8 * vR**4 + vR**2 * (24 - 16 / vl**2) + 16 * (1 / vl**2 - 1)


def find_Rayleigh_velocity(Vs, Vl):
    """Find Rayleigh velocity v_q for isotropic mode with velocities Vs and Vl."""

    vl = Vl / Vs

    vRlo = 0.001
    vRhi = 1

    def dr_rayleigh(vR):
        return disprel_rayleigh(vR, vl)

    vres = sciopt.root_scalar(dr_rayleigh, bracket=(vRlo, vRhi))

    if not vres.converged:
        raise ValueError(vres.flag)
    else:
        vR = vres.root

    return vR * Vs
