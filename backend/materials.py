# materials.py is a subroutine of NumBAT that defines Material objects,
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


import copy
import json
import math
import pathlib
import re
import subprocess
import tempfile
import traceback

import matplotlib as mpl
import matplotlib.cm as mplcm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from nbtypes import CrystalGroup, unit_x, unit_y, unit_z, SI_permittivity_eps0
import reporting
import voigt
from bulkprops import solve_christoffel
import numbattools as nbtools

import plotting.materials as plms

class BadMaterialFileError(Exception):
    pass


g_material_library = None


def make_material(s):
    global g_material_library

    if g_material_library is None:
        g_material_library = MaterialLibrary()

    return g_material_library.get_material(s)


class MaterialLibrary:
    def __init__(self):
        self._material_data_path = ""
        self._materials = {}

        # identify mat data directory:  backend/material_data
        this_dir = pathlib.Path(__file__).resolve().parent
        self._material_data_path = this_dir / "material_data"

        self._load_materials()

    def get_material(self, matname):
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
                reporting.report_and_exit(f"Material file {fname} has no 'material_name' field.")

            if mat_name in self._materials:
                reporting.report_and_exit(f"Material file {fname} has the same name as an existing material {mat_name}.")

            try:
                new_mat = Material(json_data, fname)
            except BadMaterialFileError as err:
                reporting.report_and_exit(str(err))

            self._materials[mat_name] = new_mat


class PiezoElectricProperties:
    def __init__(self, d_piezo, stiffness_cE_IJ, diel_epsT_ij):
        # cE denotes stiffness tensor at fixed electric field E
        # epsT denotes dielecric tensor at fixed stress T

        self.active =  'piezo_active' in d_piezo

        #
        self._tens_cE_IJ = stiffness_cE_IJ

        # piezoelectric strain quantities
        self._tens_d_iJ = None      # strain piezo tensor in Voigt form
        #self._tens_dbar_Ij = None   # transpose \underscore{d} in Voigt form
        self._diel_epsT_ij = None   # zero stress dielectric tensor

        # piezoelectric stress quantities
        self._tens_e_iJ = None      # stress piezo tensor in Voigt form       = d_iJ cE_JI
        #self._tens_ebar_Ij = None   # transpose \underscore{d} in Voigt form

        self._diel_epsS_ij = None   # zero strain dielectric tensor,          = epsT_ij - d_iJ cE_JI dbar_Ij

        # Units of [Coulomb/Newton], order 1e-12
        self._tens_d_iJ = voigt.VoigtTensor3_iJ('piezo_d', 'Strain piezo coefficient', 'd_iJ', ('pC/N', 1e-12), transforms_with_factor_2=True)


        self._tens_d_iJ.set_from_json(d_piezo)

        # Units of [Coulomb/m^2], order 1
        self._tens_e_iJ = voigt.VoigtTensor3_iJ('piezo_e', 'Stress piezo coefficient', 'e_iJ', ('C/m^2', 1.0))
        self._tens_e_iJ.set_from_00matrix(self._tens_d_iJ.value() @ stiffness_cE_IJ.value() )  # Auld 8.40

        self._make_eps_S(stiffness_cE_IJ, diel_epsT_ij)

        self.store_original_tensors()

        v_kap = np.array([4,3,1], dtype=np.float64)
        v_kap /= np.sqrt(np.abs(v_kap @ v_kap.T))

        self.make_piezo_stiffened_stiffness(v_kap)

    def make_piezo_stiffened_stiffness(self, v_kap):
        '''Returns piezo-stiffened stiffness (!) for propagation along direction v_kappa. See Auld 8.147'''

        e_iL = self._tens_e_iJ.value()      # [0..2, 0..5]
        e_Kj = self._tens_e_iJ.value().T    # [0..5, 0..2]

        e_Kj_lj = e_Kj @ v_kap
        li_e_iL = v_kap @ e_iL

        l_eps_l = v_kap @ self._diel_epsS_ij @ v_kap * SI_permittivity_eps0

        print('shapes 1', e_Kj.shape, v_kap.shape, e_Kj_lj.shape)
        print('shapes 2', e_iL.shape, v_kap.shape, li_e_iL.shape)

        print('e_Kj_lj', e_Kj_lj)
        print('li_e_iL', li_e_iL)
        print('l_eps_', l_eps_l)

        cE_KL_stiff = self._tens_cE_IJ.value() + np.outer(e_Kj_lj, li_e_iL) / l_eps_l

        print('cE_stiff', cE_KL_stiff)
        print('cE_diff', cE_KL_stiff-self._tens_cE_IJ.value())

        # convert to a voigt4
        return  cE_KL_stiff



    def _make_eps_S(self, cE_IJ, epsT_ij):
        d_iJ = self._tens_d_iJ
      #  mat_dbar_Ij = mat_d_iJ.T

        self._diel_epsS_ij = epsT_ij - d_iJ.value() @ cE_IJ.value() @ d_iJ.value().T   # Auld 8.41

    def rotate(self, matR):
        '''Rotate piezo matrices by SO(3) matrix matR'''

        self._tens_d_iJ.rotate(matR)
        self._tens_e_iJ.rotate(matR)

        self._diel_epsS_ij = voigt._rotate_2tensor(self._diel_epsS_ij, matR)

    def store_original_tensors(self):
        '''Save tensors to optionally restore after rotations.'''
        self._tens_d_iJ_orig = copy.deepcopy(self._tens_d_iJ)
        #self._tens_dbar_Ij_orig = self._tens_dbar_Ij.copy()
        self._diel_epsS_iJjorig = copy.deepcopy(self._diel_epsS_ij)


    def d_iJ_to_d_ijk(self):
        pass

    def d_ijk_to_d_iJ(self):
        pass

    def __str__(self):
        s='Piezoelectric properties:'
        if self.active:
            s+= '\n Piezo effects enabled'
        else:
            s+= '\n Piezo effects disabled'

        s += '\n' + str(self._tens_d_iJ)
        #s+=str(self._diel_epsT_iJ)
        s += '\n' + str(self._tens_e_iJ)

        s += '\n\n Dielectric tensor epsS_ij\n' + nbtools.indent_string(str(self._diel_epsS_ij), 4)

        return s




class Material(object):
    """Class representing a waveguide material.

    This should not be constructed directly but by calling materials.get_material()

      Materials include the following properties and corresponding units:
          -  Refractive index []
          -  Density [kg/m3]
          -  Stiffness tensor component [Pa]
          -  Photoelastic tensor component []
          -  Acoustic loss tensor component [Pa s]

    """

    def __init__(self, json_data, filename):
        # a,b,c crystal axes according to standard conventions
        self._crystal_axes = []

        self._anisotropic = False
        self.nuPoisson = 0.0
        self.EYoung = 0.0

        self.refindex_n = 0

        # current rotated tensors
        self.stiffness_c_IJ = None
        self.viscosity_eta_IJ = None
        self.photoel_p_IJ = None
        self.eps_diel_ij = None

        # home crystal orientation tensors
        self._stiffness_c_IJ_orig = None
        self._photoel_p_IJ_orig = None
        self._viscosity_eta_IJ_orig = None
        self._eps_diel_ij_orig = None

        self._piezo = None

        self._parse_json_data(json_data, filename)

    def __str__(self):
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
        return self._piezo is not None

    def copy(self):
        return copy.deepcopy(self)

    def full_str(self):
        s = str(self)
        s += '\n\nStiffness:' + str(self.stiffness_c_IJ)
        s += '\n\nViscocity:' + str(self.viscosity_eta_IJ)
        s += '\n\nPhotoelasticity:' + str(self.photoel_p_IJ)

        if self._piezo:
            s+= '\n\n' + str(self._piezo)
        return s

    def elastic_properties(self):
        """Returns a string containing key elastic properties of the material."""

        dent = "\n  "
        try:
            s = f"Elastic properties of material {self.material_name}"
            s += dent + f"Density:        {self.rho:.3f} kg/m^3"
            s += dent + f"Ref. index:     {self.refindex_n:.4f} "

            s += dent + f"Crystal class:  {self.crystal.name}"

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
            else:
                s += dent + "Stiffness c_IJ:" + str(self.stiffness_c_IJ) + "\n"

                # find wave properties for z propagation
                v_phase, v_evecs, v_vgroup = solve_christoffel( unit_z, self.stiffness_c_IJ, self.rho)

                with np.printoptions(
                    precision=4, floatmode="fixed", sign=" ", suppress=True
                ):
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

            if self._piezo:
                s += str (self._piezo)

        except Exception:
            s = "Unknown/undefined elastic parameters in material " + self.material_name
        return s

    def Vac_longitudinal(self):
        """For an isotropic material, returns the longitudinal (P-wave) elastic phase velocity."""
        assert self.is_isotropic()

        if not self.rho or self.rho == 0:  # Catch vacuum cases
            return 0.0
        else:
            return math.sqrt(self.stiffness_c_IJ[1, 1] / self.rho)

    def Vac_shear(self):
        """For an isotropic material, returns the shear (S-wave) elastic phase velocity."""
        assert self.is_isotropic()

        if not self.rho or self.rho == 0:  # Catch vacuum cases
            return 0.0
        else:
            return math.sqrt(self.stiffness_c_IJ[4, 4] / self.rho)

    def Vac_phase(self):
        '''Returns triple of phase velocities for propagation along z given current orientation of cyrstal.'''
        vkap=np.array([0,0,1], dtype=np.float64)
        v_vphase, vecs, v_vgroup = solve_christoffel(vkap, self.stiffness_c_IJ, self.rho)
        return v_vphase


    def has_elastic_properties(self):
        """Returns true if the material has at least some elastic properties defined."""
        return self.rho is not None

    def _parse_json_data(self, json_data, fname):
        """
        Load material data from json file.

        Args:
            data_file  (str): name of data file located in NumBAT/backend/material_data

        """
        self._params = json_data  # Do without this?
        self.json_file = fname

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

        Re_n = json_data["Re_n"]  # Real part of refractive index []
        # Imaginary part of refractive index []
        Im_n = json_data["Im_n"]
        self.refindex_n = Re_n + 1j * Im_n  # Complex refractive index []
        self.eps_diel_ij = np.eye(3) * self.refindex_n**2

        self.rho = json_data["s"]  # Density [kg/m3]

        if self.is_vacuum():  # no mechanical properties available
            return

        self.EYoung = None
        self.nuPoisson = None

        if "crystal_class" not in json_data:
            raise BadMaterialFileError(
                f"Material file {fname} has no 'crystal_class' field."
            )
        try:
            self.crystal = CrystalGroup[json_data["crystal_class"]]
        except ValueError as exc:
            print("Unknown crystal class in material data file")
            raise BadMaterialFileError(
                f"Unknown crystal class in material data file {fname}"
            ) from exc

        if self.crystal == CrystalGroup.Isotropic:
            self.construct_crystal_isotropic()
        else:
            self.construct_crystal_anisotropic()


        dpiezo = {k:v for k,v in json_data.items() if k.startswith('piezo')}
        if dpiezo:
            self._piezo = PiezoElectricProperties(dpiezo, self.stiffness_c_IJ, self.eps_diel_ij)

        self._store_original_tensors()

    def _store_original_tensors(self):
        self._stiffness_c_IJ_orig = self.stiffness_c_IJ
        self._photoel_p_IJ_orig = self.photoel_p_IJ
        self._viscosity_eta_IJ_orig = self.viscosity_eta_IJ
        self._eps_diel_ij_orig = self.eps_diel_ij

    def is_vacuum(self):
        """Returns True if the material is the vacuum."""
        return self.chemical == "Vacuum"

    # (don't really need this as isotropic materials are the same)
    def construct_crystal_cubic(self):
        # plain cartesian axes
        self.set_crystal_axes(unit_x, unit_y, unit_z)

        try:
            self.stiffness_c_IJ.read_from_json(1, 1)
            self.stiffness_c_IJ.read_from_json(1, 2)
            self.stiffness_c_IJ[1, 3] = self.stiffness_c_IJ[1, 2]
            self.stiffness_c_IJ[2, 1] = self.stiffness_c_IJ[1, 2]
            self.stiffness_c_IJ[2, 2] = self.stiffness_c_IJ[1, 1]
            self.stiffness_c_IJ[2, 3] = self.stiffness_c_IJ[1, 2]
            self.stiffness_c_IJ[3, 1] = self.stiffness_c_IJ[1, 2]
            self.stiffness_c_IJ[3, 2] = self.stiffness_c_IJ[1, 2]
            self.stiffness_c_IJ[3, 3] = self.stiffness_c_IJ[1, 1]
            self.stiffness_c_IJ.read_from_json(4, 4)
            self.stiffness_c_IJ[5, 5] = self.stiffness_c_IJ[4, 4]
            self.stiffness_c_IJ[6, 6] = self.stiffness_c_IJ[4, 4]

            self.viscosity_eta_IJ.read_from_json(1, 1)
            self.viscosity_eta_IJ.read_from_json(1, 2)
            self.viscosity_eta_IJ[1, 3] = self.viscosity_eta_IJ[1, 2]
            self.viscosity_eta_IJ[2, 1] = self.viscosity_eta_IJ[1, 2]
            self.viscosity_eta_IJ[2, 2] = self.viscosity_eta_IJ[1, 1]
            self.viscosity_eta_IJ[2, 3] = self.viscosity_eta_IJ[1, 2]
            self.viscosity_eta_IJ[3, 1] = self.viscosity_eta_IJ[1, 2]
            self.viscosity_eta_IJ[3, 2] = self.viscosity_eta_IJ[1, 2]
            self.viscosity_eta_IJ[3, 3] = self.viscosity_eta_IJ[1, 1]
            self.viscosity_eta_IJ.read_from_json(4, 4)
            self.viscosity_eta_IJ[5, 5] = self.viscosity_eta_IJ[4, 4]
            self.viscosity_eta_IJ[6, 6] = self.viscosity_eta_IJ[4, 4]

            self.photoel_p_IJ.read_from_json(1, 1)
            self.photoel_p_IJ.read_from_json(1, 2)

            self.photoel_p_IJ[1, 3] = self.photoel_p_IJ[1, 2]
            self.photoel_p_IJ[2, 1] = self.photoel_p_IJ[1, 2]
            self.photoel_p_IJ[2, 2] = self.photoel_p_IJ[1, 1]
            self.photoel_p_IJ[2, 3] = self.photoel_p_IJ[1, 2]
            self.photoel_p_IJ[3, 1] = self.photoel_p_IJ[1, 2]
            self.photoel_p_IJ[3, 2] = self.photoel_p_IJ[1, 2]
            self.photoel_p_IJ[3, 3] = self.photoel_p_IJ[1, 1]
            self.photoel_p_IJ.read_from_json(4, 4)

            # According to Powell, for Oh group, these are distinct elements, but no one seems to quote them
            if not self.photoel_p_IJ.read_from_json(5, 5, optional=True):
                self.photoel_p_IJ[5, 5] = self.photoel_p_IJ[4, 4]
            if not self.photoel_p_IJ.read_from_json(6, 6, optional=True):
                self.photoel_p_IJ[6, 6] = self.photoel_p_IJ[4, 4]

        except Exception:
            reporting.report_and_exit(
                f"Failed to load cubic crystal class in material data file {self.json_file}"
            )

    def construct_crystal_trigonal(self):
        # Good source for these rules is the supp info of doi:10.1364/JOSAB.482656 (Gustavo surface paper)

        self.set_crystal_axes(unit_x, unit_y, unit_z)

        try:
            for lintens in [self.stiffness_c_IJ, self.viscosity_eta_IJ]:
                for i, j in [(1, 1), (1, 2), (1, 3), (1, 4), (3, 3), (4, 4)]:
                    lintens.read_from_json(i, j)

                lintens[2, 1] = lintens[1, 2]
                lintens[2, 2] = lintens[1, 1]
                lintens[2, 3] = lintens[1, 3]
                lintens[2, 4] = -lintens[1, 4]

                lintens[3, 1] = lintens[1, 3]
                lintens[3, 2] = lintens[1, 3]

                lintens[4, 1] = lintens[1, 4]
                lintens[4, 2] = -lintens[1, 4]

                lintens[5, 5] = lintens[4, 4]
                lintens[5, 6] = lintens[1, 4]
                lintens[6, 5] = lintens[1, 4]
                lintens[6, 6] = (lintens[1, 1] - lintens[1, 2]) / 2.0

            # TODO: confirm correct symmetry properties for p.
            # PreviouslyuUsing trigonal = C3v from Powell, now the paper above
            self.photoel_p_IJ.read_from_json(1, 1)
            self.photoel_p_IJ.read_from_json(1, 2)
            self.photoel_p_IJ.read_from_json(1, 3)
            self.photoel_p_IJ.read_from_json(1, 4)
            self.photoel_p_IJ.read_from_json(3, 1)
            self.photoel_p_IJ.read_from_json(3, 3)
            self.photoel_p_IJ.read_from_json(4, 1)
            self.photoel_p_IJ.read_from_json(4, 4)

            self.photoel_p_IJ[2, 1] = self.photoel_p_IJ[1, 2]
            self.photoel_p_IJ[2, 2] = self.photoel_p_IJ[1, 1]
            self.photoel_p_IJ[2, 3] = self.photoel_p_IJ[1, 3]
            self.photoel_p_IJ[2, 4] = -self.photoel_p_IJ[1, 4]

            self.photoel_p_IJ[3, 2] = self.photoel_p_IJ[3, 1]

            self.photoel_p_IJ[4, 2] = -self.photoel_p_IJ[4, 1]

            self.photoel_p_IJ[5, 5] = self.photoel_p_IJ[4, 4]
            self.photoel_p_IJ[5, 6] = self.photoel_p_IJ[4, 1]
            self.photoel_p_IJ[6, 5] = self.photoel_p_IJ[1, 4]
            self.photoel_p_IJ[6, 6] = (
                self.photoel_p_IJ[1, 1] - self.photoel_p_IJ[1, 2]
            ) / 2

        except Exception:
            reporting.report_and_exit(
                f"Failed to load trigonal crystal class in material data file {self.json_file}"
            )

    def construct_crystal_general(self):
        try:  # full anisotropic tensor components
            for i in range(1, 7):
                for j in range(1, 7):
                    self.stiffness_c_IJ.read_from_json(i, j)
                    self.photoel_p_IJ.read_from_json(i, j)
                    self.viscosity_eta_IJ.read_from_json(i, j)

        except KeyError:
            reporting.report_and_exit(
                "Failed to load anisotropic crystal class in material data file {self.json_file}"
            )

    # TODO: make a property
    def set_refractive_index(self, nr, ni=0.0):
        self.refindex_n = nr + 1j * ni

    def is_isotropic(self):
        return not self._anisotropic

    # deprecated
    def rotate_axis(self, rotation_axis, theta, save_rotated_tensors=False):
        reporting.register_warning("rotate_axis function is depprecated. Use rotate()")
        self.rotate(rotation_axis, theta, save_rotated_tensors)

    def rotate(self, rot_axis_spec, theta, save_rotated_tensors=False):
        """Rotate crystal axis by theta radians.

        Args:
            theta  (float): Angle to rotate by in radians.

            rotate_axis  (str): Axis around which to rotate.

        Keyword Args:
            save_rotated_tensors  (bool): Save rotated tensors to csv.

        Returns:
            ``Material`` object with rotated tensor values.
        """

        rotation_axis = voigt.parse_rotation_axis(rot_axis_spec)
        matR = voigt.make_rotation_matrix(rotation_axis, theta)

        self.stiffness_c_IJ.rotate(matR)
        self.photoel_p_IJ.rotate(matR)
        self.viscosity_eta_IJ.rotate(matR)

        self._piezo.rotate(matR)

        self.stiffness_c_IJ.check_symmetries()

        caxes = self._crystal_axes.copy()
        self.set_crystal_axes(
            voigt.rotate_3vector(caxes[0], matR),
            voigt.rotate_3vector(caxes[1], matR),
            voigt.rotate_3vector(caxes[2], matR),
        )

        if save_rotated_tensors:
            np.savetxt(
                "rotated_stiffness_c_IJ.csv", self.stiffness_c_IJ.mat, delimiter=","
            )
            np.savetxt("rotated_photoel_p_IJ.csv", self.photoel_p_IJ.mat, delimiter=",")
            np.savetxt(
                "rotated_viscosity_eta_IJ.csv", self.viscosity_eta_IJ.mat, delimiter=","
            )

    # restore orientation to original axes in spec file.
    def reset_orientation(self):
        self.stiffness_c_IJ = copy.deepcopy(self._stiffness_c_IJ_orig)
        self.photoel_p_IJ = copy.deepcopy(self._photoel_p_IJ_orig)
        self.viscosity_eta_IJ = copy.deepcopy(self._viscosity_eta_IJ_orig)

        self.set_crystal_axes(unit_x, unit_y, unit_z)

    # rotate original crystal to specific named-orientation, eg x-cut, y-cut. '111' etc.
    def set_orientation(self, label):
        self.reset_orientation()

        try:
            ocode = self._params[f"orientation_{label.lower()}"]
        except KeyError:
            reporting.report_and_exit( f'Orientation "{label}" is not defined for material {self.material_name}.')

        if ocode == "ident":  # native orientation is the desired one
            return

        try:
            ux, uy, uz, rot = map(float, ocode.split(","))
        except Exception:
            reporting.report_and_exit( f"Can't parse crystal orientation code {ocode} for material {self.material_name}.")

        rot_axis = np.array((ux, uy, uz))
        theta = rot * np.pi / 180

        self.rotate(rot_axis, theta)

    def set_crystal_axes(self, va, vb, vc):
        self._crystal_axes = [va, vb, vc]

    def construct_crystal_isotropic(self):
        # ordinary Cartesian axes for the crystal axes
        self.set_crystal_axes(unit_x, unit_y, unit_z)

        self._anisotropic = False

        self.stiffness_c_IJ = voigt.VoigtTensor4( self.material_name, "c", self._params, "stiffness", ("GPa", 1.0e9))
        self.viscosity_eta_IJ = voigt.VoigtTensor4( self.material_name, "eta", self._params, "viscosity")
        self.photoel_p_IJ = voigt.VoigtTensor4( self.material_name, "p", self._params, "photoelasticity")

        # Try to read isotropic from stiffness and then from Young's modulus and Poisson ratio
        if "c_11" in self._params and "c_12" in self._params and "c_44" in self._params:
            #    self.stiffness_c_IJ = VoigtTensor4(self.material_name, 'c', self._params)
            self.stiffness_c_IJ.load_isotropic_from_json()
            mu = self.stiffness_c_IJ.mat[4, 4]
            lam = self.stiffness_c_IJ.mat[1, 2]
            r = lam / mu
            self.nuPoisson = 0.5 * r / (1 + r)
            self.EYoung = 2 * mu * (1 + self.nuPoisson)
            self.Lame_mu = mu
            self.Lame_lambda = lam

        elif "EYoung" in self._params and "nuPoisson" in self._params:
            self.EYoung = self._params["EYoung"]
            self.nuPoisson = self._params["nuPoisson"]
            c44 = 0.5 * self.EYoung / (1 + self.nuPoisson)
            c12 = ( self.EYoung * self.nuPoisson / ((1 + self.nuPoisson) * (1 - 2 * self.nuPoisson)))
            c11 = c12 + 2 * c44
            self.stiffness_c_IJ = voigt.VoigtTensor4(self.material_name, "c")
            self.stiffness_c_IJ.make_isotropic_tensor(c11, c12, c44)

            self.Lame_mu = self.stiffness_c_IJ.mat[4, 4]
            self.Lame_lambda = self.stiffness_c_IJ.mat[1, 2]
        else:
            reporting.report_and_exit( "Broken isotropic material file:" + self.json_file)

        self.photoel_p_IJ.load_isotropic_from_json()
        self.viscosity_eta_IJ.load_isotropic_from_json()

        self.stiffness_c_IJ.check_symmetries()

    # not do this unless symmetry is off?
    def construct_crystal_anisotropic(self):
        self.stiffness_c_IJ = voigt.VoigtTensor4( self.material_name, "c", self._params, "stiffness", ("GPa", 1.0e9))
        self.viscosity_eta_IJ = voigt.VoigtTensor4( self.material_name, "eta", self._params, "viscosity")
        self.photoel_p_IJ = voigt.VoigtTensor4( self.material_name, "p", self._params, "photoelasticity")

        self._anisotropic = True

        # TODO: change to match/case
        if self.crystal == CrystalGroup.Trigonal:
            self.construct_crystal_trigonal()
        elif self.crystal == CrystalGroup.Cubic:
            self.construct_crystal_cubic()
        elif self.crystal == CrystalGroup.GeneralAnisotropic:
            self.construct_crystal_general()

        self.stiffness_c_IJ.check_symmetries()


    def plot_bulk_dispersion(self, pref, label=None, show_poln=True):
        """Draw slowness surface 1/v_p(kappa) and ray surface contours in the horizontal
        (x-z) plane for the crystal axes current orientation.

        Solving the Christoffel equation: D C D^T u = -\rho v_p^2 u, for eigenvalue v_p
        and eigengector u.

        C is the Voigt form stiffness.
        D = [
        [kapx  0   0   0  kapz  kapy  ]
        [0   kapy  0   kapz 0   kapx  ]
        [0   0   kapz  kapy kapx  0]] where kap=(cos phi, 0, sin phi).

        """

        plms.plot_bulk_dispersion_2D(self, pref, label, show_poln)

    def plot_bulk_dispersion_3D(self, pref):
        """
        Generate isocontour surfaces of the bulk dispersion in 3D k-space.
        """

        plms.plot_bulk_dispersion_3D(self, pref)


    def make_crystal_axes_plot(self, pref):
        """Build crystal coordinates diagram using call to external asymptote application."""

        plms.make_crystal_axes_plot(pref)

    def plot_photoelastic_IJ(self, prefix, v_comps):
        """ Plot photoelastic tensor components as a function of rotation angle about y-axis.

        Args:
            prefix (str): Prefix for output file names.
            v_comps (list): List of strings of desired elements: "11", "12", "31" etc

        """

        plms.plot_material_photoelastic_IJ(prefix, v_comps, self)


def compare_bulk_dispersion(mat1, mat2, pref):
    plms.compare_bulk_dispersion(mat1, mat2, pref)


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




