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

from nbtypes import CrystalGroup, unit_x, unit_y, unit_z
import reporting
import voigt
from bulkprops import solve_christoffel


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

        self.stiffness_c_IJ = None
        self.viscosity_eta_IJ = None
        self.photoel_p_IJ = None

        self.refindex_n = 0

        self._stiffness_c_IJ_orig = None
        self._photoel_p_IJ_orig = None
        self._viscosity_eta_IJ_orig = None

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

    def copy(self):
        return copy.deepcopy(self)

    def full_str(self):
        s = str(self)
        s += str(self.stiffness_c_IJ)
        s += str(self.viscosity_eta_IJ)
        s += str(self.photoel_p_IJ)
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

        self._store_original_tensors()

    def _store_original_tensors(self):
        self._stiffness_c_IJ_orig = self.stiffness_c_IJ
        self._photoel_p_IJ_orig = self.photoel_p_IJ
        self._viscosity_eta_IJ_orig = self.viscosity_eta_IJ

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

        elif "EYoung" in self._params and "nuPoisson" in self._params:
            self.EYoung = self._params["EYoung"]
            self.nuPoisson = self._params["nuPoisson"]
            c44 = 0.5 * self.EYoung / (1 + self.nuPoisson)
            c12 = ( self.EYoung * self.nuPoisson / ((1 + self.nuPoisson) * (1 - 2 * self.nuPoisson)))
            c11 = c12 + 2 * c44
            self.stiffness_c_IJ = voigt.VoigtTensor4(self.material_name, "c")
            self.stiffness_c_IJ.make_isotropic_tensor(c11, c12, c44)
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






    def _add_3d_dispersion_curves_to_axes(self, ax_ivp=None, ax_vg=None):
        '''
        Draw phase and group velocity surfaces on 3D axes.

        :param ax_ivp: 3D axes for the slowness (reciprocal phase velocity).
        :param ax_vg: 3D axes for the group velocity.
        '''
        axs = []

        if ax_ivp is not None: axs.append(ax_ivp)
        if ax_vg is not None: axs.append(ax_vg)

        # Make data
        tpts = 50
        ppts = 100
        vphi = np.linspace(0, 2 * np.pi, ppts)
        vtheta = np.linspace(0, np.pi, tpts)

        ivx = np.zeros([tpts, ppts, 3])
        ivy = np.zeros([tpts, ppts, 3])
        ivz = np.zeros([tpts, ppts, 3])

        ivgx = np.zeros([tpts, ppts, 3])
        ivgy = np.zeros([tpts, ppts, 3])
        ivgz = np.zeros([tpts, ppts, 3])

        for ip, phi in enumerate(vphi):
            for itheta, theta in enumerate(vtheta):
                vkap = np.array(
                    [
                        np.sin(theta) * np.cos(phi),
                        np.sin(theta) * np.sin(phi),
                        np.cos(theta),
                    ]
                )
                v_vphase, vecs, v_vgroup = solve_christoffel( vkap, self.stiffness_c_IJ, self.rho)

                # slowness curve  eta(vkap) = 1/v_phase(vkap)
                ivx[itheta, ip, :] = vkap[0] / v_vphase
                ivy[itheta, ip, :] = vkap[1] / v_vphase
                ivz[itheta, ip, :] = vkap[2] / v_vphase

                ivgx[itheta, ip, :] = v_vgroup[:, 0]
                ivgy[itheta, ip, :] = v_vgroup[:, 1]
                ivgz[itheta, ip, :] = v_vgroup[:, 2]

        for i in range(3):
            if ax_ivp:
                ax_ivp.plot_surface( ivx[:, :, i], ivy[:, :, i], ivz[:, :, i], alpha=0.25)

            if ax_vg:
                ax_vg.plot_surface( ivgx[:, :, i], ivgy[:, :, i], ivgz[:, :, i], alpha=0.25)

        if ax_ivp:
            ax_ivp.set_xlabel(r"$1/v_x^{(p)}$ [s/km]", fontsize=8, labelpad=1)
            ax_ivp.set_ylabel(r"$1/v_y^{(p)}$ [s/km]", fontsize=8, labelpad=1)
            ax_ivp.set_zlabel(r"$1/v_z^{(p)}$ [s/km]", fontsize=8, labelpad=1)

        if ax_vg:
            ax_vg.set_xlabel(r"$v_x^{(g)}$ [km/s]", fontsize=8, labelpad=1)
            ax_vg.set_ylabel(r"$v_y^{(g)}$ [km/s]", fontsize=8, labelpad=1)
            ax_vg.set_zlabel(r"$v_z^{(g)}$ [km/s]", fontsize=8, labelpad=1)

        for ax in axs:
            for a in ("x", "y", "z"):
                ax.tick_params(axis=a, labelsize=8, pad=0)
            for t_ax in [ax.xaxis, ax.yaxis, ax.zaxis]:
                t_ax.line.set_linewidth(0.5)

            # ax.set_aspect('equal')

    def plot_bulk_dispersion_3D(self, pref):
        """
        Generate isocontour surfaces of the bulk dispersion in 3D k-space.
        """

        fig, axs = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
        ax_vp, ax_vg = axs

        self._add_3d_dispersion_curves_to_axes(ax_vp, ax_vg)

        fig.savefig(pref + "-bulkdisp3D.png")
        plt.close(fig)

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

        fig, axs = setup_bulk_dispersion_2D_plot()

        ax_sl, ax_vp, ax_vg, ax_ivp_3d = axs

        cm = "cool"  # Color map for polarisation coding
        self.add_bulk_slowness_curves_to_axes(pref, fig, ax_sl, ax_vp, ax_vg, cm, show_poln)

        #label = self.material_name
        if label is not None:
            ax_vp.text( -0.1, 1.1, label, fontsize=14, style="italic", transform=ax_sl.transAxes)

        self._add_3d_dispersion_curves_to_axes(ax_ivp_3d)

        fig.savefig(pref + "-bulkdisp.png")
        plt.close(fig)

    def add_bulk_slowness_curves_to_axes(self, pref, fig, ax_sl, ax_vp, ax_vg, cm, show_poln=True):
        npolpts = 28
        npolskip = 10  # make bigger
        npts = npolpts * npolskip  # about 1000
        v_kphi = np.linspace(0.0, np.pi * 2, npts)
        v_vel = np.zeros([npts, 3])
        v_velc = np.zeros([npts, 3])
        v_vgx = np.zeros([npts, 3])
        v_vgz = np.zeros([npts, 3])

        cmm = mpl.colormaps[cm]
        with open(pref + "-bulkdisp.dat", "w") as fout:
            fout.write(
                "#phi    kapx     kapz      vl          vs1         vs2         vlx      vly       vlz     vs1x     vs1y      vs1z       vs2x     vs2y     vs2z      k.v1    k.v2   k.v3\n"
            )

            kapcomp = np.zeros(3)
            ycomp = np.zeros(3)
            for ik, kphi in enumerate(v_kphi):
                # kapx = np.cos(kphi)
                # kapz = np.sin(kphi)
                # kapy = 0.0
                vkap = np.array([np.cos(kphi), 0.0, np.sin(kphi)])

                fout.write(f"{kphi:.4f}  {vkap[0]:+.4f}  {vkap[2]:+.4f}  ")

                # solve_christoffel returns:
                # eigvecs are sorted by phase velocity
                # v_vphase[m]:   |vphase| of modes m=1 to 3
                # vecs[:,m]:     evecs of modes m=1 to 3
                # v_vgroup[m,:]  vgroup of mode m, second index is x,y,z
                v_vphase, vecs, v_vgroup = solve_christoffel( vkap, self.stiffness_c_IJ, self.rho)

                v_vel[ik, :] = v_vphase  # phase velocity
                v_vgx[ik, :] = v_vgroup[:, 0]  # group velocity components
                v_vgz[ik, :] = v_vgroup[:, 2]

                ycomp = np.abs(vecs[1, :])  # $\unity \cdot u_i$
                kapcomp = np.abs( np.matmul(vkap, vecs))  # component of vkap along each evec
                v_velc[ik, :] = kapcomp  # phase velocity color by polarisation

                for iv in range(3):
                    fout.write(f"{v_vphase[iv]*1000:10.4f}  ")
                for iv in range(3):
                    fout.write( f"{vecs[0,iv]:7.4f}  {vecs[1,iv]:7.4f}   {vecs[2,iv]:7.4f}  ")
                fout.write(f"{kapcomp[0]:6.4f}  {kapcomp[1]:6.4f} {kapcomp[2]:6.4f}")

                fout.write("\n")

                if show_poln:
                    # Draw polarisation ball and stick notations
                    irad = 0.07 / v_vel[0, 0]  # length of polarisation sticks
                    rad = 0.07 * v_vel[0, 0]  # length of polarisation sticks
                    lwstick = 0.9
                    srad = 3  # diameter of polarisation dots
                    if ik % npolskip == 0:
                        for i in range(3):
                            radsl = 1 / v_vel[ik, i]
                            radvp = v_vel[ik, i]
                            polc = cmm(kapcomp[i])
                            polc = "gray"  # all black for now

                            ptm = radsl * np.array([np.cos(kphi), np.sin(kphi)])
                            pt0 = np.real(ptm - vecs[0:3:2, i] * irad)
                            pt1 = np.real(ptm + vecs[0:3:2, i] * irad)
                            ax_sl.plot( (pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                            ax_sl.plot( ptm[0], ptm[1], "o", c=polc, markersize=srad * ycomp[i])

                            ptm = radvp * np.array([np.cos(kphi), np.sin(kphi)])
                            pt0 = np.real(ptm - vecs[0:3:2, i] * rad)
                            pt1 = np.real(ptm + vecs[0:3:2, i] * rad)

                            ax_vp.plot( (pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                            ax_vp.plot( ptm[0], ptm[1], "o", c=polc, markersize=srad * ycomp[i])

        # the main curves for v_p, 1/v_p and v_g
        lw=.5
        for i in range(3):
            ax_vp.scatter( np.cos(v_kphi) * v_vel[:, i], np.sin(v_kphi) * v_vel[:, i],
                          c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm, lw=lw)

            ax_sl.scatter( np.cos(v_kphi) / v_vel[:, i], np.sin(v_kphi) / v_vel[:, i],
                          c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm, lw=lw)


            ax_vg.scatter( v_vgx[:, i], v_vgz[:, i], c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm, lw=lw)

        # Tick location seems to need help here
        for tax in [ax_vp.xaxis, ax_vp.yaxis, ax_vg.xaxis, ax_vg.yaxis]:
            tax.set_major_locator( ticker.MultipleLocator( 2.0 ))# , offset=0

        make_axes_square(np.abs(1 / v_vel).max(), ax_sl)
        make_axes_square(np.abs(v_vel).max(), ax_vp)
        make_axes_square(max(np.abs(v_vgx).max(), np.abs(v_vgz).max()), ax_vg)

        # Add radial speed grid
        v_theta = np.linspace(0, 2*np.pi, 300)
        for vr in range(1,8):
            ax_vp.plot( np.cos(v_theta) * vr , np.sin(v_theta) * vr, ':', c='gray', lw=.25)
            ax_vg.plot( np.cos(v_theta) * vr , np.sin(v_theta) * vr, ':', c='gray', lw=.25)
            ax_sl.plot( np.cos(v_theta) / vr , np.sin(v_theta) / vr, ':', c='gray', lw=.25)

        # fig.colorbar(mplcm.ScalarMappable(cmap=cm), ax=ax_vp, shrink=.5,
        #             pad=.025, location='top', label='$\hat{e} \cdot \hat{\kappa}$')

    def add_bulk_slowness_curves_to_axes_2x1(
        self, pref, fig, ax_sl, ax_vp, cm, mat1or2
    ):
        npolpts = 28
        npolskip = 10  # make bigger
        npts = npolpts * npolskip  # about 1000
        v_kphi = np.linspace(0.0, np.pi * 2, npts)
        v_vel = np.zeros([npts, 3])
        v_velc = np.zeros([npts, 3])
        # v_vgx = np.zeros([npts, 3])
        # v_vgz = np.zeros([npts, 3])

        cmm = mpl.colormaps[cm]

        kapcomp = np.zeros(3)
        ycomp = np.zeros(3)
        for ik, kphi in enumerate(v_kphi):
            vkap = np.array([np.cos(kphi), 0.0, np.sin(kphi)])

            # solve_christoffel returns:
            # eigvecs are sorted by phase velocity
            # v_vphase[m]:   |vphase| of modes m=1 to 3
            # vecs[:,m]:     evecs of modes m=1 to 3
            # v_vgroup[m,:]  vgroup of mode m, second index is x,y,z
            v_vphase, vecs, v_vgroup = solve_christoffel(vkap, self.stiffness_c_IJ, self.rho)

            v_vel[ik, :] = v_vphase  # phase velocity
            # v_vgx[ik, :] = v_vgroup[:,0]  # group velocity components
            # v_vgz[ik, :] = v_vgroup[:,2]

            ycomp = np.abs(vecs[1, :])  # $\unity \cdot u_i$
            kapcomp = np.abs(np.matmul(vkap, vecs))  # component of vkap along each evec
            # This causes both shear waves to have the same colour if there are two pure
            v_velc[ik, :] = kapcomp  # phase velocity color by polarisation.

            # Draw polarisation ball and stick notations
            irad = 0.07 / v_vel[0, 0]  # length of polarisation sticks
            rad = 0.07 * v_vel[0, 0]  # length of polarisation sticks
            lwstick = 0.9
            srad = 5  # diameter of polarisation dots
            if ik % npolskip == 0:
                for i in range(3):
                    radsl = 1 / v_vel[ik, i]
                    radvp = v_vel[ik, i]
                    polc = cmm(kapcomp[i])
                    polc = "k"  # all black for now

                    ptm = radsl * np.array([np.cos(kphi), np.sin(kphi)])
                    pt0 = np.real(ptm - vecs[0:3:2, i] * irad)
                    pt1 = np.real(ptm + vecs[0:3:2, i] * irad)
                    ax_sl.plot((pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                    ax_sl.plot(ptm[0], ptm[1], "o", c=polc, markersize=srad * ycomp[i])

                    ptm = radvp * np.array([np.cos(kphi), np.sin(kphi)])
                    pt0 = np.real(ptm - vecs[0:3:2, i] * rad)
                    pt1 = np.real(ptm + vecs[0:3:2, i] * rad)

                    # ax_vp.plot((pt0[0], pt1[0]), (pt0[1], pt1[1]), c=polc, lw=lwstick)
                    # ax_vp.plot(ptm[0], ptm[1], 'o', c=polc, markersize=srad*ycomp[i])

        # the main curves for 1/v_p and v_g
        for i in range(3):
            ax_sl.scatter(
                np.cos(v_kphi) / v_vel[:, i],
                np.sin(v_kphi) / v_vel[:, i],
                c=v_velc[:, i],
                vmin=0,
                vmax=1,
                s=0.5,
                cmap=cm,
            )

            # ax_vp.scatter(np.cos(v_kphi)*v_vel[:, i], np.sin(v_kphi) *
            #            v_vel[:, i], c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm)

            # ax_vg.scatter(v_vgx[:,i], v_vgz[:,i],  c=v_velc[:, i], vmin=0, vmax=1, s=0.5, cmap=cm)

        # Tick location seems to need help here
        # for tax in [ax_vp.xaxis, ax_vp.yaxis, ax_vg.xaxis, ax_vg.yaxis]:
        #   tax.set_major_locator(ticker.MultipleLocator(2.0, offset=0))

        make_axes_square(np.abs(1 / v_vel).max(), ax_sl)
        # make_axes_square(np.abs(v_vel).max(), ax_vp)
        # make_axes_square(max(np.abs(v_vgx).max(), np.abs(v_vgz).max()), ax_vg)

        cbar = fig.colorbar(
            mplcm.ScalarMappable(cmap=cm),
            ax=ax_sl,
            shrink=0.5,
            pad=0.025,
            location="right",
        )
        cbar.ax.tick_params(labelsize=6, width=0.25)

        cbar.outline.set_linewidth(1)
        cbar.set_label(
            label=f"Mat {mat1or2} " + r"$\hat{e} \cdot \hat{\kappa}$", fontsize=10
        )

    def make_crystal_axes_plot(self, pref):
        """Build crystal coordinates diagram using call to external asymptote application."""

        fn = tempfile.NamedTemporaryFile(suffix=".asy", mode="w+t", delete=False)

        asy_cmds = asy_draw_crystal_axes(self._crystal_axes)
        fn.write(asy_cmds)
        fn.close()

        # run .asy
        subprocess.run(["asy", fn.name, "-o", f"{pref}-crystal"], check=False)

    def plot_photoelastic_IJ(self, prefix, v_comps):
        """ Plot photoelastic tensor components as a function of rotation angle about y-axis.

        Args:
            prefix (str): Prefix for output file names.
            v_comps (list): List of strings of desired elements: "11", "12", "31" etc

        """

        npts = 200

        fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(4,4))
        # fig.subplots_adjust(hspace=.35, wspace=0)

        d_p_vecs = {}  # map of "IJ" strings -> (el_IJ, np.zeros(npts))

        for s_elt in v_comps:
            if len(s_elt) !=2:
                reporting.report_and_exit('Bad photoelastic tensor index: {s_elt}.')

            el_IJ=(int(s_elt[0]), int(s_elt[1]))

            if set(el_IJ) <= set(range(1,7)): # all elements are in range 1-6
            #if el_IJ[0] not in range(1,7) or el_IJ[1] not in range(1,7): #TODO express as a set operation
                reporting.report_and_exit('Bad photoelastic tensor index: {s_elt}.')

            d_p_vecs[s_elt] = (el_IJ, np.zeros(npts))

        mat0 = copy.deepcopy(self)
        haty=np.array([0,1,0])

        v_phi = np.linspace(0.0, np.pi * 2, npts)
        for iphi, phi in enumerate(v_phi): # for all angles around the circle
            t_mat = copy.deepcopy(mat0)
            t_mat.rotate(haty, phi)

            for (k,v) in d_p_vecs.items(): # extract the desired p_IJ components
                (I,J) = v[0]
                v[1][iphi] = t_mat.photoel_p_IJ[I,J]

        for (k,v) in d_p_vecs.items():
            (I,J) = v[0]
            v_pIJ = v[1]

            lab = '$p_{' + f'{I},{J}' + '}$'
            plt.polar(v_phi, v_pIJ,  label=lab, lw=.5)

        ax.set_rmax(0.4)
        ax.set_rmin(-0.2)
        ax.set_rticks(np.arange(-0.2,0.4+.0001,.1))

        ax.set_rlabel_position(90)
        ax.set_thetagrids(np.arange(0,360-1e-5,30))
        ax.grid(True)

        #ax.set_ylim(-.2,.4)
        ax.legend()


def setup_bulk_dispersion_2D_plot():
    """Plots both slowness and ray normal contours."""

    fig, axs = plt.subplots(2, 2, figsize=(7, 6), dpi=300)
    fig.subplots_adjust(hspace=0.35, wspace=0)

    ax_vp, ax_sl, ax_vg = axs[0, 0], axs[0, 1], axs[1, 0]

    axs[1, 1].set_axis_off()  # Hide axis 2,2

    axs[1, 1].remove()
    ax_ivp3d = fig.add_subplot(2, 2, 4, projection="3d")

    ax_vp.set_xlabel(r"$v^{(p)}_{x}$ [km/s]")
    ax_vp.set_ylabel(r"$v^{(p)}_{z}$ [km/s]")

    ax_sl.set_xlabel(r"$1/v^{(p)}_{x}$ [s/km]")
    ax_sl.set_ylabel(r"$1/v^{(p)}_{z}$ [s/km]")

    ax_vg.set_xlabel(r"$v^{(g)}_{x}$ [km/s]")
    ax_vg.set_ylabel(r"$v^{(g)}_{z}$ [km/s]")

    for ax in axs.flat[:3]:  # Don't write to axis 2,2
        ax.axhline(0, c="gray", lw=0.5)
        ax.axvline(0, c="gray", lw=0.5)
        ax.tick_params(width=0.5)
        for item in ( [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(10)
        for t_ax in ["top", "bottom", "left", "right"]:
            ax.spines[t_ax].set_linewidth(0.5)
    axs = ax_sl, ax_vp, ax_vg, ax_ivp3d
    return fig, axs


def setup_bulk_dispersion_2D_plot_2x1():
    """Plots both slowness and ray normal contours."""

    fig, axs = plt.subplots(1, 1, figsize=(6, 4))
    # fig.subplots_adjust(hspace=.35, wspace=0)
    axs = (axs,)
    # ax_sl, ax_vg = axs
    ax_sl = axs[0]

    # ax_sl, ax_vp, ax_vg = axs[0,0], axs[0,1], axs[1,0]

    # axs[1,1].set_axis_off()  # Hide axis 2,2

    # axs[1,1].remove()
    # ax_ivp3d = fig.add_subplot(2,2,4, projection='3d')
    ax_sl.set_xlabel(r"$1/v^{(p)}_{x}$ [s/km]")
    ax_sl.set_ylabel(r"$1/v^{(p)}_{z}$ [s/km]")
    # ax_vp.set_xlabel(r'$v^{(p)}_{x}$ [s/km]')
    # ax_vp.set_ylabel(r'$v^{(p)}_{z}$ [s/km]')
    # ax_vg.set_xlabel(r'$v^{(g)}_{x}$ [km/s]')
    # ax_vg.set_ylabel(r'$v^{(g)}_{z}$ [km/s]')

    for ax in axs:  # Don't write to axis 2,2
        ax.axhline(0, c="gray", lw=0.5)
        ax.axvline(0, c="gray", lw=0.5)
        ax.tick_params(width=0.5)
        for item in (
            [ax.title, ax.xaxis.label, ax.yaxis.label]
            + ax.get_xticklabels()
            + ax.get_yticklabels()
        ):
            item.set_fontsize(12)
        for t_ax in ["top", "bottom", "left", "right"]:
            ax.spines[t_ax].set_linewidth(0.5)
    # axs = ax_sl, ax_vp, ax_vg, ax_ivp3d
    return fig, axs


def compare_bulk_dispersion(mat1, mat2, pref):
    fig, axs = setup_bulk_dispersion_2D_plot_2x1()

    # ax_sl, ax_vg = axs
    ax_sl = axs[0]
    ax_vg = None

    cm1 = "cool"  # Color map for polarisation coding
    cm2 = "autumn"  # Color map for polarisation coding

    mat1.add_bulk_slowness_curves_to_axes_2x1(pref + "_mat1", fig, ax_sl, ax_vg, cm1, 1)
    mat2.add_bulk_slowness_curves_to_axes_2x1(pref + "_mat2", fig, ax_sl, ax_vg, cm2, 2)

    ax_sl.text(
        0.05,
        1.15,
        f"Mat 1: {mat1.material_name}",
        fontsize=14,
        style="italic",
        transform=ax_sl.transAxes,
    )
    ax_sl.text(
        0.05,
        1.05,
        f"Mat 2: {mat2.material_name}",
        fontsize=14,
        style="italic",
        transform=ax_sl.transAxes,
    )

    plt.savefig(pref + "-compare-bulkdisp.png")



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


def asy_draw_crystal_axes(crystal_axes):
    (va, vb, vc) = crystal_axes
    s_avec = "(" + ",".join(map(str, va)) + ")"
    s_bvec = "(" + ",".join(map(str, vb)) + ")"
    s_cvec = "(" + ",".join(map(str, vc)) + ")"

    s1 = r"""
settings.outformat='png';
settings.render=8;
import three;
import graph3;

size(2cm,0);
defaultpen(fontsize(7pt));
defaultpen(.2);

real axlen=1.25;
int arrsize=3;
real blen=.5;

//currentprojection=orthographic(1,1,1);
currentprojection=oblique;


draw(O--2X, black, Arrow3(arrsize), L=Label("$\hat{x}$", position=EndPoint));
draw(O--2Y, black, Arrow3(arrsize), L=Label("$\hat{y}$", position=EndPoint));
draw(O--3Z, black, Arrow3(arrsize), L=Label("$\hat{z}$", position=EndPoint));

draw(O-- -2X, gray);
draw(O-- -2Y, gray);
draw(O-- -2Z, gray);


//label("$\hat{x}$", 3X*1.1);
//label("$\hat{y}$", 3Y*1.1);
//label("$\hat{z}$", 3Z*1.1);

draw(box((-1,-.5,-2)*blen,(1,.5,2)*blen),blue);
"""

    s2 = f"""triple avec={s_avec};
triple bvec={s_bvec};
triple cvec={s_cvec};
"""

    s3 = """triple corig=(0,.5,2)*blen;
draw(corig--avec+corig, red, Arrow3(arrsize), L=Label("$c_x$", position=EndPoint));
draw(corig--bvec+corig, red, Arrow3(arrsize), L=Label("$c_y$", position=EndPoint));
draw(corig--cvec+corig, red, Arrow3(arrsize), L=Label("$c_z$", position=EndPoint));

triple k0=(1,-1,-1);
triple k1=k0+(0,0,2);

draw(k0--k1,green, Arrow3(arrsize), L=Label("$k$"));
"""

    return s1 + s2 + s3


def make_axes_square(ext0, ax):
    ext = 1.1 * ext0
    ax.set_xlim(-ext, ext)
    ax.set_ylim(-ext, ext)
    ax.set_aspect("equal")
