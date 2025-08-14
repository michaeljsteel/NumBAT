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


from nbtypes import CrystalGroup, unit_x, unit_y, unit_z

def construct_crystal_for_symmetry_class(crystal_class, symmetry_group=None):
    match crystal_class:

        case CrystalGroup.Trigonal:
            return _construct_crystal_trigonal()

        case CrystalGroup.Hexagonal:
            return _construct_crystal_hexagonal()

        case CrystalGroup.Cubic:
            return _construct_crystal_cubic()




def _construct_crystal_cubic():
    """
    Return lists of independent and dependent elements for cubic crystals.
    """

    # plain cartesian axes
    C_indep_elts = ["11", "12", "44"]
    C_dep_elts = {
        "13": [("12", 1)],
        "21": [("12", 1)],
        "22": [("11", 1)],
        "23": [("12", 1)],
        "31": [("12", 1)],
        "32": [("12", 1)],
        "33": [("11", 1)],
        "55": [("44", 1)],
        "66": [("44", 1)],
    }

    pIJ_indep_elts = ["11", "12", "44"]
    pIJ_dep_elts = {
        "13": [("12", 1)],
        "21": [("12", 1)],
        "22": [("11", 1)],
        "23": [("12", 1)],
        "31": [("12", 1)],
        "32": [("12", 1)],
        "33": [("11", 1)],
        "55": [
            ("44", 1)
        ],  # According to Powell, for Oh group, these are distinct elements, but no one seems to quote them
        "66": [("44", 1)],
    }

    return (unit_x, unit_y, unit_z), C_indep_elts, C_dep_elts, pIJ_indep_elts, pIJ_dep_elts


def _construct_crystal_trigonal():
    """
    Return lists of independent and dependent elements for trigonal crystals.
    """

    # Good source for these rules is the supp info of doi:10.1364/JOSAB.482656 (Gustavo surface paper)

    C_indep_elts = ["11", "12", "13", "14", "33", "44"]
    C_dep_elts = {
        "21": [("12", 1)],
        "22": [("11", 1)],
        "23": [("13", 1)],
        "24": [("14", -1)],
        "31": [("13", 1)],
        "32": [("13", 1)],
        "41": [("14", 1)],
        "42": [("14", -1)],
        "55": [("44", 1)],
        "56": [("14", 1)],
        "65": [("14", 1)],
        "66": [("11", 0.5), ("12", -0.5)],
    }

    pIJ_indep_elts = ["11", "12", "13", "14", "31", "33", "41", "44"]
    pIJ_dep_elts = {
        "21": [("12", 1)],
        "22": [("11", 1)],
        "23": [("13", 1)],
        "24": [("14", -1)],
        "32": [("31", 1)],
        "42": [("41", -1)],
        "55": [("44", 1)],
        "56": [("41", 1)],
        "65": [("14", 1)],
        "66": [("11", 0.5), ("12", -0.5)],
    }

    return (unit_x, unit_y, unit_z), C_indep_elts, C_dep_elts, pIJ_indep_elts, pIJ_dep_elts



def _construct_crystal_hexagonal():
    """
    Return lists of independent and dependent elements for hexagonal crystals.
    """

    C_indep_elts = ["11", "12", "13", "33", "44"]
    C_dep_elts = {
        "21": [("12", 1)],
        "22": [("11", 1)],
        "23": [("13", 1)],
        "31": [("13", 1)],
        "32": [("13", 1)],
        "55": [("44", 1)],
        "66": [("11", 0.5), ("12", -0.5)],
    }

    pIJ_indep_elts = ["11", "12", "13", "14", "31", "33", "41", "44"]
    pIJ_dep_elts = {
        "21": [("12", 1)],
        "22": [("11", 1)],
        "23": [("13", 1)],
        "24": [("14", -1)],
        "32": [("31", 1)],
        "42": [("41", -1)],
        "55": [("44", 1)],
        "56": [("41", 1)],
        "65": [("14", 1)],
        "66": [("11", 0.5), ("12", -0.5)],
    }
    return (unit_x, unit_y, unit_z), C_indep_elts, C_dep_elts, pIJ_indep_elts, pIJ_dep_elts


def construct_piezo_elt_mappings(sym_group):

    # See Auld vol 1, appendix A for the rules for piezoelectric elements.

    # first index of elts_dep is d_ijk multipliers, second is e_ijk multipliers
    elts_dep = {}

    match sym_group:


        # Triclinic systems
        case "1":
            elts_indep = ("x1", "x2", "x3", "x4", "x5", "x6",
                          "y1", "y2", "y3", "y4", "y5", "y6",
                          "z1", "z2", "z3", "z4", "z5", "z6"
                          )


        # Monoclinic systems
        case "2":
            elts_indep = ("x4", "x6", "y1", "y2", "y3", "y5", "z4", "z6")

        case "m":
            elts_indep = ("x1", "x2", "x3", "x5", "y4", "y6", "z1", "z2", "z3", "z5")


        # Orthorhombic systems
        case "222":
            elts_indep = ("x4", "y5", "z6")

        case "2mm":
            elts_indep = ("x5", "y4", "z1", "z2", "z3")


        # Tetragonal systems
        case "4'":
            elts_indep = ("x4", "x5", "z1", "z6")
            elts_dep = {
                "y4": ("x5", -1.0, -1.0),
                "y5": ("x4", 1.0, 1.0),
                "z2": ("z1", -1.0, -1.0),
             }

        case "4":
            elts_indep = ("x4", "x5", "z1", "z3")
            elts_dep = {
                "y4": ("x5", 1.0, 1.0),
                "y5": ("x4", -1.0, -1.0),
                "z2": ("z1", 1.0, 1.0),
             }

        case "4'2m":
            elts_indep = ("x4","z6")
            elts_dep = {
                "y4": ("x4", 1.0, 1.0),
             }

        case "422":
            elts_indep = ("x4",)
            elts_dep = {
                "y5": ("x4", -1.0, -1.0),
             }

        case "4mm":
            elts_indep = ("x5", "z1", "z3")
            elts_dep = {
                "y4": ("x5", 1.0, 1.0),
                "y4": ("x5", 1.0, 1.0),
                "z2": ("z1", 1.0, 1.0),
             }


        # Trigonal systems
        case "3":
            elts_indep = ("x1", "x4", "x5", "y2", "z1", "z3")
            elts_dep = {
                "x2": ("x1", -1.0, -1.0),
                "x6": ("y2", -2.0, -1.0),
                "y1": ("y2", -1.0, -1.0),
                "y4": ("x5", 1.0, 1.0),
                "y5": ("x4", -1.0, -1.0),
                "y6": ("x1", -2.0, -1.0),
                "z2": ("z1", 1.0, 1.0),
        }

        case "32":
            elts_indep = ("x1", "x4")
            elts_dep = {
                "x2": ("x1", -1.0, -1.0),
                "y5": ("x4", -1.0, -1.0),
                "y6": ("x1", -2.0, -1.0),
            }

        case "3m":
            elts_indep = ("x5", "y2", "z1", "z3")
            elts_dep = {
                "x6": ("y2", -2.0, -1.0),
                "y1": ("y2", -1.0, -1.0),
                "y4": ("x5", 1.0, 1.0),
                "z2": ("z1", 1.0, 1.0),
        }


        # Hexagonal systems
        case "6":
            elts_indep = ("x4", "x5", "z1", "z3")
            elts_dep = {"y4": ("x5", 1.0, 1.0),
                        "y5": ("x4", -1.0, -1.0),
                        "z2": ("z1", 1.0, 1.0)}

        case "622":
            elts_indep = ("x4",)
            elts_dep = {"y5": ("x4", -1.0, -1.0)}


        case "6mm":
            elts_indep = ("x5", "z1", "z3")
            elts_dep = {"y4": ("x5", 1.0, 1.0), "z2": ("z1", 1.0, 1.0)}

        case "6'":
            elts_indep = ("x1", "y2")
            elts_dep = {"x2": ("x1", -1.0, -1.0),
                        "x6": ("y2", -2.0, -1.0),
                        "y1": ("y2", -1.0, -1.0),
                        "y6": ("x1", -2.0, -1.0)
                        }

        case "6'm2":
            elts_indep = ("x1",)
            elts_dep = {"x2": ("x1", -1.0, -1.0),
                        "y6": ("x1", -2.0, -1.0)
                       }

        # Cubic systems
        case "23" | "4'3m":
            elts_indep = ("x4",)
            elts_dep = {"y5": ("x4", 1.0, 1.0), "z6": ("x4", 1.0, 1.0)}

        case _:
            raise ValueError(f"Unknown symmetry group in piezo matrix construction: {sym_group}")


    return elts_indep, elts_dep
