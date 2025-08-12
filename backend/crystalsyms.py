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

def construct_crystal_for_symmetry_class(crystal_class):
    match crystal_class:
        case CrystalGroup.Trigonal:
            return construct_crystal_trigonal()
        case CrystalGroup.Hexagonal:
            return construct_crystal_hexagonal()
        case CrystalGroup.Cubic:
            return construct_crystal_cubic()




def construct_crystal_cubic():
    """
    Return lists of independent and dependent elements for cubic crystals.
    """

    # plain cartesian axes
    #self.set_crystal_axes(unit_x, unit_y, unit_z)
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


def construct_crystal_trigonal():
    """
    Return lists of independent and dependent elements for trigonal crystals.
    """

    # Good source for these rules is the supp info of doi:10.1364/JOSAB.482656 (Gustavo surface paper)

    #self.set_crystal_axes(unit_x, unit_y, unit_z)

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



def construct_crystal_hexagonal():
    """
    Return lists of independent and dependent elements for hexagonal crystals.
    """

    #self.set_crystal_axes(unit_x, unit_y, unit_z)

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
