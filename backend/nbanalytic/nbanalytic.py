# Copyright (C) 2017-2025  Michael Steel

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


import math
import sys

import matplotlib.pyplot as plt
import matplotlib.colors as mpcol
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

import numpy as np
#import numpy.linalg as npla

import scipy.optimize as sciopt
import scipy.signal
import scipy.special as sp

from plotting.plottools import get_rgb_for_poln
import reporting
import numtools.rootbracket as ntrb
from nbtypes import SI_um

twopi = 2 * math.pi

import IPython.display as ipydisp




#######################################################################################


# class ElasticSlab:
#     """Elastic slab waveguide solver for isotropic materials.
#     Finds the dispersion of the Lamb modes of isolated flat plate
#     All units are in SI .
#     """

#     def __init__(self, mat_s, mat_f, mat_c, wid):
#         self._mats = mat_s  # substrate material
#         self._matf = mat_f  # film material
#         self._matc = mat_c  # cover material
#         self.width = wid  # width


