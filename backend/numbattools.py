
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


import numpy as np
import scipy.integrate as sciint


import reporting


def almost_zero(x, tol=1e-10):
    return abs(x) < tol


def almost_unity(x, tol=1e-10):
    return abs(1 - x) < tol


def np_min_max(v):
    return np.min(v), np.max(v)


def int2d(mat, dx=1.0, dy=1.0):
    return np.sum(np.sum(mat)) * dx * dy

def int2d_trapz(mat, dx=1.0, dy=1.0):
    return np.trapz(np.trapz(mat)) * dx * dy

# This does not work well on non-rectangular domains where there are sudden
# jumps at outer boundaries
def int2d_simp(mat, dx=1.0, dy=1.0):
    return sciint.simpson(sciint.simpson(mat)) * dx * dy

def process_fortran_return(resm, msg):

    fort_err, fort_mesg = resm[-2:]
    if fort_err:
        fort_mesg = str(fort_mesg, "utf-8")  # fort_mesg comes back as a byte string.
        reporting.report_and_exit(
            f"Fortran error in {msg}: \n"
            f" NumBAT Fortran error code = {fort_err}. \n Message: \n {fort_mesg}"
        )
        return None  # dummy
    else:  # everything is fine
        return resm[:-2]


