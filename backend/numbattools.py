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
from scipy.integrate import simpson
import math
import subprocess
import numbers
import collections.abc as abc

import reporting

def is_real_number(x):
    """Return True if x is a real number (including numpy types)."""
    return isinstance(x, numbers.Real)

def is_float_pair(x):
    """Return True if x is a sequence or array of two real numbers."""
    is_arr = isinstance(x, abc.Sequence) or isinstance(x, np.ndarray)
    return is_arr and len(x)==2 and is_real_number(x[0]) and is_real_number(x[1])

def almost_zero(x, tol=1e-10):
    """Return True if x is close to zero within absolute tolerance tol."""
    return math.isclose(x, 0, abs_tol=tol)

def almost_unity(x, tol=1e-10):
    """Return True if x is close to one within absolute tolerance tol."""
    return math.isclose(x, 1, abs_tol=tol)

def np_min_max(v):
    """Return the minimum and maximum of a numpy array v as a tuple."""
    return np.min(v), np.max(v)

def int2D(mat, dx=1.0, dy=1.0):
    """Return the 2D integral of mat using the rectangle rule with spacings dx, dy."""
    return np.sum(mat) * dx * dy

def int2D_trapz(mat, dx=1.0, dy=1.0):
    """Return the 2D integral of mat using the trapezoidal rule with spacings dx, dy."""
    return np.trapz(np.trapz(mat)) * dx * dy

# This does not work well on non-rectangular domains where there are sudden
# jumps at outer boundaries
def int2D_simp(mat, dx=1.0, dy=1.0):
    """Return the 2D integral of mat using Simpson's rule with spacings dx, dy."""
    return simpson(simpson(mat)) * dx * dy

def signed_log10(mat, logmax=5, tol=1e-14, power=1.0):
    """Return a signed log10 of mat, shifted by logmax, with small values clipped to zero."""
    lmat = np.log10(np.abs(mat)+tol) + logmax
    lmat = np.where(lmat>0, lmat, 0)**power
    lmat *= np.sign(mat)
    return lmat


def chopmat(mat, chop=True, rtol=1e-10, atol=1e-12):
    """Set all values below an absolute or relative tolerance to zero."""

    chopmat = mat.copy()

    if chop:
        maxval = np.max(np.abs(chopmat))
        chopmat[np.abs(chopmat) < atol] = 0
        chopmat[np.abs(chopmat) < rtol*maxval] = 0

    return chopmat


def polish_float(x, atol=1e-10):
    for v in range(-10,11,1):
        if np.isclose(x, v, atol):
            return v

    for v in (.1,.2,.3,.4,.5,.6,.7,.9):
        if np.isclose(x, v, atol):
            return v
        if np.isclose(x, -v, atol):
            return -v
    return x

def vchopmattoint(vec, atol=1e-10):
    """Set all values with an absolute tolerance of one to one."""
    chopvec = vec.copy()
    for r in range(len(vec)):
        chopvec[r] = polish_float(chopvec[r])
    return chopvec

def chopmattoint(mat, atol=1e-10):
    """Set all values with an absolute tolerance of one to one."""

    chopmat = mat.copy()
    rs,cs = mat.shape

    for r in range(rs):
        for c in range(cs):
            for val in range(-10,11,1):
                if np.isclose(chopmat[r,c], val, atol):
                    chopmat[r,c]=val

            for val in (.1,.2,.3,.4,.5,.6,.7,.9):
                if np.isclose(chopmat[r,c], val, atol):
                    chopmat[r,c]=val
                if np.isclose(chopmat[r,c], -val, atol):
                    chopmat[r,c]=-val


            #if np.isclose(chopmat[r,c], -1.0, atol):
            #    chopmat[r,c]=-1.0

    return chopmat



def process_fortran_return(resm, msg):
    """
    Check return values of any Fortran function with errco, emsg style return codes.
    If an error is detected, report and exit. Otherwise, return the result tuple without error info.
    """
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

def indent_string(s_in, indent=2):
    """Indent every line in s_in by the specified number of spaces."""
    s_ind = indent * ' '
    s_out= s_ind + s_in
    s_out=s_out.replace('\n', '\n'+s_ind)
    return s_out

def run_subprocess(cmd, proc_name, cwd='', exit_on_fail=True):
    """
    Run a subprocess with the given command and working directory.
    If exit_on_fail is True, report and exit on failure.
    Returns the process return code.
    """
    try:
        comp_stat = subprocess.run(cmd, cwd=cwd)
        if comp_stat.returncode and exit_on_fail:
            tcmd = ' '.join(cmd)
            reporting.report_and_exit(f'{proc_name} call failed executing: "{tcmd}".')
    except subprocess.CalledProcessError as e:
        reporting.report_and_exit(f"Error occurred while running {proc_name}: {e}")

    return comp_stat.returncode

def f2f_with_subs(fn_in, fn_out, d_subs):
    """
    Read a file, apply string substitutions from d_subs, and write to a new file.
    Args:
        fn_in (str): Input filename.
        fn_out (str): Output filename.
        d_subs (dict): Dictionary of string substitutions {old: new}.
    """
    with open(fn_in, 'r') as fin:
        conv = fin.read()

    for k, v in d_subs.items():
        conv = conv.replace(k, v)

    with open(fn_out, 'w') as fout:
        fout.write(conv)
