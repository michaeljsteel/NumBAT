
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
import math
import subprocess
import numbers
import collections.abc as abc

import reporting

def is_real_number(x):
    return isinstance(x, numbers.Real)  # need numpy.int32, int36, float64, etc

def is_float_pair(x):
    #print('float pair test', x,
    #isinstance(x, abc.Sequence) ,len(x)==2, is_real_number(x[0]), is_real_number(x[1]))

    is_arr = isinstance(x, abc.Sequence) or isinstance(x, np.ndarray)
    return is_arr and len(x)==2 and is_real_number(x[0]) and is_real_number(x[1])

def almost_zero(x, tol=1e-10):
    return math.isclose(x, 0, abs_tol=tol)

def almost_unity(x, tol=1e-10):
    return math.isclose(x, 1, abs_tol=tol)

def np_min_max(v):
    return np.min(v), np.max(v)


def int2D(mat, dx=1.0, dy=1.0):
    return np.sum(np.sum(mat)) * dx * dy

def int2D_trapz(mat, dx=1.0, dy=1.0):
    return np.trapz(np.trapz(mat)) * dx * dy

# This does not work well on non-rectangular domains where there are sudden
# jumps at outer boundaries
def int2D_simp(mat, dx=1.0, dy=1.0):
    return sciint.simpson(sciint.simpson(mat)) * dx * dy

def process_fortran_return(resm, msg):
    """Check return values of any Fortran function with errco, emsg style return codes"""

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
    s_ind = indent * ' '
    s_out= s_ind + s_in
    s_out=s_out.replace('\n', '\n'+s_ind)
    return s_out

def run_subprocess(cmd, proc_name, cwd='', exit_on_fail=True):
    try:
        comp_stat = subprocess.run(cmd, cwd=cwd)
        if comp_stat.returncode and exit_on_fail:
            tcmd = ' '.join(cmd)
            reporting.report_and_exit(f'{proc_name} call failed executing: "{tcmd}".')
    except subprocess.CalledProcessError as e:
        reporting.report_and_exit(f"Error occurred while running {proc_name}: {e}")

    return comp_stat.returncode


def f2f_with_subs(fn_in, fn_out, d_subs):

    with open(fn_in, 'r') as fin:
        conv_tmp = fin.read()

    for k, v in d_subs.items():
        conv = conv_tmp.replace(k, v)

    with open(fn_out, 'w') as fout:
        fout.write(conv)
