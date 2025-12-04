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
from typing import Any, Sequence, Mapping, Optional, Tuple, Union
from numpy.typing import NDArray

import reporting

def is_real_number(x: Any) -> bool:
    """Return True if ``x`` is a real scalar (Python or numpy real)."""
    return isinstance(x, numbers.Real)

def is_float_pair(x: Any) -> bool:
    """Return True if ``x`` is a length-2 sequence/array of real numbers."""
    is_arr = isinstance(x, abc.Sequence) or isinstance(x, np.ndarray)
    return bool(is_arr and len(x) == 2 and is_real_number(x[0]) and is_real_number(x[1]))

def almost_zero(x: float, tol: float = 1e-10) -> bool:
    """Return True if ``x`` is close to zero within absolute tolerance ``tol``."""
    return math.isclose(x, 0.0, abs_tol=tol)

def almost_unity(x: float, tol: float = 1e-10) -> bool:
    """Return True if ``x`` is close to one within absolute tolerance ``tol``."""
    return math.isclose(x, 1.0, abs_tol=tol)

def np_min_max(v: NDArray[Any]) -> Tuple[Any, Any]:
    """Return ``(min(v), max(v))`` for a numpy array ``v``."""
    return np.min(v), np.max(v)

def int2D(mat: NDArray[np.number], dx: float = 1.0, dy: float = 1.0) -> np.number:
    """Return 2D integral of ``mat`` via rectangle rule with spacings ``dx`` and ``dy``."""
    return np.sum(mat) * dx * dy

def int2D_trapz(mat: NDArray[np.number], dx: float = 1.0, dy: float = 1.0) -> np.number:
    """Return 2D integral of ``mat`` using nested trapezoidal rule (spacing ``dx``, ``dy``)."""
    # NOTE: numpy.trapezoid deprecated; code retained per original intent.
    return np.trapezoid(np.trapezoid(mat)) * dx * dy  # type: ignore[attr-defined]

def check_magnitude(x: float, name: str, minval: float, maxval: float) -> None:
    """Validate that the magnitude of a value lies within bounds.

    Parameters
    ----------
    x
        Input value to be checked.
    name
        Name of the quantity (used in the error message).
    minval, maxval
        Allowed bounds on |x| (inclusive). Units must match the caller's context.

    Raises
    ------
    ValueError
        If |x| is outside [minval, maxval].
    """
    if not (minval <= abs(x) <= maxval):
        raise ValueError(f"{name} value {x} out of expected range [{minval}, {maxval}]")


# This does not work well on non-rectangular domains where there are sudden
# jumps at outer boundaries
def int2D_simp(mat: NDArray[np.number], dx: float = 1.0, dy: float = 1.0) -> float:
    """Return 2D integral of ``mat`` using Simpson's rule (spacing ``dx``, ``dy``)."""
    return float(simpson(simpson(mat)) * dx * dy)

def signed_log10(mat: Union[NDArray[np.number], float], logmax: float = 5, tol: float = 1e-14, power: float = 1.0) -> Union[NDArray[np.number], float]:
    """Return signed log10 transform of ``mat`` with shift and clipping.

    ``mat`` may be a scalar or array; result mirrors the input shape.
    """
    lmat = np.log10(np.abs(mat) + tol) + logmax
    lmat = np.where(lmat > 0, lmat, 0) ** power
    lmat = lmat * np.sign(mat)
    return lmat


def chopmat(mat: NDArray[np.number], chop: bool = True, rtol: float = 1e-10, atol: float = 1e-12) -> NDArray[np.number]:
    """Return a copy of ``mat`` with small values zeroed (absolute or relative tolerance)."""
    out = mat.copy()
    if chop:
        maxval = np.max(np.abs(out))
        out[np.abs(out) < atol] = 0
        out[np.abs(out) < rtol * maxval] = 0
    return out


def polish_float(x: float, atol: float = 1e-10) -> float:
    for v in range(-10,11,1):
        if np.isclose(x, v, atol):
            return v

    for v in (.1,.2,.3,.4,.5,.6,.7,.9):
        if np.isclose(x, v, atol):
            return v
        if np.isclose(x, -v, atol):
            return -v
    return x

def vchopmattoint(vec: NDArray[np.number], atol: float = 1e-10) -> NDArray[np.number]:
    """Return a copy with elements rounded to nearby canonical integers/fractions."""
    out = vec.copy()
    for r in range(len(vec)):
        out[r] = polish_float(float(out[r]), atol)
    return out

def chopmattoint(mat: NDArray[np.number], atol: float = 1e-10) -> NDArray[np.number]:
    """Return a copy with elements rounded to nearby canonical integers/fractions."""
    out = mat.copy()
    rs, cs = out.shape
    for r in range(rs):
        for c in range(cs):
            for val in range(-10, 11, 1):
                if np.isclose(out[r, c], val, atol):
                    out[r, c] = val
            for val in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9):
                if np.isclose(out[r, c], val, atol):
                    out[r, c] = val
                if np.isclose(out[r, c], -val, atol):
                    out[r, c] = -val
    return out



def process_fortran_return(resm: Sequence[Any], msg: str) -> Optional[Sequence[Any]]:
    """Validate Fortran-style return tuple with (.., err_code, err_msg).

    Returns the data portion (all but last two entries) or None if error and handled.
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

def indent_string(s_in: str, indent: int = 2) -> str:
    """Return ``s_in`` with each line prefixed by ``indent`` spaces."""
    s_ind = indent * ' '
    s_out = s_ind + s_in
    s_out = s_out.replace('\n', '\n' + s_ind)
    return s_out

def run_subprocess(cmd: Sequence[str], proc_name: str, cwd: str = '', exit_on_fail: bool = True) -> int:
    """Run a subprocess command list and return its return code.

    If ``exit_on_fail`` is True, reports and exits on non-zero status.
    """
    try:
        comp = subprocess.run(cmd, cwd=cwd)
    except subprocess.CalledProcessError as e:
        reporting.report_and_exit(f"Error occurred while running {proc_name}: {e}")
        return -1  # unreachable if report_and_exit terminates, but satisfies type checker
    if comp.returncode and exit_on_fail:
        tcmd = ' '.join(cmd)
        reporting.report_and_exit(f'{proc_name} call failed executing: "{tcmd}".')
        return comp.returncode  # same note as above
    return comp.returncode

def f2f_with_subs(fn_in: str, fn_out: str, d_subs: Mapping[str, str]) -> None:
    """Read file ``fn_in``, apply substitutions, write result to ``fn_out``."""
    with open(fn_in, 'r') as fin:
        conv = fin.read()

    for k, v in d_subs.items():
        conv = conv.replace(k, v)

    with open(fn_out, 'w') as fout:
        fout.write(conv)
