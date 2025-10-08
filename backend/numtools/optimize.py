import scipy.optimize as sciopt

class BracketError(Exception):
    pass

def bracket_minima(v_x, v_fx):
    """
    Find all local minima in a discrete function by identifying bracketed intervals.

    Searches through the discrete function values to find points where a function value
    is smaller than both its left and right neighbors, creating brackets around each
    local minimum for subsequent optimization.

    Parameters
    ----------
    v_x : array-like
        Array of x-coordinates (independent variable values).
    v_fx : array-like
        Array of function values f(x) corresponding to each x-coordinate.
        Must have the same length as v_x.

    Returns
    -------
    list of tuples
        List of bracket tuples, where each tuple contains three (x, f(x)) pairs:
        ((x_left, f_left), (x_center, f_center), (x_right, f_right))
        representing a bracketed minimum with f_center < f_left and f_center < f_right.

    Notes
    -----
    This function only identifies brackets around local minima; it does not perform
    the actual minimization. Use find_bracketed_minima() to refine the minima locations.

    Examples
    --------
    >>> x = np.array([0, 1, 2, 3, 4, 5])
    >>> f = np.array([5, 2, 1, 3, 0, 2])
    >>> brackets = bracket_minima(x, f)
    >>> len(brackets)  # Should find minima around x=2 and x=4
    2
    """
    l_braks = []
    for iv in range(1, len(v_x) - 1):
        fm1 = v_fx[iv - 1]
        f0 = v_fx[iv]
        fp1 = v_fx[iv + 1]

        if f0 < fm1 and f0 < fp1:
            l_braks.append(
                ((v_x[iv - 1], fm1), (v_x[iv], f0), (v_x[iv + 1], fp1))
            )
    return l_braks


def find_bracketed_minima(minfunc, l_braks):
    """
    Refine bracketed minima using scalar optimization.

    Takes brackets identified by bracket_minima() and uses scipy's minimize_scalar
    to find the precise location and value of each minimum within its bracket.

    Parameters
    ----------
    minfunc : callable
        The function to minimize. Should accept a scalar argument and return a scalar.
    l_braks : list of tuples
        List of bracket tuples from bracket_minima(), where each tuple contains
        three (x, f(x)) pairs defining a bracketed minimum.

    Returns
    -------
    list of tuples
        List of refined minima as (x_min, f_min) tuples, where x_min is the
        precise location of the minimum and f_min is the function value at that point.

    Notes
    -----
    Uses scipy.optimize.minimize_scalar with the 'bracket' method for each identified
    bracket. If optimization fails for a particular bracket, an error message is printed
    and that minimum is skipped.

    Examples
    --------
    >>> def f(x): return (x-2)**2 + 1
    >>> x = np.linspace(0, 4, 100)
    >>> fx = f(x)
    >>> brackets = bracket_minima(x, fx)
    >>> minima = find_bracketed_minima(f, brackets)
    >>> minima[0][0]  # Should be close to 2.0
    """
    l_mins = []

    for ibr, br in enumerate(l_braks):
        (xm1, fm1), (x0, f0), (xp1, fp1) = br
        assert fm1 > f0 < fp1, f"Invalid bracket: f0 is not less than neighbors for bracket {ibr}: {br}."

        res = sciopt.minimize_scalar(minfunc, bracket=(xm1, x0, xp1))
        if res.success:
            #print(f"Found minimum in bracket {xm1, xp1}: {res.x}, {res.fun}")
            #if not (xm1 <= res.x <= xp1):
             #   print(f"Minimum {res.x} found outside bracket {xm1, xp1}.", res)
            l_mins.append((res.x, res.fun))
        else:
            print(f"Failed to find minimum in bracket {xm1, xp1}.")
            print(f"  Message: {res.message}")

    return l_mins
