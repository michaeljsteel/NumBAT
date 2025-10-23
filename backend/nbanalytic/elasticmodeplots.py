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


"""elasticmodeplots.py
----------------------
Small plotting utilities for visualising analytic elastic mode solutions.

This module provides simple 1D and 2D field plotting helpers and a
small animation helper class `SlabMode2DAnimator` that renders time-dependent
displacements for slab-like mode fields.

Functions and classes are intentionally lightweight and designed for
interactive use in notebooks and small scripts.
"""

import math
import numpy as np

from typing import Optional, Sequence, Tuple
from numpy.typing import NDArray


import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes as MplAxes
from matplotlib.figure import Figure

import numbat
from nbtypes import SI_um, SI_GHz, SI_kmps
from enum import IntEnum

import IPython.display as ipydisp
import plotting.plottools as nbpt

twopi = 2 * math.pi


class QuantLabel:
    def __init__(self, name: str, symbol: str, abs_symbol: str, units: str) -> None:
        self.name = name
        self.symbol = symbol
        self.abs_symbol = abs_symbol
        self.units = units

    def name_sym(self) -> str:
        return f"{self.name} {self.symbol}"

    def name_sym_units(self) -> str:
        return f"{self.name} {self.symbol} [{self.units}]"


class ModeNormalisationType(IntEnum):
    """Types of normalisation for mode functions."""

    MAX_AMPLITUDE = 1
    MAX_MAGNITUDE = 2
    TOTAL_ENERGY = 3


class ModeFunction1D:
    """Wrapper for 1D mode functions."""

    def __init__(
        self,
        Omega_SI: float,
        Vp_SI: float,
        v_x: NDArray,
        m_uxyz: NDArray,
        Vbulks: Optional[NDArray] = None,
        use_SI_dimensions=False,
    ) -> None:

        self.Omega_SI = Omega_SI
        self.Vp_SI = Vp_SI
        self.q_SI = Omega_SI / Vp_SI

        self.Omega_norm = Omega_SI
        self.Vp_norm = Vp_SI
        self.q_norm = Omega_SI / Vp_SI

        self.v_x = v_x.copy()  # Nelts
        self.m_uxyz = m_uxyz.copy()  # Nelts x 3
        self.Vbulks = Vbulks
        self.indep_var_is_SI = True
        # self.indep_var_is_y = True # v_x is y coordinate
        self.mode_id = ""

        self.set_labels()

        if not use_SI_dimensions:
            self.set_indep_var_to_microns()

    def __str__(self) -> str:
        s = f"ModeFunction1D: Omega={self.Omega_SI:.4e} rad/s, Vp={self.Vp_SI:.4f} m/s, npts={len(self.v_x)}, vx={self.v_x[0]} ... {self.v_x[-1]}\n"
        return s

    def set_mode_id(self, mode_id: int) -> None:
        """Set the mode ID."""
        self.mode_id = mode_id

    def set_indep_var_to_microns(self) -> None:
        """Set independent variable to microns."""
        if self.indep_var_is_SI:
            self.v_x /= SI_um
            self.indep_var_is_SI = False

            self.Omega_norm /= SI_GHz
            self.Vp_norm /= SI_kmps
            self.q_norm *= SI_um

    def set_labels(
        self,
        ivar_x=QuantLabel("Position", r"$y$", r"$|y|$", "µm"),
        ivar_y=QuantLabel("Position", r"$z$", r"$|z|$", "µm"),
        dvar=QuantLabel("Displacement", r"$\vec u$", r"$|\vec u|$", "µm"),
    ) -> None:
        """Set labels for plotting."""

        self.ivar_x = ivar_x
        self.ivar_y = ivar_y
        self.dvar = dvar

    def normalise_mode(
        self, norm_type: ModeNormalisationType, value: float = 1.0
    ) -> None:

        m_uxyz = self.m_uxyz

        match norm_type:
            case ModeNormalisationType.MAX_AMPLITUDE:
                maxel = m_uxyz.flat[np.abs(m_uxyz).argmax()]
                m_uxyz *= value / maxel

            case ModeNormalisationType.MAX_MAGNITUDE:
                maxmag = np.max(np.abs(m_uxyz))
                m_uxyz *= value / maxmag

            case (
                ModeNormalisationType.TOTAL_ENERGY
            ):  # should really use int of S.T, or  Omega^2 rho |u|^2
                pass

    def plot_mode_profile_1d(
        self,
        comps: Tuple[str, ...] = (),
        Vbulks: Optional[NDArray] = None,
        ax: Optional[MplAxes] = None,
        legend: bool = True,
        prefix: str = "",
    ) -> None:
        """Plot 1D mode profile."""

        nutil = self.Omega_norm / (2 * np.pi)
        Vstil = self.Vp_norm
        title = rf"Mode {self.mode_id} for $\Omega/2\pi$={nutil:.3f} GHz, $V$={Vstil:.3f} km/s"

        if Vbulks is not None:
            title += f"\nBulk Velocities: {Vbulks[0]:.3f}, {Vbulks[1]:.3f}, {Vbulks[2]:.3f} km/s"

        return _field_plot_1d(
            self.v_x,
            self.m_uxyz,
            comps=comps,
            title=title,
            ivar_x=self.ivar_x,
            dvar=self.dvar,
            ax=ax,
            legend=legend,
            prefix=prefix,
        )

    def plot_mode_profile_2d(
        self,
        zperiods: int = 1,
        npts_per_period: int = 15,
        displacement_scale: float = 0.3,
        Vbulks: Optional[NDArray] = None,
        phi: float = 0.0,
        ax: Optional[MplAxes] = None,
        prefix: str = "",
        style: str = "comps",
    ) -> None:
        """Plot 2D mode profile."""

        kw = locals()
        kw.pop("self")
        kw.pop("style")
        if style == "morph":
            return self.plot_mode_profile_2d_morph(**kw)

        elif style == "comps":
            kw.pop("displacement_scale")
            kw.pop("phi")
            return self.plot_mode_profile_2d_comps(**kw)

        elif style == "quiver":
            kw.pop("displacement_scale")
            kw.pop("phi")
            return self.plot_mode_profile_2d_quiver(**kw)

        else:
            raise ValueError(f"Unknown style '{style}' for 2D mode profile plotting.")


    def plot_mode_profile_2d_morph(
        self,
        zperiods: int = 1,
        npts_per_period: int = 15,
        displacement_scale: float = 0.05,
        Vbulks: Optional[NDArray] = None,
        phi: float = 0.0,
        ax: Optional[MplAxes] = None,
        prefix: str = "",
    ) -> None:
        """Plot 2D mode profile."""

        title = make_title(self.Omega_norm, self.Vp_norm, self.mode_id, Vbulks)

        npts_z = npts_per_period * zperiods
        zlo = 0
        zhi = zlo + zperiods * 2 * np.pi / self.q_norm
        v_z = np.linspace(zlo, zhi, npts_z)
        q = self.Omega_norm / self.Vp_norm

        return _field_plot_2d(
            self.v_x,
            v_z,
            self.m_uxyz,
            q=q,
            phi=phi,
            ivar_x=self.ivar_x,
            ivar_y=self.ivar_y,
            dvar=self.dvar,
            title=title,
            ax=ax,
            displacement_scale=displacement_scale,
            prefix=prefix,
        )

    def plot_mode_profile_2d_comps(
        self,
        zperiods: int = 1,
        npts_per_period: int = 15,
        Vbulks: Optional[NDArray] = None,
        figsize: Tuple[float, float] = (6, 4),
        ax: Optional[MplAxes] = None,
        aspect: float = 1.0,
        prefix: str = "",
    ) -> None:
        """Plot 2D mode profile."""

        title = make_title(self.Omega_norm, self.Vp_norm, self.mode_id, Vbulks)

        q = self.Omega_norm / self.Vp_norm
        v_z = np.linspace(0, 2 * np.pi / q * zperiods, zperiods * npts_per_period)
        v_Y, v_Z = np.meshgrid(self.v_x, v_z)

        m_eiphiz = np.exp(1j * q * v_Z).T

        m_fx = np.real(self.m_uxyz[:, 0][:, np.newaxis] * m_eiphiz)
        m_fy = np.real(self.m_uxyz[:, 1][:, np.newaxis] * m_eiphiz)
        m_fz = np.real(self.m_uxyz[:, 2][:, np.newaxis] * m_eiphiz)

        exts = (self.v_x[0], self.v_x[-1], v_z[0], v_z[-1])
        fs = 6
        plotprefs = numbat.NumBATPlotPrefs()
        cm = plotprefs.cmap_ac_field_signed
        kwimgs = dict(
            extent=exts, aspect="auto", origin="lower", cmap=cm  # aspect=aspect,
        )

        fig, axs = plt.subplots(3, 1, figsize=figsize)
        fig.set_tight_layout(True)
        fig.suptitle(title, fontsize=8)

        for iax, (m_T, mlab) in enumerate(
            (
                [m_fx, rf"{self.dvar.symbol}$_x]$"],
                [m_fy, rf"{self.dvar.symbol}$_y]$"],
                [m_fz, rf"{self.dvar.symbol}$_z]$"],
            )
        ):
            im = axs[iax].imshow(m_T.T, **kwimgs)
            axs[iax].set_title(mlab, fontsize=fs)

            cbar = fig.colorbar(im, ax=axs[iax])
            cbar.ax.tick_params(labelsize=fs)
            cbar.set_label(rf"[{self.dvar.units}]", fontsize=fs)

        for tax in axs.flatten():
            tax.tick_params(axis="both", labelsize=fs)
            tax.set_xlabel(self.ivar_x.name_sym_units(), fontsize=fs)
            tax.set_ylabel(self.ivar_y.name_sym_units(), fontsize=fs)

        fname = prefix + "_profile_2d_comps.png"
        nbpt.save_and_close_figure(fig, fname)
        return fname

    def plot_mode_profile_2d_quiver(
        self,
        zperiods: int = 1,
        npts_per_period: int = 15,
        Vbulks: Optional[NDArray] = None,
        figsize: Tuple[float, float] = (6, 4),
        ax: Optional[MplAxes] = None,
        aspect: float = 1.0,
        prefix: str = "",
    ) -> None:
        """Plot 2D mode profile."""

        title = make_title(self.Omega_norm, self.Vp_norm, self.mode_id, Vbulks)

        q = self.Omega_norm / self.Vp_norm
        v_z = np.linspace(0, 2 * np.pi / q * zperiods, zperiods * npts_per_period)
        v_Y, v_Z = np.meshgrid(self.v_x, v_z)

        m_eiphiz = np.exp(1j * q * v_Z).T

        # m_fx = np.real(self.m_uxyz[:, 0][:, np.newaxis] * m_eiphiz)
        m_fy = np.real(self.m_uxyz[:, 1][:, np.newaxis] * m_eiphiz)
        m_fz = np.real(self.m_uxyz[:, 2][:, np.newaxis] * m_eiphiz)

        fs = 6
        # plotprefs = numbat.NumBATPlotPrefs()

        fig, ax = plt.subplots(figsize=figsize)
        fig.set_tight_layout(True)
        fig.suptitle(title, fontsize=8)

        ax.quiver(
            v_Y[::2, ::2],
            v_Z[::2, ::2],
            m_fy[::2, ::2].T,
            m_fz[::2, ::2].T,
            # scale=1/displacement_scale,
            # angles='xy',
            # scale_units='xy',
            width=0.002,
            headwidth=3,
            headlength=4,
            minlength=0.1,
            color="b",
        )
        ax.set_title(self.dvar.name_sym_units(), fontsize=fs)

        for tax in [ax]:
            tax.tick_params(axis="both", labelsize=fs)
            tax.set_xlabel(self.ivar_x.name_sym_units(), fontsize=fs)
            tax.set_ylabel(self.ivar_y.name_sym_units(), fontsize=fs)

        fname = prefix + "_profile_2d_vec.png"
        nbpt.save_and_close_figure(fig, fname)
        return fname


def make_title(
    Omega_norm: str, Vp_norm: str, mode_id: str, Vbulks: Optional[NDArray] = None
) -> str:
    nutil = Omega_norm / (2 * np.pi)
    Vstil = Vp_norm
    title = rf"Mode {mode_id} for $\Omega/2\pi$={nutil:.3f} GHz, $V$={Vstil:.3f} km/s"

    if Vbulks is not None:
        title += (
            f"\nBulk Velocities: {Vbulks[0]:.3f}, {Vbulks[1]:.3f}, {Vbulks[2]:.3f} km/s"
        )

    return title


def _field_plot_1d(
    v_x: Sequence[float],
    m_u: NDArray[np.complexfloating],
    comps: Sequence[str] = (),
    title: str = "",
    ivar_x=QuantLabel("Position", r"$y$", r"$|y|$", "[µm]"),
    dvar=QuantLabel("Displacement", r"$\vec u(y)$", r"$|\vec u(y)|$", "[µm]"),
    ax: Optional[MplAxes] = None,
    legend: bool = True,
    prefix: str = "tmpmode",
) -> None:
    """Plot 1D real/imag components of a complex displacement field.

    Parameters
    ----------
    v_x : array-like
        1D coordinate values for the horizontal axis.
    m_u : ndarray
        Nx3 array of complex displacement components [ux, uy, uz].
    comps : sequence, optional
        Which components to plot (e.g. 'xr','yi'). If empty, all components are plotted.
    labx, laby : str, optional
        Axis labels.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. If None a new figure and axes are created.
    legend : bool, optional
        Show a legend when True.
    """

    ax_in = ax
    if ax is None:
        fig, axs = plt.subplots(2, 1, figsize=(6, 8))
        fig.set_tight_layout(True)
        ax_comp, ax_abs = axs
        fig.suptitle(title, fontsize=12)
    else:
        ax_comp, ax_abs = ax, None
        axs = [ax_comp]

    cols_comp = {"x": "red", "y": "green", "z": "blue"}
    ls_reim = {"r": "-", "i": ":"}

    if not comps or "xr" in comps:
        ax_comp.plot(
            v_x, np.real(m_u[:, 0]), c=cols_comp["x"], ls=ls_reim["r"], label="$u_{x}$"
        )
    if not comps or "xi" in comps:
        ax_comp.plot(
            v_x, np.imag(m_u[:, 0]), c=cols_comp["x"], ls=ls_reim["i"]
        )  #  label='$u_{x,i}$')
    if not comps or "yr" in comps:
        ax_comp.plot(
            v_x, np.real(m_u[:, 1]), c=cols_comp["y"], ls=ls_reim["r"], label="$u_{y}$"
        )
    if not comps or "yi" in comps:
        ax_comp.plot(
            v_x, np.imag(m_u[:, 1]), c=cols_comp["y"], ls=ls_reim["i"]
        )  #  label='$u_{y,i}$')
    if not comps or "zr" in comps:
        ax_comp.plot(
            v_x, np.real(m_u[:, 2]), c=cols_comp["z"], ls=ls_reim["r"], label="$u_{z}$"
        )
    if not comps or "zi" in comps:
        ax_comp.plot(
            v_x, np.imag(m_u[:, 2]), c=cols_comp["z"], ls=ls_reim["i"]
        )  # , label='$u_{z,i}$')

    if ax_abs:
        ax_abs.plot(
            v_x, np.abs(m_u[:, 0]), c=cols_comp["x"], ls=ls_reim["r"], label="$|u_{x}|$"
        )
        ax_abs.plot(
            v_x, np.abs(m_u[:, 1]), c=cols_comp["y"], ls=ls_reim["r"], label="$|u_{y}|$"
        )
        ax_abs.plot(
            v_x, np.abs(m_u[:, 2]), c=cols_comp["z"], ls=ls_reim["r"], label="$|u_{z}|$"
        )

    for ax in axs:
        ax.set_xlabel(ivar_x.name_sym_units())

    ax_comp.set_ylabel(f"{dvar.name_sym_units()}")
    if ax_abs:
        ax_abs.set_ylabel(f"{dvar.name_sym_units()}")

    if legend:
        for ax in axs:
            ax.legend()

    fname = ""
    if ax_in is None:
        fname = prefix + "_profile_1d.png"
        fig.savefig(fname)
    return fname


def _field_plot_2d(
    v_y: Sequence[float],
    v_z: Sequence[float],
    m_uxyz: NDArray[np.complexfloating],
    q: float,
    phi: float = 0.0,
    title='',
    ivar_x=QuantLabel("Position", r"$y$", r"$|y|$", "[µm]"),
    ivar_y=QuantLabel("Position", r"$z$", r"$|z|$", "[µm]"),
    dvar=QuantLabel("Displacement", r"$\vec u(y)$", r"$|\vec u(y)|$", "[µm]"),
    ax: Optional[MplAxes] = None,
    displacement_scale: float = 0.05,
    prefix: str = "tmpmode",
) -> None:
    """Plot 2D displacement field as either a quiver or scatter field.

    Parameters
    ----------
    v_y, v_z : array-like
        1D coordinate arrays defining the plotting grid in physical units.
    m_uxyz : ndarray
        3-element complex vector [ux, uy, uz] representing modal profile across the cross-section.
    q : float
        Axial wavenumber used to add phase variation along the z-direction.
    phi : float
        Additional phase to apply to the field.
    labx, laby : str
        Axis labels.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. If None a new figure and axes are created.
    displacement_scale : float
        Scaling factor applied to displacements for visualization.
    """
    ax_in = ax
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
        fig.set_tight_layout(True)

    if fig: fig.suptitle(title, fontsize=10)

    ax.set_title(dvar.name_sym_units())
    ax.set_xlabel(ivar_x.name_sym_units())
    ax.set_ylabel(ivar_y.name_sym_units())

    npts_y = len(v_y)
    npts_z = len(v_z)

    m_Y, m_Z = np.meshgrid(v_y, v_z)

    m_uy = np.full((npts_z, npts_y), m_uxyz[:, 1])
    m_uz = np.full((npts_z, npts_y), m_uxyz[:, 2])

    m_exp_iqZ = np.exp(1j * q * m_Z) * np.exp(1j * phi)
    m_uy = np.real(m_uy * m_exp_iqZ)
    m_uz = np.real(m_uz * m_exp_iqZ)

    m_abs = np.sqrt(np.abs(m_uy) ** 2 + np.abs(m_uz) ** 2)
    dscal = .1*displacement_scale * (v_y[-1] - v_y[0])  # fraction of the y domain
    ax.scatter(m_Y + dscal * m_uy, m_Z + dscal * m_uz, s=0.5, c=m_abs)

    fname = ""
    if ax_in is None:
        fname = prefix + "_profile_2d_morph.png"
        fig.savefig(fname)
    return fname


class SlabMode2DAnimator:
    """Animate 2D slab-like mode displacements over one temporal period.

    This helper wraps matplotlib's FuncAnimation and produces an animation
    showing the real part of the in-plane displacements as they oscillate in time.

    Parameters
    ----------
    v_y : array-like
        Coordinates along the lateral direction (µm).
    m_uxyz : ndarray
        3-element complex displacement vector describing the modal profile.
    q : float
        Axial wavenumber (1/µm).
    fig, ax : matplotlib Figure/Axes, optional
        If provided, animation is drawn into the supplied axes.
    zperiods : int
        Number of z-periods to display across the plotted z-range.
    n_frames : int
        Number of frames in the animation.
    interval : int
        Frame interval in milliseconds.
    flip_axes : bool
        If True, swap y/z axes for plotting convenience.
    use_arrows : bool
        Use quiver arrows instead of scatter points.
    displacement_scale : float
        Scaling factor for visual displacement amplitude.
    label : str
        Optional plot title.
    """

    def __init__(
        self,
        v_y: Sequence[float],
        m_uxyz: NDArray[np.complexfloating],
        q: float,
        fig: Optional[Figure] = None,
        ax: Optional[MplAxes] = None,
        zperiods: int = 1,
        n_frames: int = 20,
        interval: int = 100,
        flip_axes: bool = False,
        use_arrows: bool = False,
        displacement_scale: float = 0.05,
        label: str = "",
    ) -> None:

        self.n_frames = n_frames
        self.interval = interval

        if ax is None:
            fig, ax = plt.subplots()

        self.fig = fig
        self.ax = ax

        self.ani = None
        self._use_arrows = use_arrows
        self._displacement_scale = displacement_scale
        self._flip_axes = flip_axes

        npts_y = len(v_y)
        npts_z = npts_y * zperiods
        zlo = 0
        zhi = zlo + zperiods * 2 * np.pi / q
        v_z = np.linspace(zlo, zhi, npts_z)

        v_Y = v_y / SI_um
        v_Z = v_z / SI_um

        cmp = cm.inferno
        dotsz = 0.5
        ds = self._displacement_scale

        qnorm = q * SI_um

        self.m_uy_phi0 = np.full((npts_z, npts_y), m_uxyz[:, 1])
        self.m_uz_phi0 = np.full((npts_z, npts_y), m_uxyz[:, 2])
        if flip_axes:
            self.m_uy_phi0 = self.m_uy_phi0.T
            self.m_uz_phi0 = self.m_uz_phi0.T

        if not flip_axes:
            xlab, ylab = r"$y$ [µm]", r"$z$ [µm]"
            m_Y, m_Z = np.meshgrid(v_Y, v_Z)
            m_exp_iqz = np.exp(1j * qnorm * m_Z)

            self.m_uy_phi0 *= m_exp_iqz
            self.m_uz_phi0 *= m_exp_iqz

            if use_arrows:
                self.the_plot = self.ax.quiver(v_Y, v_Z, m_Y * 0, m_Z * 0, cmap=cmp)
            else:
                self.m_Y = m_Y
                self.m_Z = m_Z
                m_abs = np.sqrt(
                    np.abs(self.m_uy_phi0) ** 2 + np.abs(self.m_uz_phi0) ** 2
                ).flatten()
                v_dY = (m_Y + ds * np.real(self.m_uy_phi0)).flatten()
                v_dZ = (m_Z + ds * np.real(self.m_uz_phi0)).flatten()
                self.the_plot = self.ax.scatter(v_dY, v_dZ, c=m_abs, cmap=cmp, s=dotsz)

        else:
            xlab, ylab = r"$z$ [µm]", r"$y$ [µm]"
            m_Z, m_Y = np.meshgrid(v_Z, v_Y)
            m_exp_iqz = np.exp(1j * qnorm * m_Z)

            self.m_uy_phi0 *= m_exp_iqz
            self.m_uz_phi0 *= m_exp_iqz

            if use_arrows:
                self.the_plot = self.ax.quiver(v_Z, v_Y, m_Z * 0, m_Y * 0, cmap=cmp)
            else:
                self.m_Y = m_Y
                self.m_Z = m_Z
                m_abs = np.sqrt(
                    np.abs(self.m_uy_phi0) ** 2 + np.abs(self.m_uz_phi0) ** 2
                ).flatten()
                v_dY = (m_Y + ds * np.real(self.m_uy_phi0)).flatten()
                v_dZ = (m_Z + ds * np.real(self.m_uz_phi0)).flatten()
                self.the_plot = self.ax.scatter(v_dZ, v_dY, c=m_abs, cmap=cmp, s=dotsz)

        # self.ax.set_aspect('equal')

        self.ax.set_xlabel(xlab)
        self.ax.set_ylabel(ylab)
        if label:
            self.ax.set_title(label)

        self._apply_frame_for_phase(0.0)

    def _apply_frame_for_phase(self, phi: float) -> None:
        m_uy_r = np.real(self.m_uy_phi0 * np.exp(-1j * phi))
        m_uz_r = np.real(self.m_uz_phi0 * np.exp(-1j * phi))

        m_abs = np.sqrt(m_uy_r**2 + m_uz_r**2).flatten()
        if self._use_arrows:
            if self._flip_axes:
                self.the_plot.set_UVC(m_uz_r, m_uy_r, m_abs)
            else:
                self.the_plot.set_UVC(m_uy_r, m_uz_r, m_abs)
        else:
            v_dY = (self.m_Y + self._displacement_scale * (m_uy_r)).flatten()
            v_dZ = (self.m_Z + self._displacement_scale * (m_uz_r)).flatten()
            if self._flip_axes:
                self.the_plot.set_offsets(np.column_stack([v_dZ, v_dY]))
            else:
                self.the_plot.set_offsets(np.column_stack([v_dY, v_dZ]))
            # self.the_plot.set_array(m_abs)  # update colors

    def init_func(self) -> Tuple[object, ...]:
        self._apply_frame_for_phase(0.0)
        return (self.the_plot,)

    def update_func(self, iphi: int) -> Tuple[object, ...]:
        phi = iphi / self.n_frames * 2 * np.pi
        self._apply_frame_for_phase(phi)
        return (self.the_plot,)

    def build_animation(self):
        self.ani = FuncAnimation(
            self.fig,
            self.update_func,
            frames=self.n_frames,
            init_func=self.init_func,
            blit=True,
            interval=self.interval,
            repeat=True,
        )
        return self.ani, self.fig

    def html_animate(self):
        anim, fig = self.build_animation()
        plt.close(self.fig)

        html_anim = ipydisp.HTML(self.ani.to_jshtml())
        ipydisp.display(html_anim)
