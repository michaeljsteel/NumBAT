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

from typing import Optional, Tuple
from numpy.typing import NDArray


import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes as MplAxes
from matplotlib.figure import Figure

from nbtypes import SI_um, SI_GHz, SI_kmps
from enum import IntEnum

import IPython.display as ipydisp
import plotting.plottools as nbpt
import plotting.plotprefs as nbpp

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

    def name_abs_sym_units(self) -> str:
        return f"{self.name} |{self.symbol}| [{self.units}]"

    def sym_units(self) -> str:
        return f"{self.abs_symbol} [{self.units}]"


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
        v_x: NDArray,  # Nelt 1D array of x coordinates
        m_uxyz: NDArray,  # Nelts x 3 array of some scalar complex values
        Vbulks: Optional[NDArray] = None,
        use_scaled_units=True,
        relative_sing_val: Optional[float] = None,
    ) -> None:

        self.Omega_SI = Omega_SI
        self.Vp_SI = Vp_SI
        self.q_SI = Omega_SI / Vp_SI

        self.Omega_norm = self.Omega_SI
        self.Vp_norm = self.Vp_SI
        self.q_norm = self.q_SI

        self.v_x = v_x.copy()  # Nelts
        self.m_uxyz = m_uxyz.copy()  # Nelts x 3
        self.Vbulks = Vbulks
        self.vars_are_scaled = False

        self.mode_id = ''
        self.extra_lab = ''

        self.ivar_x = QuantLabel("Position", r"$y$", r"$|y|$", "µm")
        self.ivar_y = QuantLabel("Position", r"$z$", r"$|z|$", "µm")
        self.dvar = QuantLabel("Displacement", r"$\vec u$", r"$|\vec u|$", "µm")


        gpp = nbpp.PlotPrefs()
        self.fs_suptitle = gpp.d_multi['title_fs'] -3
        self.fs_axtitle = gpp.d_multi['title_fs'] -4
        self.fs_axlab = gpp.d_multi['ax_label_fs']
        self.fs_ticklab = gpp.d_multi['ax_ticklabel_fs']
        self.fs_cbar = gpp.d_multi['cb_label_fs']
        self.cm_cbar = gpp.cmap_ac_field_signed


        self.show_singval = False
        self.relative_sing_val = relative_sing_val

        self.bc_evals=None

        if use_scaled_units:
            self.switch_to_scaled_units()

    def __str__(self) -> str:
        s = f"ModeFunction1D: Omega={self.Omega_SI:.4e} rad/s, Vp={self.Vp_SI:.4f} m/s, npts={len(self.v_x)}, vx={self.v_x[0]} ... {self.v_x[-1]}\n"
        return s

    def interpolate_to_new_ivar_grid(self, npts: int) -> None:
        """Interpolate the mode function onto a new grid with npts points."""
        from scipy.interpolate import interp1d

        new_v_x = np.linspace(self.v_x[0], self.v_x[-1], npts)

        interp_real_f = interp1d(self.v_x, np.real(self.m_uxyz), axis=0, kind="cubic")
        interp_imag_f = interp1d(self.v_x, np.imag(self.m_uxyz), axis=0, kind="cubic")

        self.m_uxyz = interp_real_f(new_v_x) + 1j * interp_imag_f(new_v_x)
        self.v_x = new_v_x.copy()


    def set_Vbulks(self, Vbulks: NDArray) -> None:
        """Set bulk velocities."""
        self.Vbulks = Vbulks

    def set_mode_id(self, mode_id: str) -> None:
        """Set the mode ID."""
        self.mode_id = mode_id

    def switch_to_scaled_units(self) -> None:
        """Set independent variable to microns."""
        if not self.vars_are_scaled:
            self.v_x /= SI_um
            self.vars_are_scaled = True

            self.Omega_norm /= SI_GHz
            self.Vp_norm /= SI_kmps
            self.q_norm *= SI_um

    def set_labels(
        self,
        ivar_x: Optional[QuantLabel] = None,
        ivar_y: Optional[QuantLabel] = None,
        dvar: Optional[QuantLabel] = None,
    ) -> None:
        """Set labels for plotting."""
        if ivar_x is not None:
            self.ivar_x = ivar_x
        if ivar_y is not None:
            self.ivar_y = ivar_y
        if dvar is not None:
            self.dvar = dvar

    def normalise_mode(
        self,
        norm_type: ModeNormalisationType = ModeNormalisationType.MAX_AMPLITUDE,
        value: float = 1.0,
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
                raise NotImplementedError(
                    "TOTAL_ENERGY normalisation not implemented yet."
                )

    def plot_mode_profile_1d(
        self,
        ax: Optional[MplAxes] = None,
        prefix: str = "",
        comps: Tuple[str, ...] = (),
        legend: bool = True,
    ) -> str:
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

        title = make_title(self.Omega_norm, self.Vp_norm, self.mode_id, self.Vbulks, self.extra_lab)

        ax_in = ax
        if ax is None:
            fig, axs = plt.subplots(2, 1, figsize=(6, 8))
            fig.set_tight_layout(True)
            ax_comp, ax_abs = axs
            fig.suptitle(title, fontsize=self.fs_suptitle)
        else:
            ax_comp, ax_abs = ax, None
            axs = [ax_comp]

        cols_comp = {"x": "red", "y": "green", "z": "blue"}
        ls_reim = {"r": "-", "i": ":"}

        d_comps = {"x": 0, "y": 1, "z": 2}

        for c, ic in d_comps.items():
            if comps and ic not in comps:
                continue
            ax_comp.plot(
                self.v_x,
                np.real(self.m_uxyz[:, ic]),
                c=cols_comp[c],
                ls=ls_reim["r"],
                label=f"{self.dvar.symbol}$_{c}$"
            )
            ax_comp.plot(
                self.v_x, np.imag(self.m_uxyz[:, ic]), c=cols_comp[c], ls=ls_reim["i"]
            )

            if ax_abs:
                ax_abs.plot(
                    self.v_x,
                    np.abs(self.m_uxyz[:, ic]),
                    c=cols_comp[c],
                    ls=ls_reim["r"],
                    label=f"$|${self.dvar.symbol}$_{c}|$",
                )

        for ax in axs:
            ax.set_xlabel(self.ivar_x.name_sym_units())

        ax_comp.set_ylabel(f"{self.dvar.name_sym_units()}")
        if ax_abs:
            ax_abs.set_ylabel(f"{self.dvar.name_abs_sym_units()}")


        if self.show_singval and self.relative_sing_val is not None:
            ax_comp.text(
                0.05,
                0.95,
                f"Sing. val.: {self.relative_sing_val:.3e}", fontsize=self.fs_axlab,
                transform=ax_comp.transAxes,
            )
        if self.bc_evals is not None:
            ax_comp.text( 0.05, 0.90,
                f"lam 1.: {np.real(self.bc_evals[0]):+.5f}+{np.imag(self.bc_evals[0]):+.5f}j",
                fontsize=self.fs_axlab, transform=ax_comp.transAxes,)
            ax_comp.text( 0.05, 0.85,
                f"lam 2.: {np.real(self.bc_evals[1]):+.5f}+{np.imag(self.bc_evals[1]):+.5f}j",
                fontsize=self.fs_axlab, transform=ax_comp.transAxes,)
            ax_comp.text( 0.05, 0.80,
                f"lam 3.: {np.real(self.bc_evals[2]):+.5f}+{np.imag(self.bc_evals[2]):+.5f}j",
                fontsize=self.fs_axlab, transform=ax_comp.transAxes,)


        if legend:
            for ax in axs:
                ax.legend(loc='lower left')

        fname = ""
        if ax_in is None:
            fname = prefix + "_profile_1d.png"
            nbpt.save_and_close_figure(fig, fname)

        return fname

    def plot_mode_profile_2d(
        self,
        ax: Optional[MplAxes] = None,
        prefix: str = "",
        zperiods: int = 1,
        npts_per_period: int = 30,
        displacement_scale: float = 0.3,
        style: str = "comps",
    ) -> str:
        """Plot a 2D representation of the modal displacement field.

        Parameters
        ----------
        ax : matplotlib Axes, optional
            Axes to plot into. If None, a new figure is created and saved.
        prefix : str
            Filename prefix for saved figure when ax is None.
        zperiods : int
            Number of axial periods to cover in the plot domain.
        npts_per_period : int
            Sample points per period along the axial direction.
        displacement_scale : float
            Visual displacement scaling factor.
        style : {"comps","morph","quiver"}
            Plot style: component panels, morph scatter, or quiver arrows.

        Returns
        -------
        str
            Path to saved figure (empty string if ax provided).
        """
        kw = locals().copy()
        kw.pop("self")
        kw.pop("style")
        if style == "morph":
            return self._plot_mode_profile_2d_morph(**kw)
        elif style == "comps":
            kw.pop("displacement_scale")
            return self._plot_mode_profile_2d_comps(**kw)
        elif style == "quiver":
            kw.pop("displacement_scale")
            return self._plot_mode_profile_2d_quiver(**kw)
        else:
            raise ValueError(f"Unknown style '{style}' for 2D mode profile plotting.")

    def animate_mode_profile_2d(
        self,
        ax: Optional[MplAxes] = None,
        fname: str = "",
        zperiods: int = 1,
        npts_per_period: int = 30,
        displacement_scale: float = 0.3,
        hide_decorations: bool = False,
        shear_signed: bool = False,
    ) -> str | None:
        """Create and save a simple 2D slab mode animation as HTML.

        Returns the output HTML filename. Currently uses the 'morph' style scatter
        animation over two axial periods.
        """
        anim = SlabMode2DAnimator(
            self.v_x * SI_um,
            self.m_uxyz,
            self.q_norm / SI_um,
            zperiods=2,
            displacement_scale=displacement_scale,
            hide_decorations=hide_decorations,
            shear_signed=shear_signed,
        )
        return anim.html_animate(fname=fname)


    def _make_2d_field_arrays(self, zperiods: int, npts_per_period: int):

        v_z = np.linspace(0, 2 * np.pi / self.q_norm * zperiods, zperiods * npts_per_period)
        m_Y, m_Z = np.meshgrid(self.v_x, v_z)
        m_eiphiz = np.exp(1j * self.q_norm * m_Z).T

        m_Fx = np.real(self.m_uxyz[:, 0][:, np.newaxis] * m_eiphiz)
        m_Fy = np.real(self.m_uxyz[:, 1][:, np.newaxis] * m_eiphiz)
        m_Fz = np.real(self.m_uxyz[:, 2][:, np.newaxis] * m_eiphiz)

        return v_z, m_Y, m_Z, m_Fx, m_Fy, m_Fz

    def _plot_mode_profile_2d_morph(
        self,
        ax: Optional[MplAxes],
        prefix: str,
        zperiods: int,
        npts_per_period: int,
        displacement_scale: float,
    ) -> str:
        """Plot 2D displacement field as either a quiver or scatter field.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw into. If None a new figure and axes are created.
        displacement_scale : float
            Scaling factor applied to displacements for visualization.
        """

        title = make_title(self.Omega_norm, self.Vp_norm, self.mode_id, self.Vbulks, self.extra_lab)

        v_z, m_Y, m_Z, m_Fx, m_Fy, m_Fz = self._make_2d_field_arrays(zperiods, npts_per_period)

        ax_in = ax
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 5))
            fig.set_tight_layout(True)

        if fig:
            fig.suptitle(title, fontsize=self.fs_suptitle)

        ax.set_title(self.dvar.name_sym_units(), fontsize=self.fs_axtitle)
        self.dress_axes([ax], 'ivx', 'ivy')

        m_Fabs = np.sqrt(np.abs(m_Fx) ** 2 + np.abs(m_Fy) ** 2 + np.abs(m_Fz) ** 2)
        dscal = displacement_scale * (self.v_x[-1] - self.v_x[0])  # fraction of the y domain
        ax.scatter(m_Y + dscal * m_Fy.T, m_Z + dscal * m_Fz.T, s=0.5, c=m_Fabs.T,
                   cmap=self.cm_cbar)  # transposes are confusing here!

        fname = ''
        if ax_in is None:
            fname = prefix + "_profile_2d_morph.png"
            nbpt.save_and_close_figure(fig, fname)
        return fname


    def _plot_mode_profile_2d_comps(
        self,
        zperiods: int ,
        npts_per_period: int ,
        ax: Optional[MplAxes] ,
        prefix: str ,
        aspect: float =1.0,
        figsize: Tuple[float, float] = (6, 4),
    ) -> str:
        """Plot 2D mode profile."""

        title = make_title(self.Omega_norm, self.Vp_norm, self.mode_id, self.Vbulks, self.extra_lab)
        v_z, m_Y, m_Z, m_Fx, m_Fy, m_Fz = self._make_2d_field_arrays(zperiods, npts_per_period)


        exts = (self.v_x[0], self.v_x[-1], v_z[0], v_z[-1])
        cm = self.cm_cbar
        kwimgs = dict(
            extent=exts, aspect="auto", origin="lower", cmap=cm  # aspect=aspect,
        )

        fig, axs = plt.subplots(3, 1, figsize=figsize)
        fig.set_tight_layout(True)
        fig.suptitle(title, fontsize=self.fs_suptitle)

        for iax, (m_T, mlab) in enumerate(
            (
                [m_Fx, rf"{self.dvar.symbol}$_x$"],
                [m_Fy, rf"{self.dvar.symbol}$_y$"],
                [m_Fz, rf"{self.dvar.symbol}$_z$"],
            )
        ):
            im = axs[iax].imshow(m_T.T, **kwimgs)
            axs[iax].set_title(mlab, fontsize=self.fs_axlab+2)

            cbar = fig.colorbar(im, ax=axs[iax])
            cbar.ax.tick_params(labelsize=self.fs_ticklab)
            cbar.set_label(rf"[{self.dvar.units}]", fontsize=self.fs_cbar)




        self.dress_axes(axs, 'ivx', 'ivy')

        fname = prefix + "_profile_2d_comps.png"
        nbpt.save_and_close_figure(fig, fname)
        return fname

    def _plot_mode_profile_2d_quiver(
        self,
        ax: Optional[MplAxes],
        prefix: str ,
        zperiods: int ,
        npts_per_period: int ,
        aspect: float = 1.0,
        figsize: Tuple[float, float] = (6, 4),
    ) -> str:
        """Plot 2D mode profile."""

        title = make_title(self.Omega_norm, self.Vp_norm, self.mode_id, self.Vbulks, self.extra_lab)

        v_z, m_Y, m_Z, m_Fx, m_Fy, m_Fz = self._make_2d_field_arrays(zperiods, npts_per_period)

        fig, ax = plt.subplots(figsize=figsize)
        fig.set_tight_layout(True)
        fig.suptitle(title, fontsize=self.fs_suptitle)

        ax.quiver(
            m_Y,
            m_Z,
            m_Fy.T,
            m_Fz.T,
            scale=50,
            width=0.002,
            headwidth=3,
            headlength=4,
            minlength=0.1,
            color="b",
        )
        ax.set_title(self.dvar.name_sym_units(), fontsize=self.fs_axtitle)

        self.dress_axes([ax], 'ivx', 'ivy')


        fname = prefix + "_profile_2d_vec.png"
        nbpt.save_and_close_figure(fig, fname)
        return fname

    def dress_axes(self, axs, vx: str, vy: str, vz: str='') -> None:
        d_ql = {'ivx': self.ivar_x, 'ivy': self.ivar_y, 'dv': self.dvar}
        for tax in axs:
            tax.tick_params(axis="both", labelsize=self.fs_ticklab)
            tax.set_xlabel(d_ql[vx].name_sym_units(), fontsize=self.fs_axlab)
            tax.set_ylabel(d_ql[vy].name_sym_units(), fontsize=self.fs_axlab)
            if vz:
                tax.set_xlabel(d_ql[vx].name_sym_units(), fontsize=self.fs_axlab)


def make_title(
    Omega_norm: float, Vp_norm: float, mode_id: str, Vbulks: Optional[NDArray] = None,
    extralab: str='') -> str:
    nutil = Omega_norm / (2 * np.pi)
    Vstil = Vp_norm
    title = rf"Mode {mode_id} for {extralab} $\Omega/2\pi$={nutil:.3f} GHz, $V$={Vstil:.3f} km/s"

    if Vbulks is not None:
        title += (
            f"\nBulk Velocities: {Vbulks[0]:.3f}, {Vbulks[1]:.3f}, {Vbulks[2]:.3f} km/s"
        )

    return title



class SlabMode2DAnimator:
    """Animate 2D slab-like mode displacements over one temporal period.

    This helper wraps matplotlib's FuncAnimation and produces an animation
    showing the real part of the in-plane displacements as they oscillate in time.

    Parameters
    ----------
    v_y_SI : array-like
        Coordinates along the lateral direction (m).
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
        v_y_SI: NDArray,
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
        hide_decorations: bool = False,
        shear_signed: bool = False,
    ) -> None:

        pp = nbpp.PlotPrefs()


        self.n_frames = n_frames
        self.interval = interval
        self.shear_signed = shear_signed

        if ax is None:
            fig, ax = plt.subplots()

        self.fig = fig
        self.ax = ax

        self.ani = None
        self._use_arrows = use_arrows
        self._displacement_scale = displacement_scale
        self._flip_axes = flip_axes

        npts_y = len(v_y_SI)
        npts_z = npts_y #* zperiods
        zlo = 0
        zhi = zlo + zperiods * 2 * np.pi / q
        v_z = np.linspace(zlo, zhi, npts_z)

        v_Y = v_y_SI / SI_um
        v_Z = v_z/ SI_um
        #cmp = cm.inferno
        if self.shear_signed:
            cmp = pp.cmap_em_field_signed
        else:
            cmp = pp.cmap_em_field_unsigned

        dotsz = 0.5
        ds = self._displacement_scale

        qnorm = q * SI_um

        self.m_ux_phi0 = np.full((npts_z, npts_y), m_uxyz[:, 0])
        self.m_uy_phi0 = np.full((npts_z, npts_y), m_uxyz[:, 1])
        self.m_uz_phi0 = np.full((npts_z, npts_y), m_uxyz[:, 2])

        if self.shear_signed:
            vmin, vmax = -1,1
        else:
            vmin, vmax = 0, None

        if flip_axes:
            self.m_ux_phi0 = self.m_ux_phi0.T
            self.m_uy_phi0 = self.m_uy_phi0.T
            self.m_uz_phi0 = self.m_uz_phi0.T

        if not flip_axes:
            xlab, ylab = r"$y$ [µm]", r"$z$ [µm]"
            m_Y, m_Z = np.meshgrid(v_Y, v_Z)
            m_exp_iqz = np.exp(1j * qnorm * m_Z)


            self.m_ux_phi0 *= m_exp_iqz
            self.m_uy_phi0 *= m_exp_iqz
            self.m_uz_phi0 *= m_exp_iqz

            if use_arrows:
                self.the_plot = self.ax.quiver(v_Y, v_Z, m_Y * 0, m_Z * 0, cmap=cmp)
            else:
                self.m_Y = m_Y
                self.m_Z = m_Z
                m_abs = np.sqrt(
                    np.abs(self.m_ux_phi0) ** 2 +
                    np.abs(self.m_uy_phi0) ** 2 + np.abs(self.m_uz_phi0) ** 2
                ).flatten()
                v_dY = (m_Y + ds * np.real(self.m_uy_phi0)).flatten()
                v_dZ = (m_Z + ds * np.real(self.m_uz_phi0)).flatten()
                self.the_plot = self.ax.scatter(v_dY, v_dZ, c=m_abs, cmap=cmp, s=dotsz, vmin=vmin, vmax=vmax)

        else:
            xlab, ylab = r"$z$ [µm]", r"$y$ [µm]"
            m_Z, m_Y = np.meshgrid(v_Z, v_Y)
            m_exp_iqz = np.exp(1j * qnorm * m_Z)

            self.m_ux_phi0 *= m_exp_iqz
            self.m_uy_phi0 *= m_exp_iqz
            self.m_uz_phi0 *= m_exp_iqz

            if use_arrows:
                self.the_plot = self.ax.quiver(v_Z, v_Y, m_Z * 0, m_Y * 0, cmap=cmp)
            else:
                self.m_Y = m_Y
                self.m_Z = m_Z
                m_abs = np.sqrt(
                    np.abs(self.m_ux_phi0) ** 2 +
                    np.abs(self.m_uy_phi0) ** 2 + np.abs(self.m_uz_phi0) ** 2
                ).flatten()
                v_dY = (m_Y + ds * np.real(self.m_uy_phi0)).flatten()
                v_dZ = (m_Z + ds * np.real(self.m_uz_phi0)).flatten()
                self.the_plot = self.ax.scatter(v_dZ, v_dY, c=m_abs, cmap=cmp, s=dotsz, vmin=vmin, vmax=vmax)

        # self.ax.set_aspect('equal')

        if hide_decorations:
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.ax.set_xticklabels([])
            self.ax.set_yticklabels([])
            for spine in self.ax.spines.values():
                spine.set_visible(False)

        else:
            self.ax.set_xlabel(xlab)
            self.ax.set_ylabel(ylab)
            if label:
                self.ax.set_title(label)


        self._apply_frame_for_phase(0.0)

    def _apply_frame_for_phase(self, phi: float) -> None:
        m_ux_r = np.real(self.m_ux_phi0 * np.exp(-1j * phi))
        m_uy_r = np.real(self.m_uy_phi0 * np.exp(-1j * phi))
        m_uz_r = np.real(self.m_uz_phi0 * np.exp(-1j * phi))

        m_abs = np.sqrt(m_ux_r**2 + m_uy_r**2 + m_uz_r**2).flatten()
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
            if self.shear_signed:
                m_abs = m_ux_r.flatten()
            else :
                m_abs = np.sqrt(m_uy_r**2 + m_uz_r**2).flatten()
            self.the_plot.set_array(m_abs)  # update colors

    def init_func(self) -> Tuple[object, ...]:
        self._apply_frame_for_phase(0.0)
        return (self.the_plot,)

    def update_func(self, iphi: int) -> Tuple[object, ...]:
        phi = iphi / self.n_frames * 2 * np.pi
        self._apply_frame_for_phase(phi)
        return (self.the_plot,)

    def build_animation(self):
        self.ani = FuncAnimation(
            self.fig,                  # type: ignore
            self.update_func,          # type: ignore
            frames=self.n_frames,
            init_func=self.init_func,  # type:ignore
            blit=True,
            interval=self.interval,
            repeat=True,
        )
        return self.ani, self.fig

    def html_animate(self, fname: Optional[str] = None) -> Optional[str]:
        """Render the animation as HTML/JS. If ``output_path`` is provided, write to file.

        Parameters
        ----------
        output_path : str, optional
            Path to write the HTML fragment. If omitted, display inline.

        Returns
        -------
        Optional[str]
            The path written if ``output_path`` was supplied, else None.
        """
        anim, fig = self.build_animation()
        plt.close(self.fig)

        if not fname:
            anim_out = self.ani.to_jshtml()  # type: ignore
            html_anim = ipydisp.HTML(anim_out)
            ipydisp.display(html_anim)
            return None


        fmode = fname.split('.')[-1]
        anim_out = None
        #fmode = 'jshtml'
        #fmode = 'html'
        #fmode = 'mp4'
        #fmode = 'gif'
        #fname = f"{prefix}_anim.{fmode}"

        if fmode == 'jshtml': #display in pynotebook
            anim_out = self.ani.to_jshtml()  # type: ignore
            with open(fname, "w", encoding="utf-8") as f:
                f.write(anim_out)
        elif fmode == 'html':
            anim_out = self.ani.to_html5_video()  # type: ignore
            with open(fname, "w", encoding="utf-8") as f:
                f.write(anim_out)

        elif fmode == 'mp4':
            self.ani.save(fname, writer='ffmpeg', dpi=120)  # type: ignore
        else:  # gif
            self.ani.save(fname, writer='imagemagick', dpi=80)  # type: ignore

        return fname