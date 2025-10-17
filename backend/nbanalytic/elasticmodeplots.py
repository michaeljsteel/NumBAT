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
from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes as MplAxes
from matplotlib.figure import Figure

import numpy as np
from numpy.typing import NDArray
from nbtypes import SI_um

twopi = 2 * math.pi

import IPython.display as ipydisp


def _field_plot_1d(v_x: Sequence[float], m_u: NDArray[np.complexfloating], comps: Sequence[str] = (),
                   labx: str = '', laby: str = '', ax: Optional[MplAxes] = None,
                   legend: bool = False) -> None:
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
    if ax is None:
        fig,ax=plt.subplots(figsize=(5,5))

    cols_comp={'x':'blue', 'y':'red', 'z':'green'}
    ls_reim={'r':'-', 'i':':'}

    if not comps or 'xr' in comps:
        ax.plot(v_x, np.real(m_u[:,0]), c=cols_comp['x'], ls=ls_reim['r'], label='$u_{x,r}$')
    if not comps or 'xi' in comps:
        ax.plot(v_x, np.imag(m_u[:,0]), c=cols_comp['x'], ls=ls_reim['i'], label='$u_{x,i}$')
    if not comps or 'yr' in comps:
        ax.plot(v_x, np.real(m_u[:,1]), c=cols_comp['y'], ls=ls_reim['r'], label='$u_{y,r}$')
    if not comps or 'yi' in comps:
        ax.plot(v_x, np.imag(m_u[:,1]), c=cols_comp['y'], ls=ls_reim['i'], label='$u_{y,i}$')
    if not comps or 'zr' in comps:
        ax.plot(v_x, np.real(m_u[:,2]), c=cols_comp['z'], ls=ls_reim['r'], label='$u_{z,r}$')
    if not comps or 'zi' in comps:
        ax.plot(v_x, np.imag(m_u[:,2]), c=cols_comp['z'], ls=ls_reim['i'], label='$u_{z,i}$')

    if labx: ax.set_xlabel(labx)
    if laby: ax.set_ylabel(laby)

    if legend: ax.legend()



def _field_plot_2d(v_y: Sequence[float], v_z: Sequence[float], m_uxyz: NDArray[np.complexfloating],
                   q: float, phi: float = 0.0, labx: str = '', laby: str = '',
                   ax: Optional[MplAxes] = None, displacement_scale: float = 0.05,
                   use_arrows: bool = False) -> None:
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
    use_arrows : bool
        If True, draw quiver arrows; otherwise show displaced scatter points coloured by amplitude.
    """
    if ax is None:
        fig,ax=plt.subplots(figsize=(5,5))

    if labx: ax.set_xlabel(labx)
    if laby: ax.set_ylabel(laby)


    npts_y = len(v_y)
    npts_z = len(v_z)

    m_Y, m_Z = np.meshgrid(v_y, v_z)

    m_uy = np.full((npts_z, npts_y), m_uxyz[:,1])
    m_uz = np.full((npts_z, npts_y), m_uxyz[:,2])

    m_exp_iqz = np.exp(1j*q*m_Z)
    m_uy *= m_exp_iqz * np.exp(1j*phi)
    m_uz *= m_exp_iqz * np.exp(1j*phi)

    m_abs = np.sqrt(np.abs(m_uy)**2+np.abs(m_uz)**2)
    if use_arrows:
        ax.quiver(v_y, v_z, np.real(m_uy), np.real(m_uz), m_abs)
    else:
        dscal = displacement_scale * (v_y[-1]-v_y[0]) # fraction of the y domain
        ax.scatter(m_Y+dscal*np.real(m_uy), m_Z+dscal*np.real(m_uz), s=.5, c=m_abs)



class SlabMode2DAnimator():
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
    def __init__(self, v_y: Sequence[float], m_uxyz: NDArray[np.complexfloating], q: float,
                 fig: Optional[Figure] = None, ax: Optional[MplAxes] = None, zperiods: int = 1,
                 n_frames: int = 20, interval: int = 100, flip_axes: bool = False,
                 use_arrows: bool = False, displacement_scale: float = 0.05,
                 label: str = '') -> None:

        self.n_frames=n_frames
        self.interval = interval

        if ax is None:
            fig,ax = plt.subplots()

        self.fig = fig
        self.ax = ax

        self.ani=None
        self._use_arrows = use_arrows
        self._displacement_scale = displacement_scale
        self._flip_axes = flip_axes

        npts_y=len(v_y)
        npts_z=npts_y*zperiods
        zlo=0
        zhi=zlo + zperiods * 2*np.pi/q
        v_z = np.linspace(zlo, zhi, npts_z)

        v_Y = v_y/SI_um
        v_Z = v_z/SI_um

        cmp = cm.inferno
        dotsz = 0.5
        ds = self._displacement_scale

        qnorm = q * SI_um

        self.m_uy_phi0 = np.full((npts_z, npts_y), m_uxyz[:,1])
        self.m_uz_phi0 = np.full((npts_z, npts_y), m_uxyz[:,2])
        if flip_axes:
            self.m_uy_phi0 = self.m_uy_phi0.T
            self.m_uz_phi0 = self.m_uz_phi0.T

        if not flip_axes:
            xlab, ylab = r'$y$ [µm]',r'$z$ [µm]'
            m_Y, m_Z = np.meshgrid(v_Y, v_Z)
            m_exp_iqz = np.exp(1j*qnorm*m_Z)

            self.m_uy_phi0 *= m_exp_iqz
            self.m_uz_phi0 *= m_exp_iqz

            if use_arrows:
                self.the_plot = self.ax.quiver(v_Y, v_Z, m_Y*0, m_Z*0, cmap=cmp)
            else:
                self.m_Y=m_Y
                self.m_Z=m_Z
                m_abs = np.sqrt(np.abs(self.m_uy_phi0)**2+ np.abs(self.m_uz_phi0)**2).flatten()
                v_dY = (m_Y + ds * np.real(self.m_uy_phi0)).flatten()
                v_dZ = (m_Z + ds * np.real(self.m_uz_phi0)).flatten()
                self.the_plot = self.ax.scatter(v_dY, v_dZ, c=m_abs, cmap=cmp, s=dotsz)

        else:
            xlab, ylab = r'$z$ [µm]',r'$y$ [µm]'
            m_Z, m_Y = np.meshgrid(v_Z, v_Y)
            m_exp_iqz = np.exp(1j*qnorm*m_Z)

            self.m_uy_phi0 *= m_exp_iqz
            self.m_uz_phi0 *= m_exp_iqz

            if use_arrows:
                self.the_plot = self.ax.quiver(v_Z, v_Y, m_Z*0, m_Y*0, cmap=cmp)
            else:
                self.m_Y=m_Y
                self.m_Z=m_Z
                m_abs = np.sqrt(np.abs(self.m_uy_phi0)**2+ np.abs(self.m_uz_phi0)**2).flatten()
                v_dY = (m_Y + ds * np.real(self.m_uy_phi0)).flatten()
                v_dZ = (m_Z + ds * np.real(self.m_uz_phi0)).flatten()
                self.the_plot = self.ax.scatter(v_dZ, v_dY, c=m_abs, cmap=cmp, s=dotsz)

        #self.ax.set_aspect('equal')

        self.ax.set_xlabel(xlab)
        self.ax.set_ylabel(ylab)
        if label:
            self.ax.set_title(label)

        self._apply_frame_for_phase(0.0)


    def _apply_frame_for_phase(self, phi: float) -> None:
        m_uy_r = np.real(self.m_uy_phi0 * np.exp(-1j*phi))
        m_uz_r = np.real(self.m_uz_phi0 * np.exp(-1j*phi))

        m_abs = np.sqrt(m_uy_r**2+ m_uz_r**2).flatten()
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
            #self.the_plot.set_array(m_abs)  # update colors

    def init_func(self) -> Tuple[object, ...]:
        self._apply_frame_for_phase(0.0)
        return self.the_plot,

    def update_func(self, iphi: int) -> Tuple[object, ...]:
        phi = iphi/self.n_frames * 2*np.pi
        self._apply_frame_for_phase(phi)
        return self.the_plot,

    def build_animation(self):
        self.ani = FuncAnimation(self.fig, self.update_func, frames=self.n_frames,
                            init_func=self.init_func, blit=True,
                                 interval=self.interval, repeat=True)
        return self.ani, self.fig

    def html_animate(self):
        anim, fig = self.build_animation()
        plt.close(self.fig)

        html_anim = ipydisp.HTML(self.ani.to_jshtml())
        ipydisp.display(html_anim)
