"""Field rendering over NumBAT finite-element meshes."""

import copy
from collections.abc import Sequence
from numbers import Real
from pathlib import Path
from typing import Any, Optional, Union
import itertools

import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np

from numbat_core import runtime
from numbat_core.mesh.fem import FemMesh, validate_1D_cut_format
from numbat_core.plotting import tools as plottools
from numbat_core.utils.numeric import np_min_max


class FEMScalarFieldPlotter:
    """Class to make 1D and 2D plots of a scalar quantity represented on a FEM mesh."""


    def __init__(self, wguide: Any, n_points: int = 500) -> None:
        """Build structure to hold a scalar field represented on a FEM mesh.

        The FEM mesh is constructed from the mesh file mesh_mail_fname.
        n_points is the nominal linear resolution of the rectangular grid output,
        partitioned so that x-y pixels are roughly square.
        That is, in the final grid, nx * ny ~= n_points^2, and dx ~= dy.

        """

        self.fem_mesh = FemMesh()
        self.fem_mesh.build_from_gmsh_mail(wguide)

        self.quantity_name = ""
        self.file_suffix = ""
        self.n_points = n_points

        self.scalar_fields = []

        # get approx square pixel rectangular grid
        self.v_x, self.v_y = self.fem_mesh.rect_grid(n_points)

        # m_x, m_y = np.meshgrid(self.v_x, self.v_y)
        # v_x_flat = m_x.flatten('F')
        # v_y_flat = m_y.flatten('F')
        # self.interper = self.fem_mesh.make_interpolator_for_grid(v_x_flat, v_y_flat, len(self.v_x), len(self.v_y))

        self.x_lab = r"$x$ [μm]"
        self.y_lab = r"$y$ [μm]"
        self.d_lab = r"$d$ [μm]"

    def setup_scalar_properties(
        self, nm_eng: str, unit: str, nm_math: str, fname_suffix: str
    ) -> None:
        """Configure labels and output names for a scalar quantity.

        Args:
            nm_eng: Human-readable quantity name.
            unit: Unit label appended to plotted values.
            nm_math: Mathematical quantity label.
            fname_suffix: Suffix used for generated filenames.
        """
        self.dim = 1
        self.nm_eng = nm_eng
        self.nm_math = nm_math
        self.unit = unit
        self.fname_suffix = fname_suffix

    def setup_vector_properties(
        self,
        dim: int,
        nm_eng: str,
        unit: str,
        nm_math: str,
        nm_math_comps: Sequence[str],
        fname_suffix: str,
        fname_suffix_comps: Sequence[str],
    ) -> None:
        """Configure labels and output names for a vector quantity.

        Args:
            dim: Number of vector components.
            nm_eng: Human-readable quantity name.
            unit: Unit label appended to plotted values.
            nm_math: Mathematical vector label.
            nm_math_comps: Mathematical labels for each vector component.
            fname_suffix: Base suffix used for generated filenames.
            fname_suffix_comps: Suffixes used for each component filename.
        """
        self.dim = dim
        self.nm_eng = nm_eng
        self.nm_math = nm_math
        self.unit = unit
        self.nm_math_comps = nm_math_comps
        self.fname_suffix = fname_suffix
        self.fname_suffix_comps = fname_suffix_comps

    def set_quantity_name(self, nm: str, suf: str) -> None:
        """Set the human-readable quantity name and output filename suffix."""
        self.quantity_name = nm
        self.file_suffix = suf

    def n_elts(self) -> int:  # NAME PROBLEMATIC
        """Return the number of finite elements in the underlying mesh."""
        return self.fem_mesh.n_msh_el

    def fill_quantity_by_material_index(self, mati_to_scalars: np.ndarray) -> None:
        """Populate scalar fields from values indexed by material.

        Args:
            mati_to_scalars: Scalar or vector values indexed by zero-based
                material index, with vacuum at index zero.
        """

        # make both scalar and vector inputs usable
        if self.dim > 1:
            mati_to_vals = mati_to_scalars
        else:
            mati_to_vals = np.expand_dims(copy.deepcopy(mati_to_scalars), axis=1)

        n_elts = self.fem_mesh.n_msh_el
        self.scalar_fields = []
        for i in range(self.dim):
            t_scalar_field = np.zeros([6, n_elts], dtype=np.float64)

            for i_el in range(n_elts):
                matel = self.element_to_material_index(i_el)
                t_scalar_field[:, i_el] = np.real(mati_to_vals[matel, i])

            t_scalar_field = t_scalar_field.flatten("F")
            self.scalar_fields.append(t_scalar_field)

    def element_to_material_index(self, i_el: int) -> int:
        """Return the zero-based material index for one mesh element."""
        return self.fem_mesh.element_to_material_index(i_el)

    def make_plot_1D_cut(
        self,
        pref: Union[str, Path],
        s_cut: str,
        val1: Union[Real, Sequence[Real]],
        val2: Optional[Sequence[Real]] = None,
    ) -> list[str]:
        """Make a 1D plot of the scalar along a line.

        Args:
            pref: Output filename prefix relative to the NumBAT output directory.
            s_cut: Cut direction: ``"x"``, ``"y"``, or ``"line"``.
            val1: Fixed coordinate for an axis-aligned cut, or the starting
                ``(x, y)`` point for a line cut.
            val2: Ending ``(x, y)`` point for a line cut.

        Returns:
            Paths of the generated cut plot files.

        For ``s_cut="x"``, the function is plotted along y at fixed x=val1.
        For ``s_cut="y"``, the function is plotted along x at fixed y=val1.
        For ``s_cut="line"``, the function is plotted from ``val1`` to ``val2``.
        """

        # full out-directory parth
        prefix = str(Path(runtime.outdir(), pref))

        validate_1D_cut_format(s_cut, val1, val2)

        fn=''
        if s_cut.lower() == "x":
            fn=self._make_plot_xcut(prefix, val1)
        elif s_cut.lower() == "y":
            fn=self._make_plot_ycut(prefix, val1)
        else:
            fn=self._make_plot_1D(prefix, val1, val2)

        return [fn]

    def _make_plot_ycut(self, prefix, y0):
        v_x_flat = self.v_x
        v_y_flat = y0 + np.zeros(len(self.v_x))
        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(self.v_x), 1
        )

        ## TODO: explain the need for the transpose in the imshow call below
        # epslo, epshi = np_min_max(m_scalar)

        fig, ax = plt.subplots()
        ax.set_xlabel(self.x_lab)
        ylab = self.nm_eng + " " + self.nm_math + " " + self.unit
        ax.set_ylabel(ylab)

        for i in range(self.dim):
            v_scalar = self.interper(self.scalar_fields[i])
            ax.plot(self.v_x, v_scalar)

        fname = prefix + "-" + self.fname_suffix + "_ycut.png"
        plottools.save_and_close_figure( fig, fname )
        return fname

    def _make_plot_xcut(self, prefix, x0):
        v_y_flat = self.v_y
        v_x_flat = x0 + np.zeros(len(self.v_y))

        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(self.v_y), 1
        )

        fig, ax = plt.subplots()
        ax.set_xlabel(self.y_lab)

        ylab = self.nm_eng + " " + self.nm_math + " " + self.unit
        ax.set_ylabel(ylab)

        for i in range(self.dim):
            v_scalar = self.interper(self.scalar_fields[i])
            ax.plot(self.v_y, v_scalar)

        fname = prefix + "-" + self.fname_suffix + "_ycut.png"
        plottools.save_and_close_figure( fig, fname)
        return fname


    def _make_plot_1D(self, prefix, pt0, pt1):

        x0, y0 = pt0
        x1, y1 = pt1
        v_x_flat = np.linspace(x0, x1, self.n_points)
        v_y_flat = np.linspace(y0, y1, self.n_points)
        v_d = np.sqrt((v_x_flat - v_x_flat[0]) ** 2 + (v_y_flat - v_y_flat[0]) ** 2)

        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(v_x_flat), 1
        )

        fig, ax = plt.subplots()
        ax.set_xlabel(self.d_lab)
        ylab = self.nm_eng + " " + self.nm_math + " " + self.unit
        ax.set_ylabel(ylab)

        for i in range(self.dim):
            v_scalar = self.interper(self.scalar_fields[i])
            ax.plot(v_d, v_scalar)

        fname = prefix + "-" + self.fname_suffix + "_linecut.png"
        plottools.save_and_close_figure(fig, fname)
        return fname

    def make_plot_2D(
        self, pref: Union[str, Path], aspect: float = 1.0, with_cb: bool = True
    ) -> list[str]:
        """Render every scalar-field component as a two-dimensional image.

        Args:
            pref: Output filename prefix relative to the NumBAT output directory.
            aspect: Matplotlib aspect ratio for the image axes.
            with_cb: Whether to add a colour bar to every image.

        Returns:
            Paths of the generated component image files.
        """

        fnames = []
        # full out-directory path
        prefix = str(Path(runtime.outdir(), pref))

        m_x, m_y = np.meshgrid(self.v_x, self.v_y)
        v_x_flat = m_x.flatten("F")
        v_y_flat = m_y.flatten("F")
        self.interper = self.fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, len(self.v_x), len(self.v_y)
        )

        for i in range(self.dim):
            m_scalar = self.interper(self.scalar_fields[i])

            # TODO: explain the need for the transpose in the imshow call below
            epslo, epshi = np_min_max(m_scalar)

            fig, ax = plt.subplots()
            cmap = "cool"
            im = ax.imshow(
                m_scalar.T,
                cmap=cmap,
                vmin=1.0,
                vmax=epshi,
                origin="lower",
                extent=[
                    np.min(self.v_x),
                    np.max(self.v_x),
                    np.min(self.v_y),
                    np.max(self.v_y),
                ],
            )
            # cf=ax.contourf(m_regx, m_regy, v_regindex, cmap=cmap, vmin=1.0, vmax=np.nanmax(v_regindex))
            ax.set_xlabel(self.x_lab)
            ax.set_ylabel(self.y_lab)
            ax.set_aspect(aspect)
            im.set_clim(1, np.nanmax(m_scalar))
            if with_cb:
                ticks = np.linspace(epslo, epshi, 5)
                fmts = (
                    (f"{epslo:.4f}",),
                    map(lambda x: f"{x:.1f}", ticks[1:-1]),
                    (f"{epshi:.4f}",),
                )
                fmts = list(itertools.chain.from_iterable(fmts))

                fmt = mticker.FixedFormatter(fmts)
                cax = ax.inset_axes([1.04, 0.1, 0.03, 0.8])
                cb = fig.colorbar(im, cax=cax, ticks=ticks, format=fmt)

                if self.dim == 1:
                    fname = prefix + "-" + self.fname_suffix + ".png"
                    qlab = self.nm_eng + " " + self.nm_math + " " + self.unit
                else:
                    fname = (
                        prefix
                        + "-"
                        + self.fname_suffix
                        + f"_{self.fname_suffix_comps[i]}.png"
                    )
                    qlab = self.nm_eng + " " + self.nm_math_comps[i] + " " + self.unit

                fig.suptitle(qlab)
                cb.ax.tick_params(labelsize=12)
                cb.outline.set_linewidth(0.5)
                cb.outline.set_color("gray")

            plottools.save_and_close_figure(fig, fname)
            fnames.append(fname)
        return fnames
