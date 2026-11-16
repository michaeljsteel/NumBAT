"""Sampling and interpolation of solver mode fields."""

import numpy as np

from numbat.field_selection import storage_key
from numbat.nbtypes import (
    ComplexPart,
    FieldComponent,
    FieldSelection,
    FieldType,
    SI_um,
    SI_vacuum_impedance_Z0,
)


def _field_storage_key(field_type, component, complex_part=None):
    return storage_key(FieldSelection(field_type, component, complex_part))


class ModeInterpolator:
    """Cache interpolation grids and sample fields for a simulation result."""

    def __init__(self, simresult):
        self.sim_result = simresult

        self.setup_for_2D_n_pts = 0
        self.setup_for_1D = ()

        self.xy_raw = {}
        self.xy_out = {}
        self.v_x_axis = None

        self.interper_f_2D = None
        self.interper_f_1D = None

        self.zero_arrays()

    def zero_arrays(self):
        self.interper_f_2D = None
        self.interper_f_1D = None

        self.xy_raw = {}
        self.xy_out = {}

    def cleanup(self):
        if self.interper_f_2D is not None:
            del self.interper_f_2D
        if self.interper_f_1D is not None:
            del self.interper_f_1D

        self.zero_arrays()

    def interpolate_mode_i(self, md, field_type, dims=2):
        """Interpolate one mode onto the configured 1D or 2D sampling grid."""
        simres = self.sim_result

        v_Fx6p, v_Fy6p, v_Fz6p, v_Fa6p = simres.get_mode_fields_on_fem_mesh(
            md, field_type
        )

        interper_f = self.interper_f_2D if dims == 2 else self.interper_f_1D

        m_ReFx = interper_f(v_Fx6p.real)
        m_ReFy = interper_f(v_Fy6p.real)
        m_ImFz = interper_f(v_Fz6p.imag)
        m_AbsF = interper_f(v_Fa6p)

        m_ImFx = interper_f(v_Fx6p.imag)
        m_ImFy = interper_f(v_Fy6p.imag)
        m_ReFz = interper_f(v_Fz6p.real)

        d_fields = {
            _field_storage_key(
                field_type, FieldComponent.X, ComplexPart.PRIMARY
            ): m_ReFx,
            _field_storage_key(
                field_type, FieldComponent.X, ComplexPart.SECONDARY
            ): m_ImFx,
            _field_storage_key(
                field_type, FieldComponent.Y, ComplexPart.PRIMARY
            ): m_ReFy,
            _field_storage_key(
                field_type, FieldComponent.Y, ComplexPart.SECONDARY
            ): m_ImFy,
            _field_storage_key(
                field_type, FieldComponent.Z, ComplexPart.PRIMARY
            ): m_ImFz,
            _field_storage_key(
                field_type, FieldComponent.Z, ComplexPart.SECONDARY
            ): m_ReFz,
            _field_storage_key(field_type, FieldComponent.MAGNITUDE): m_AbsF,
        }

        if field_type == FieldType.EM_H:
            for field in d_fields.values():
                field *= SI_vacuum_impedance_Z0

        d_fields[_field_storage_key(field_type, FieldComponent.INTENSITY)] = (
            m_AbsF**2
        )

        return d_fields

    def _choose_plot_points(self, n_pts):
        self.setup_for_2D_n_pts = n_pts

        fem_mesh = self.sim_result.fem_mesh
        x_min, x_max, y_min, y_max = fem_mesh.get_xy_limits()

        area = abs((x_max - x_min) * (y_max - y_min))
        n_pts_x = int(n_pts * abs(x_max - x_min) / np.sqrt(area))
        n_pts_y = int(n_pts * abs(y_max - y_min) / np.sqrt(area))

        shiftx, shifty = self.sim_result.get_xyshift()

        v_x = np.linspace(x_min, x_max, n_pts_x)
        v_y = np.linspace(y_min, y_max, n_pts_y)
        m_x, m_y = np.meshgrid(v_x, v_y)

        self.xy_raw = {"v_x": v_x, "v_y": v_y, "m_x": m_x, "m_y": m_y}

        v_x_out = (v_x + shiftx) / SI_um
        v_y_out = (v_y + shifty) / SI_um
        m_x_out, m_y_out = np.meshgrid(v_x_out, v_y_out)

        self.xy_out = {
            "v_x": v_x_out,
            "v_y": v_y_out,
            "m_x": m_x_out,
            "m_y": m_y_out,
        }

        print(
            "Structure has raw domain (x,y) = "
            f" [{v_x[0] / SI_um:.5f}, {v_x[-1] / SI_um:.5f}] x"
            f" [{v_y[0] / SI_um:.5f}, {v_y[-1] / SI_um:.5f}] (μm),"
            "\n             mapped to (x',y') = "
            f" [{v_x_out[0]:.5f}, {v_x_out[-1]:.5f}] x"
            f" [{v_y_out[0]:.5f}, {v_y_out[-1]:.5f}] (μm)"
        )
        return shiftx, shifty

    def define_plot_grid_2D(self, n_pts=501):
        """Define an interpolation grid for a nominal ``n_pts**2`` samples."""
        if self.setup_for_2D_n_pts == n_pts:
            return

        fem_mesh = self.sim_result.fem_mesh
        shiftx, shifty = self._choose_plot_points(n_pts)

        v_x_flat = self.xy_raw["m_x"].flatten("F") - shiftx
        v_y_flat = self.xy_raw["m_y"].flatten("F") - shifty

        nx, ny = len(self.xy_out["v_x"]), len(self.xy_out["v_y"])
        self.interper_f_2D = fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, nx, ny
        )

    def define_plot_grid_1D(self, s_cut, val1, val2, n_pts):
        """Define an interpolation grid for an axis-aligned or arbitrary cut."""
        setup_key = str((s_cut, val1, val2, n_pts))
        if self.setup_for_1D == setup_key:
            return
        self.setup_for_1D = setup_key

        fem_mesh = self.sim_result.fem_mesh
        x_min, x_max, y_min, y_max = fem_mesh.get_xy_limits()

        match s_cut.lower():
            case "x":
                v_y_flat = np.linspace(y_min, y_max, n_pts)
                v_x_flat = val1 * SI_um + np.zeros(n_pts)
            case "y":
                v_x_flat = np.linspace(x_min, x_max, n_pts)
                v_y_flat = val1 * SI_um + np.zeros(n_pts)
            case "line":
                v_x_flat = np.linspace(val1[0] * SI_um, val2[0] * SI_um, n_pts)
                v_y_flat = np.linspace(val1[1] * SI_um, val2[1] * SI_um, n_pts)

        nx = len(v_x_flat)
        self.interper_f_1D = fem_mesh.make_interpolator_for_grid(
            v_x_flat, v_y_flat, nx, 1
        )

        shiftx, shifty = self.sim_result.get_xyshift()
        v_x_flat -= shiftx
        v_y_flat -= shifty

        v_z_1D = np.sqrt(
            (v_x_flat - v_x_flat[0]) ** 2 + (v_y_flat - v_y_flat[0]) ** 2
        )

        match s_cut.lower():
            case "x":
                v_x_axis = v_y_flat
            case "y":
                v_x_axis = v_x_flat
            case "line":
                v_x_axis = v_z_1D
        self.v_x_axis = v_x_axis / SI_um
