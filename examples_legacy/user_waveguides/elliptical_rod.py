"""A source-tree NumBAT waveguide module containing an elliptical rod."""

from pathlib import Path

import matplotlib.patches as mplpatches

from usermesh import UserGeometryBase


class EllipticalRod(UserGeometryBase):
    """A single elliptical core suspended in a background material."""

    def init_geometry(self):
        self.set_properties("elliptical_rod", is_curvi=True)
        self.set_required_parameters(["inc_a_x", "inc_a_y"], 2)
        self.set_allowed_parameters(["lc_refine_1", "lc_refine_2"])
        self.set_parameter_help(
            {
                "inc_a_x": "horizontal diameter of the elliptical core",
                "inc_a_y": "vertical diameter of the elliptical core",
                "lc_bkg": "background mesh size",
                "lc_refine_1": "mesh refinement factor at the core boundary",
                "lc_refine_2": "mesh refinement factor at the core centre",
            }
        )

    def gmsh_template_path(self, p_dir_templates):
        return Path(__file__).with_name("elliptical_rod_msh_template.geo")

    def apply_parameters(self):
        return [
            ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
            ("a1 = 20;", "a1 = %f;", "inc_a_x"),
            ("a1y = 10;", "a1y = %f;", "inc_a_y"),
            ("lc = 0.1;", "lc = %f;", "lc_bkg"),
            ("lc_refine_1 = lc/1;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc/1;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
        ]

    def draw_mpl_frame(self, ax, styles):
        ax.add_patch(
            mplpatches.Ellipse(
                (0, 0),
                self.get_param("inc_a_x") / 1000,
                self.get_param("inc_a_y") / 1000,
                facecolor="none",
                edgecolor=styles["edgecolor"],
                linewidth=styles["linewidth"],
            )
        )
