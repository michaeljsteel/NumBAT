import matplotlib.patches as mplpatches
import numpy as np

from usermesh import UserGeometryBase

from ._common import _process_one_and_two_incls, _process_one_and_two_incls_subs, nmtoum


class Circular(UserGeometryBase):
    """NumBAT geometry template for a circular waveguide."""

    def init_geometry(self):
        self.set_properties("circular")
        self.set_allowed_parameters(
            [
                "inc_a_y",
                "inc_b_x",
                "inc_b_y",
                "lc_refine_1",
                "lc_refine_2",
            ],
            3,
        )

        gmshfile, nelts = _process_one_and_two_incls(self._d_params)
        self._gmsh_template_filename = gmshfile  # special case where Circular and Rectangular share common gmshfile, so geom name and geom file are different

    def apply_parameters(self):

        subs = _process_one_and_two_incls_subs(self._gmsh_template_filename)
        subs.append(("rect = 1;", "rect = %d;", "0"))  # apply circularness

        return subs

    def draw_mpl_frame(self, ax, styles):

        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']

        rad = self.get_param("inc_a_x") * 0.5

        circ = mplpatches.Circle(
            (0, 0),
            rad * nmtoum,
            facecolor=None,
            fill=False,
            edgecolor=EDGE_COLOR,
            linewidth=LINEWIDTH,
        )
        ax.add_patch(circ)


class Rectangular(UserGeometryBase):
    """NumBAT geometry template for a rectangular waveguide."""

    def init_geometry(self):
        self.set_properties("rectangular")

        self.set_required_parameters(["inc_a_x"], 2)

        self.set_allowed_parameters(
            [
                "inc_a_y",
                "inc_b_x",
                "inc_b_y",
                "lc_refine_1",
                "lc_refine_2",
            ],
            3,
        )

        gmshfile, nelts = _process_one_and_two_incls(self._d_params)

        self._gmsh_template_filename = gmshfile  # special case where Circular and Rectangular share common gmshfile, so geom name and geom file are different

    def apply_parameters(self):
        subs = _process_one_and_two_incls_subs(self._gmsh_template_filename)

        return subs

    def draw_mpl_frame(self, ax, styles):

        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']

        wid = self.get_param("inc_a_x") * nmtoum
        hgt = self.get_param("inc_a_y") * nmtoum

        ax.add_patch(
            mplpatches.Rectangle(
                (-wid / 2, -hgt / 2),
                wid,
                hgt,
                facecolor=None,
                fill=False,
                edgecolor=EDGE_COLOR,
                linewidth=LINEWIDTH
            )
        )

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        wid = self.get_param("inc_a_x")
        hgt = self.get_param("inc_a_y")

        msg = ""

        if wid >= dom_x:
            msg += "Waveguide width (inc_a_x) is larger than domain width (domain_x).\n"
        if hgt >= dom_y:
            msg += (
                "Waveguide height (inc_a_y) is larger than domain height (domain_y).\n"
            )

        dims_ok = not len(msg)
        return dims_ok, msg


class Triangular(UserGeometryBase):
    """NumBAT geometry template for a triangular waveguide."""

    def init_geometry(self):
        self.set_properties("triangular")
        self.set_required_parameters(
            ["base_width", "peak_xoff", "peak_height"], num_req_mats=1
        )
        self.set_allowed_parameters(
            ["lc_bkg", "lc_refine_1", "lc_refine_2"], num_allowed_mats=2
        )
        self.set_parameter_help(
            {
                "base_width": "length of base of triangle along x-axis",
                "peak_xoff": "horizontal offset of peak along x-axis from left vertex",
                "peak_height": "perpendicular height of triangle",
                "material_bkg": "background material",
                "material_a": "material of triangular core",
            }
        )

    def apply_parameters(self):
        subs = [
            ("dx_in_nm = 1000.0;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 1000.0;", "dy_in_nm = %f;", "domain_y"),
            ("base_width = 600.0;", "base_width = %f;", "base_width"),
            ("peak_xoff = 200.0;", "peak_xoff = %f;", "peak_xoff"),
            ("peak_height = 400.0;", "peak_height = %f;", "peak_height"),
            ("lc = 0.1;", "lc = %f;", "lc_bkg"),
            ("lc_refine_1 = lc/1.0;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc/1.0;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
        ]

        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        base_wid = self.get_param("base_width")
        peak_xoff = self.get_param("peak_xoff")
        peak_height = self.get_param("peak_height")

        peak_locx = -base_wid / 2 + peak_xoff

        msg = ""

        if base_wid >= dom_x:
            msg += "Waveguide base width (base_width) is larger than the domain width (domain_x).\n"

        if peak_locx < -dom_x / 2 or peak_locx > dom_x / 2:
            msg += "Waveguide peak is outside the x-domain (domain_x).\n"
        if peak_height > dom_y / 2 or peak_height < -dom_y / 2:
            msg += "Waveguide height (peak_height) is too large for the domain height (domain_y).\n"

        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']

        wid = self.get_param("base_width") * nmtoum
        xoff = self.get_param("peak_xoff") * nmtoum
        hgt = self.get_param("peak_height") * nmtoum
        vertices = np.array(
            [[-wid / 2, -hgt / 2], [-wid / 2 + xoff, hgt / 2], [wid / 2, -hgt / 2]]
        )

        ax.add_patch(
            mplpatches.Polygon(
                vertices, facecolor=None, fill=False, edgecolor=EDGE_COLOR, linewidth=LINEWIDTH
            )
        )
