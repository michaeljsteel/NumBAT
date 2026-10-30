import matplotlib.patches as mplpatches

from numbat_core.mesh.geometry import UserGeometryBase

from ._common import nmtoum


def make_onion_subs():
    subs = [
        ("dx_in_nm = 2000;", "dx_in_nm = %f;", "domain_x"),
        ("dy_in_nm = 2000;", "dy_in_nm = %f;", "domain_y"),
        ("a1 = 100;", "a1 = %f;", "inc_a_x"),
        ("a2 = 100;", "a2 = %f;", "inc_b_x"),
        ("a3 = 100;", "a3 = %f;", "inc_c_x"),
        ("a4 = 100;", "a4 = %f;", "inc_d_x"),
        ("a5 = 100;", "a5 = %f;", "inc_e_x"),
        ("a6 = 100;", "a6 = %f;", "inc_f_x"),
        ("a7 = 100;", "a7 = %f;", "inc_g_x"),
        ("a8 = 100;", "a8 = %f;", "inc_h_x"),
        ("a9 = 100;", "a9 = %f;", "inc_i_x"),
        ("a10 = 100;", "a10 = %f;", "inc_j_x"),
        ("a11 = 100;", "a11 = %f;", "inc_k_x"),
        ("a12 = 100;", "a12 = %f;", "inc_l_x"),
        ("a13 = 100;", "a13 = %f;", "inc_m_x"),
        ("a14 = 100;", "a14 = %f;", "inc_n_x"),
        ("a15 = 100;", "a15 = %f;", "inc_o_x"),
        ("lc = 0.1;", "lc = %f;", "lc_bkg"),
        ("lc_refine_1 = lc/1;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
        ("lc_refine_2 = lc/1;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
    ]

    return subs


def draw_onion_frame(ax, umb, styles):
    EDGE_COLOR = styles['edgecolor']
    LINEWIDTH = styles['linewidth']


    layers = (
        "inc_a_x",
        "inc_b_x",
        "inc_c_x",
        "inc_d_x",
        "inc_e_x",
        "inc_f_x",
        "inc_g_x",
        "inc_h_x",
        "inc_i_x",
        "inc_j_x",
        "inc_k_x",
        "inc_l_x",
        "inc_m_x",
        "inc_n_x",
        "inc_o_x",
    )

    rad = 0
    for sl in layers:
        lwid = umb.get_param(sl)
        if lwid is not None:
            if sl == "inc_a_x":
                rad += lwid / 2  # inc_a_x is diameter
            else:
                rad += lwid
            ax.add_patch(
                mplpatches.Circle(
                    (0, 0),
                    rad * nmtoum,
                    facecolor=None,
                    fill=False,
                    edgecolor=EDGE_COLOR,
                    linewidth=LINEWIDTH,
                )
            )


class Onion(UserGeometryBase):
    """NumBAT geometry template for a many-layer circular waveguide in a square domain."""

    def init_geometry(self):
        self.set_properties("onion", is_curvi=True)
        self.set_required_parameters(["inc_a_x"], num_req_mats=2)
        self.set_allowed_parameters(
            [
                "lc_refine_1",
                "lc_refine_2",
                "inc_b_x",
                "inc_c_x",
                "inc_d_x",
                "inc_e_x",
                "inc_f_x",
                "inc_g_x",
                "inc_h_x",
                "inc_i_x",
                "inc_j_x",
                "inc_k_x",
                "inc_l_x",
                "inc_m_x",
                "inc_n_x",
                "inc_o_x",
            ],
            num_allowed_mats=16,
        )

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class Onion1(UserGeometryBase):
    """NumBAT geometry template for a one-layer circular waveguide in a square domain."""

    def init_geometry(self):

        self.set_properties("onion1", is_curvi=True)

        self.set_required_parameters(["inc_a_x"], num_req_mats=2)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_2"])
        self.set_parameter_help(
            {
                "inc_a_x": "diameter of central cylinder",
                "material_a": "material of central cylinder",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_2": "mesh refinement on cylinder 1",
            }
        )

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        diam_a = self.get_param("inc_a_x")

        msg = ""

        if diam_a >= dom_x:
            msg += "Waveguide cylinder a (inc_a_x) has diameter larger than domain width (domain_x).\n"
        if diam_a >= dom_y:
            msg += "Waveguide cylinder a (inc_a_x) has diameter larger than domain height (domain_y).\n"

        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class Onion2(UserGeometryBase):
    """NumBAT geometry template for a two-layer circular waveguide in a square domain."""

    def init_geometry(self):

        self.set_properties("onion2", is_curvi=True)

        self.set_required_parameters(["inc_a_x", "inc_b_x"], num_req_mats=3)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_2"])

        self.set_parameter_help(
            {
                "inc_a_x": "diameter of central (a) cylinder",
                "inc_b_x": "annular radius of second (b) ring",
                "material_a": "material of central (a) cylinder",
                "material_a": "material of second (b) ring",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_2": "mesh refinement on cylinders",
            }
        )

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_a = self.get_param("inc_a_x") / 2.0
        rad_ann_b = self.get_param("inc_b_x")

        msg = ""

        diam_outer = 2 * (rad_a + rad_ann_b)

        if diam_outer >= dom_x:
            msg += "Outer cylinder has total diameter larger than domain width (domain_x).\n"
        if diam_outer >= dom_y:
            msg += "Outer cylinder has total diameter larger than domain height (domain_y).\n"

        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class Onion3(UserGeometryBase):
    """NumBAT geometry template for a three-layer circular waveguide in a square domain."""

    def init_geometry(self):

        self.set_properties("onion3", is_curvi=True)

        self.set_required_parameters(["inc_a_x", "inc_b_x", "inc_c_x"], num_req_mats=4)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_2"])

        self.set_parameter_help(
            {
                "inc_a_x": "diameter of central (a) cylinder",
                "inc_b_x": "annular radius of second (b) ring",
                "inc_c_x": "annular radius of third (c) ring",
                "material_a": "material of central (a) cylinder",
                "material_b": "material of second (b) ring",
                "material_c": "material of third (c) ring",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_2": "mesh refinement on cylinders",
            }
        )

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_a = self.get_param("inc_a_x") / 2.0
        rad_ann_b = self.get_param("inc_b_x")
        rad_ann_c = self.get_param("inc_c_x")

        msg = ""

        diam_outer = 2 * (rad_a + rad_ann_b + rad_ann_c)

        if diam_outer >= dom_x:
            msg += "Outer cylinder has total diameter larger than domain width (domain_x).\n"
        if diam_outer >= dom_y:
            msg += "Outer cylinder has total diameter larger than domain height (domain_y).\n"
        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class CircOnion(UserGeometryBase):
    """NumBAT geometry template for a many-layer circular waveguide in a circular domain."""

    def init_geometry(self):
        self.set_properties("circ_onion", is_curvi=True)
        self.set_required_parameters([
            "inc_a_x", "inc_b_x", "inc_c_x", "inc_d_x", "inc_e_x",
            "inc_f_x", "inc_g_x", "inc_h_x", "inc_i_x", "inc_j_x",
            "inc_k_x", "inc_l_x", "inc_m_x", "inc_n_x", "inc_o_x"], num_req_mats=16)


        self.set_allowed_parameters(["lc_bkg", "lc_refine_2"])

        self.set_parameter_help(
            {
                "inc_a_x": "diameter of central (a) cylinder",
                "inc_b_x": "annular radius of second (b) ring",
                "inc_c_x": "annular radius of third (c) ring",
                "inc_d_x": "annular radius of fourth (d) ring",
                "inc_e_x": "annular radius of fifth (e) ring",
                "inc_f_x": "annular radius of sixth (f) ring",
                "inc_g_x": "annular radius of seventh (g) ring",
                "inc_h_x": "annular radius of eigth (h) ring",
                "inc_i_x": "annular radius of ninth (i) ring",
                "inc_j_x": "annular radius of tenth (j) ring",
                "inc_k_x": "annular radius of eleventh (k) ring",
                "inc_l_x": "annular radius of twelth (l) ring",
                "inc_m_x": "annular radius of thirteenth (m) ring",
                "inc_n_x": "annular radius of fourteenth (n) ring",
                "inc_o_x": "annular radius of fifteenth (o) ring",

                "material_a_x": "material central (a) cylinder",
                "material_b_x": "material of second (b) ring",
                "material_c_x": "material of third (c) ring",
                "material_d_x": "material of fourth (d) ring",
                "material_e_x": "material of fifth (e) ring",
                "material_f_x": "material of sixth (f) ring",
                "material_g_x": "material of seventh (g) ring",
                "material_h_x": "material of eigth (h) ring",
                "material_i_x": "material of ninth (i) ring",
                "material_j_x": "material of tenth (j) ring",
                "material_k_x": "material of eleventh (k) ring",
                "material_l_x": "material of twelth (l) ring",
                "material_m_x": "material of thirteenth (m) ring",
                "material_n_x": "material of fourteenth (n) ring",
                "material_o_x": "material of fifteenth (o) ring",

                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_2": "mesh refinement on rings",
            }
        )

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_a = self.get_param("inc_a_x") / 2.0

        diam_outer = 2* rad_a
        for c in ('b','c','d','e','f','g','h','i','j','k','l','m','n','o'):
            t_rad = self.get_param('inc_%s_x' % c)
            diam_outer += t_rad * 2

        msg = ""

        if diam_outer >= dom_x:
            msg += f"Outer cylinder has total diameter of {diam_outer} which is larger than domain width (domain_x={dom_x}).\n"
        if diam_outer >= dom_y:
            msg += f"Outer cylinder has total diameter of {diam_outer} which is larger than domain height (domain_y={dom_y}).\n"

        dims_ok = not len(msg)
        return dims_ok, msg


    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class CircOnion1(UserGeometryBase):
    """NumBAT geometry template for a one-layer circular waveguide in a circular domain."""

    def init_geometry(self):
        self.set_properties("circ_onion1", is_curvi=True)

        self.set_required_parameters(["inc_a_x"], num_req_mats=2)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_2"])

        self.set_parameter_help(
            {
                "inc_a_x": "diameter of central (a) cylinder",
                "material_a": "material of central (a) cylinder",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_2": "mesh refinement on cylinders",
            }
        )

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_a = self.get_param("inc_a_x") / 2.0

        msg = ""

        diam_outer = 2 * (rad_a )

        if diam_outer >= dom_x:
            msg += "Outer cylinder has total diameter larger than domain width (domain_x).\n"
        if diam_outer >= dom_y:
            msg += "Outer cylinder has total diameter larger than domain height (domain_y).\n"
        dims_ok = not len(msg)
        return dims_ok, msg


class CircOnion2(UserGeometryBase):
    """NumBAT geometry template for a two-layer circular waveguide in a circular domain."""

    def init_geometry(self):
        self.set_properties("circ_onion2", is_curvi=True)

        self.set_required_parameters(["inc_a_x", "inc_b_x"], num_req_mats=3)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_2"])

        self.set_parameter_help(
            {
                 "inc_a_x": "diameter of central (a) cylinder",
                "inc_b_x": "annular radius of second (b) ring",
                "material_a": "material of central (a) cylinder",
                "material_b": "material of second (b) ring",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_2": "mesh refinement on cylinders",
            }
        )
    def apply_parameters(self):
        subs = make_onion_subs()
        return subs


    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_a = self.get_param("inc_a_x") / 2.0
        rad_ann_b = self.get_param("inc_b_x")

        msg = ""

        diam_outer = 2 * (rad_a + rad_ann_b )

        if diam_outer >= dom_x:
            msg += "Outer cylinder has total diameter larger than domain width (domain_x).\n"
        if diam_outer >= dom_y:
            msg += "Outer cylinder has total diameter larger than domain height (domain_y).\n"
        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class CircOnion3(UserGeometryBase):
    """NumBAT geometry template for a three-layer circular waveguide in a circular domain."""

    def init_geometry(self):

        self.set_properties("circ_onion3", is_curvi=True)
        self.set_required_parameters(["inc_a_x", "inc_b_x", "inc_c_x"], num_req_mats=4)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_2"])

        self.set_parameter_help(
            {
                "inc_a_x": "diameter of central (a) cylinder",
                "inc_b_x": "annular radius of second (b) ring",
                "inc_c_x": "annular radius of third (c) ring",
                "material_a": "material of central (a) cylinder",
                "material_b": "material of second (b) ring",
                "material_c": "material of third (c) ring",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_2": "mesh refinement on cylinders",
            }
        )

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_a = self.get_param("inc_a_x") / 2.0
        rad_ann_b = self.get_param("inc_b_x")
        rad_ann_c = self.get_param("inc_c_x")

        msg = ""

        diam_outer = 2 * (rad_a + rad_ann_b + rad_ann_c)

        if diam_outer >= dom_x:
            msg += "Outer cylinder has total diameter larger than domain width (domain_x).\n"
        if diam_outer >= dom_y:
            msg += "Outer cylinder has total diameter larger than domain height (domain_y).\n"
        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)
