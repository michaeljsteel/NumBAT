# User NumBAT mesh implementation file

import matplotlib.patches as mplpatches
import numpy as np

from usermesh import UserGeometryBase
import reporting

nmtoum = 0.001  # template radii are in nm but matplotlib plots are in microns


#LINEWIDTH=0.75

#TODO:
# expand allowed_/reqd_ parameters to take a 'is_dim' field.
# have objects.Structure determine whether units are microns or nm.
# multiply by nm_to_um  in get_param if required.



def _process_one_and_two_incls(params):
    nelts = 0
    gmshfile = ""

    if "slab_b_x" in params:
        raise ValueError(
            f"NumBAT doesn't understand your geometry: with shape {params['inc_shape']}, I did not expect values for slab_b."
        )

    if "slab_a_x" in params:
        raise ValueError(
            f"NumBAT doesn't understand your geometry: with shape {params['inc_shape']}, I did not expect values for slab_a."
        )

    if "inc_a_x" in params:
        if (
            "coat_y" not in params and "inc_b_x" not in params
        ):  # One inclusion, no coating
            gmshfile = "oneincl"  # used to be just '1'
            nelts = 2  # bkg, core (mat_a)

        elif (
            "coat_y" not in params and "inc_b_x" in params
        ):  # Two inclusions, no coating
            gmshfile = "twoincl"  # used to be just '2'
            nelts = 3  # bkg, core 1 (mat_a), core 2 (mat_b)

        # Two inclusions, with coating # TODO:implement
        elif "coat_y" in params and "inc_b_x" in params:
            raise NotImplementedError("Have not implemented 2 coated inclusions.")

        elif (
            "coat_y" in params and "inc_b_x" not in params
        ):  # One inclusion, with coating # TODO:implement
            raise NotImplementedError("Have not implemented 1 coated inclusions.")

        else:
            raise ValueError("NumBAT doesn't understand your geometry.")
    else:
        raise ValueError("must have at least one nonzero inclusion.")

    return gmshfile, nelts


def _process_one_and_two_incls_subs(msh_template):
    # TODO: these are crazy small defaults
    subs = [
        ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
        ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
        ("a1 = 20;", "a1 = %f;", "inc_a_x"),
        ("a1y = 10;", "a1y = %f;", "inc_a_y"),
        ("lc = 0.1;", "lc = %f;", "lc_bkg"),
        ("lc_refine_1 = lc/1;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
        ("lc_refine_2 = lc/1;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
    ]

    if msh_template in ["twoincl", "2", "2_on_s", "2_on_2s"]:
        subs.append(("a2 = 10;", "a2 = %f;", "inc_b_x"))
        subs.append(("a2y = 20;", "a2y = %f;", "inc_b_y"))
        subs.append(("sep = 10;", "sep = %f;", "two_inc_sep"))

        # geo = geo.replace('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;' % umb.lc_refine_3)
    if msh_template == "2":
        subs.append(("yoff = -5;", "yoff = %f;", "incs_y_offset"))

    if msh_template in ["1_on_slab", "1_on_2slabs", "1_on_slab", "2_on_2slabs"]:
        subs.append(("slab_width = dx_in_nm;", "slab_width = %f;", "slab_a_x"))
        subs.append(("slab_height = 10;", "slab_height = %f;", ".slab_a_y"))
        subs.append(("lc_refine_3 = lc/1;", "lc_refine_3 = lc/%f;", "lc_refine_3"))
        subs.append(("lc_refine_4 = lc/1;", "lc_refine_4 = lc/%f;", "lc_refine_4"))

    if msh_template in ["1_on_2slabs", "2_on_2slabs"]:
        subs.append(("slab2_width = dx_in_nm;", "slab2_width = %f;", "slab_b_x"))
        subs.append(("slab2_height = 5;", "slab2_height = %f;", "slab_b_y"))
        subs.append(("lc_refine_3 = lc/1;", "lc_refine_3 = lc/%f;", ".lc_refine_3"))
        subs.append(("lc_refine_4 = lc/1;", "lc_refine_4 = lc/%f;", "lc_refine_4"))

    return subs


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


class TwoIncl(UserGeometryBase):
    """NumBAT geometry template for a double inclusion waveguide."""

    def init_geometry(self):
        gmshfile, nelts = _process_one_and_two_incls(self._d_params)
        self.set_properties("twoincl", is_curvi=True)
        self._gmsh_template_filename = gmshfile  # special case where Circular and Rectangular share common gmshfile, so geom name and geom file are different

    def apply_parameters(self):

        subs = _process_one_and_two_incls_subs(self._gmsh_template_filename, self)

        return subs

    def draw_mpl_frame(self, ax, styles):
        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']


        widl = self.get_param("inc_a_x") * nmtoum
        hgtl = self.get_param("inc_a_y") * nmtoum
        widr = self.get_param("inc_b_x") * nmtoum
        hgtr = self.get_param("inc_b_y") * nmtoum
        sep = self.get_param("two_inc_sep") * nmtoum
        yoff = self.get_param("yoff") * nmtoum

        shape = self.get_param("inc_shape")

        if shape == "circular":
            ax.add_patch(
                mplpatches.Circle(
                    (-sep / 2, 0),
                    widl,
                    facecolor=None,
                    fill=False,
                    edgecolor=EDGE_COLOR,
                    linewidth=LINEWIDTH,
                )
            )

            ax.add_patch(
                mplpatches.Circle(
                    (sep / 2, yoff),
                    widr,
                    facecolor=None,
                    fill=False,
                    edgecolor=EDGE_COLOR,
                    linewidth=LINEWIDTH,
                )
            )

        else:
            ax.add_patch(
                mplpatches.Rectangle(
                    (-sep / 2 - widl / 2, -hgtl / 2),
                    widl,
                    hgtl,
                    facecolor=None,
                    fill=False,
                    edgecolor=EDGE_COLOR,
                    linewidth=LINEWIDTH,
                )
            )

            ax.add_patch(
                mplpatches.Rectangle(
                    (sep / 2 - widr / 2, yoff - hgtr / 2),
                    widr,
                    hgtr,
                    facecolor=None,
                    fill=False,
                    edgecolor=EDGE_COLOR,
                    linewidth=LINEWIDTH,
                )
            )


class TwoInclVert(UserGeometryBase):
    """NumBAT geometry template for a rectangular double inclusion waveguide arranged vertically."""

    def init_geometry(self):
        # gmshfile, nelts = _process_one_and_two_incls(self._d_params)
        self.set_properties("twoinclvert", is_curvi=True)

        self.set_required_parameters(
            ["inc_a_w", "inc_a_h", "inc_b_w", "inc_b_h", "inc_sep_x", "inc_sep_y"],
            num_req_mats=3,
        )
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 100.0;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50.0;", "dy_in_nm = %f;", "domain_y"),
            ("inc_a_w = 20.0;", "inc_a_w = %f;", "inc_a_w"),
            ("inc_a_h = 10.0;", "inc_a_h = %f;", "inc_a_h"),
            ("inc_b_w = 20.0;", "inc_b_w = %f;", "inc_b_w"),
            ("inc_b_h = 10.0;", "inc_b_h = %f;", "inc_b_h"),
            ("inc_sep_x = 5.0;", "inc_sep_x = %f;", "inc_sep_x"),
            ("inc_sep_y = 15.0;", "inc_sep_y = %f;", "inc_sep_y"),
            ("lc_bkg = 0.05;", "lc_bkg = %f;", "lc_bkg"),
            ("lc_refine_1 = lc/2.0;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc/2.0;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
        ]

        return subs

    def check_dimensions(self):
        """The first box must be higher. Which means yoff must be positive."""

        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        inc_a_w = self.get_param("inc_a_w")
        inc_a_h = self.get_param("inc_a_h")
        inc_b_w = self.get_param("inc_b_w")
        inc_b_h = self.get_param("inc_b_h")
        inc_sep_x = self.get_param("inc_sep_x")
        inc_sep_y = self.get_param("inc_sep_y")

        # waveguides can't touch so either they are fully separated in y, or in x
        # either way the upper one must be higher

        msg = ""

        if (inc_a_w + inc_sep_x) / 2 > dom_x / 2 or (
            -inc_a_w + inc_sep_x
        ) / 2 < -dom_x / 2:
            msg += "Waveguide 1 is falling outside the x-domain (Check parameters: inc_a_w, inc_sep_x, domain_x). \n"
        if (inc_a_h + inc_sep_y) / 2 > dom_y / 2 or (
            -inc_a_h + inc_sep_y
        ) / 2 < -dom_y / 2:
            msg += "Waveguide 1 is falling outside the x-domain (Check parameters: inc_a_h, inc_sep_y, domain_y). \n"

        if (inc_b_w + inc_sep_x) / 2 > dom_x / 2 or (
            -inc_b_w + inc_sep_x
        ) / 2 < -dom_x / 2:
            msg += "Waveguide 1 is falling outside the x-domain (Check parameters: inc_b_w, inc_sep_x, domain_x). \n"
        if (inc_b_h + inc_sep_y) / 2 > dom_y / 2 or (
            -inc_b_h + inc_sep_y
        ) / 2 < -dom_y / 2:
            msg += "Waveguide 1 is falling outside the x-domain (Check parameters: inc_b_h, inc_sep_y, domain_y). \n"

        yoverlap = inc_sep_y - (inc_a_h + inc_b_h) / 2
        minysep = inc_sep_y - max(inc_a_h, inc_b_h) / 2

        xoverlap = abs(inc_sep_x) - (inc_a_w + inc_b_w) / 2
        if inc_sep_y <= 0 or minysep <= 0:
            msg += "Vertical separation of the two guides must be positive (Check parameter: inc_sep_y) \n"

        if yoverlap <= 0 and xoverlap <= 0:
            msg += "The two waveguides are overlapping (Check parameters: inc_a_w, inc_a_h, inc_b_w, inc_b_h, inc_sep_x, inc_sep_y). \n"

        dims_ok = not len(msg)

        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']


        widu = self.get_param("inc_a_w") * nmtoum
        hgtu = self.get_param("inc_a_h") * nmtoum
        widl = self.get_param("inc_b_w") * nmtoum
        hgtl = self.get_param("inc_b_h") * nmtoum
        xsep = self.get_param("inc_sep_x") * nmtoum
        ysep = self.get_param("inc_sep_y") * nmtoum

        xu = -xsep / 2 - widu / 2
        yu = ysep / 2 - hgtu / 2
        xl = xsep / 2 - widl / 2
        yl = -ysep / 2 - hgtu / 2

        ax.add_patch(
            mplpatches.Rectangle(
                (xu, yu),
                widu,
                hgtu,
                facecolor=None,
                fill=False,
                edgecolor=EDGE_COLOR,
                linewidth=LINEWIDTH,
            )
        )

        ax.add_patch(
            mplpatches.Rectangle(
                (xl, yl),
                widl,
                hgtl,
                facecolor=None,
                fill=False,
                edgecolor=EDGE_COLOR,
                linewidth=LINEWIDTH,
            )
        )


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

        self.set_properties("onion3", desc, is_curvi=True)

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

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class CircOnion1(UserGeometryBase):
    """NumBAT geometry template for a one-layer circular waveguide in a circular domain."""

    def init_geometry(self):
        self.set_properties("circ_onion1", is_curvi=True)

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

    def draw_mpl_frame(self, ax, styles):
        draw_onion_frame(ax, self, styles)


class CircOnion2(UserGeometryBase):
    """NumBAT geometry template for a two-layer circular waveguide in a circular domain."""

    def init_geometry(self):
        self.set_properties("circ_onion2", 3, True)

    def apply_parameters(self):
        subs = make_onion_subs()
        return subs

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


class Pedestal(UserGeometryBase):
    """NumBAT geometry template for a pedestal-type waveguide."""

    def init_geometry(self):
        self.set_properties("pedestal")

        self.set_required_parameters(
            [
                "inc_a_x",
                "inc_b_x",
                "inc_c_x",
                "slab_a_x",
                "slab_b_x",
                "pillar_x",
                "pillar_y",
            ],
            num_req_mats=4,
        )
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])

        self.set_parameter_help(
            {
                "inc_a_x": "width of the top of the rib",
                "inc_a_y": "height of the top of the rib",
                "slab_a_x": "width of the middle of the rib",
                "slab_a_y": "height of the middle of the rib",
                "material_bkg": "material of background",
                "material_a": "material of rib",
                "material_b": "material of slab",
                "material_c": "material of pillar",
                "lc": "grid points arounds boundary as fraction of domain_x",
                "lc_refine1": "refined density along upper rib",
                "lc_refine2": "refined density along buried rib",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
            ("a1 = 20;", "a1 = %f;", "inc_a_x"),
            ("a1y = 10;", "a1y = %f;", "inc_a_y"),
            ("a1top = 15;", "a1top = %f;", "inc_b_x"),
            ("slabx = 80;", "slabx = %f;", "slab_a_x"),
            ("slaby = 10;", "slaby = %f;", "slab_a_y"),
            ("slabxtop = 60;", "slabxtop = %f;", "slab_b_x"),
            ("px = 2;", "px = %f;", "pillar_x"),
            ("py = 5;", "py = %f;", "pillar_y"),
            ("lc = 0.1;", "lc = %f;", "lc_bkg"),
            ("lc_refine_1 = lc/1;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc/1;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
        ]
        return subs


class TrapezoidalRib(UserGeometryBase):
    """NumBAT geometry template for a trapezoidal_rib waveguide."""

    def init_geometry(self):

        self.set_properties("trapezoidal_rib")

        self.set_required_parameters(
            [
                "rib_top_width",
                "rib_height",
                "rib_base_width",
                "slab_width",
                "slab_thickness",
            ],
            num_req_mats=4,
        )
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])

        self.set_parameter_help(
            {
                "rib_top_width": "width of the top of the rib",
                "rib_height": "height of the top of the rib",
                "rib_base_width": "width at base of the rib",
                "slab_width": "width of the slab",
                "slab_thickness": "thickness of the slab",
                "material_bkg": "material of background",
                "material_a": "material of elevated rib",
                "material_b": "material of buried rib",
                "material_c": "material of substrate",
                "lc": "grid points arounds boundary as fraction of domain_x",
                "lc_refine1": "refined density along upper rib",
                "lc_refine2": "refined density along buried rib",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 4000.0;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 2000.0;", "dy_in_nm = %f;", "domain_y"),
            ("rib_top_width = 600.0;", "rib_top_width = %f;", "rib_top_width"),
            ("rib_base_width = 900.0;", "rib_base_width = %f;", "rib_base_width"),
            ("rib_height = 500.0;", "rib_height = %f;", "rib_height"),
            ("slab_width = 1800.0;", "slab_width = %f;", "slab_width"),
            ("slab_thickness = 300.0;", "slab_thickness = %f;", "slab_thickness"),
            ("lc = 0.020000;", "lc = %f;", "lc_bkg"),
            ("lc_refine_1 = lc/10.0;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc/5.0;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
        ]

        return subs

    def check_dimensions(self):

        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rib_wbot = self.get_param("rib_base_width")

        # TODO: more  checks

        msg = ""

        if rib_wbot >= dom_x:
            msg = "Slab width (slab_b_x) is wider than domain width (domain_x).\n"

        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']

        rib_top_w = self.get_param("rib_top_width") * nmtoum
        rib_base_w = self.get_param("rib_base_width") * nmtoum
        rib_h = self.get_param("rib_height") * nmtoum
        slab_w = self.get_param("slab_width") * nmtoum
        slab_h = self.get_param("slab_thickness") * nmtoum


        vertices = np.array(
            [
                [-slab_w / 2, -slab_h],
                [-slab_w / 2, 0],
                [-rib_base_w / 2, 0],
                [-rib_top_w / 2, rib_h],
                [rib_top_w / 2, rib_h],
                [rib_base_w / 2, 0],
                [slab_w / 2, 0],
                [slab_w / 2, -slab_h],
                [-slab_w / 2, -slab_h],
            ]
        )

        ax.add_patch(
            mplpatches.Polygon(
                vertices, facecolor=None, fill=False, edgecolor=EDGE_COLOR, linewidth=LINEWIDTH
            )
        )

class Rib(UserGeometryBase):
    """NumBAT geometry template for a rib waveguide."""

    def init_geometry(self):
        if self._d_materials[
            "c"
        ].is_vacuum():  # TODO: perhaps a better test is whether bkg = mat_c
            nt = 3
        else:
            nt = 4

        nt = 3  # including bkg

        self.set_properties("rib")

        self.set_required_parameters(
            ["rib_w", "rib_h", "slab_w", "slab_h"], num_req_mats=nt
        )

        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])

        self.set_parameter_help(
            {
                "rib_w": "width of raised rib region",
                "rib_h": "height of raised rib region",
                "slab_w": "width of slab substrate region",
                "slab_h": "height of slab substrate region",
                "material_bkg": "background material",
                "material_a": "material of raised rib region",
                "material_b": "material of slab substrate region",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
            ("a1 = 20;", "a1 = %f;", "rib_w"),
            ("a1y = 10;", "a1y = %f;", "rib_h"),
            ("slabx = 80;", "slabx = %f;", "slab_w"),
            ("slaby = 10;", "slaby = %f;", "slab_h"),
            ("lc = 0.1;", "lc = %f;", "lc_bkg"),
            ("lc_refine_1 = lc/1;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc/1;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
        ]

        return subs

    def draw_mpl_frame(self, ax, styles):
        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']

        rib_w = self.get_param("rib_w") * nmtoum
        rib_h = self.get_param("rib_h") * nmtoum
        slab_w = self.get_param("slab_w") * nmtoum
        slab_h = self.get_param("slab_h") * nmtoum

        vertices = np.array(
            [
                [-slab_w / 2, -slab_h],
                [-slab_w / 2, 0],
                [-rib_w / 2, 0],
                [-rib_w / 2, rib_h],
                [rib_w / 2, rib_h],
                [rib_w / 2, 0],
                [slab_w / 2, 0],
                [slab_w / 2, -slab_h],
                [-slab_w / 2, -slab_h],
            ]
        )

        ax.add_patch(
            mplpatches.Polygon(
                vertices, facecolor=None, fill=False, edgecolor=EDGE_COLOR, linewidth=LINEWIDTH
            )
        )


class RibCoated(UserGeometryBase):
    """NumBAT geometry template for a coated rib waveguide."""

    def init_geometry(self):

        self.set_properties("rib_coated")

        self.set_required_parameters(
            ["rib_w", "rib_h", "slab_w", "slab_h", "coat_w", "coat_h"], num_req_mats=4
        )
        self.set_allowed_parameters(
            ["lc_bkg", "lc_refine_1", "lc_refine_2", "lc_refine_3"]
        )
        self.set_parameter_help(
            {
                "rib_w": "width of raised rib region",
                "rib_h": "height of raised rib region",
                "slab_w": "width of slab substrate region",
                "slab_h": "height of slab substrate region",
                "coat_w": "horizontal thickness of coating layer",
                "coat_h": "vertical thickness of coating layer",
                "material_bkg": "background material",
                "material_a": "material of raised rib core",
                "material_b": "material of substrate slab",
                "material_c": "material of coating layer",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
            ("a1 = 20;", "a1 = %f;", "rib_w"),
            ("a1y = 10;", "a1y = %f;", "rib_h"),
            ("slabx = 80;", "slabx = %f;", "slab_w"),
            ("slaby = 10;", "slaby = %f;", "slab_h"),
            ("coatx = 2;", "coatx = %f;", "coat_w"),
            ("coaty = 2;", "coaty = %f;", "coat_h"),
            ("lc = 0.1;", "lc = %f;", "lc_bkg"),
            ("lc_refine_1 = lc/1;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc/1;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
            ("lc_refine_3 = lc/1;", "lc_refine_3 = lc/%f;", "lc_refine_3"),
        ]

        return subs


class RibDoubleCoated(UserGeometryBase):
    """NumBAT geometry template for a double coated rib waveguide."""

    def init_geometry(self):
        self.set_properties("rib_double_coated")

        self.set_required_parameters(
            [
                "rib_w",
                "rib_h",
                "slab_w",
                "slab_h",
                "coat_w",
                "coat_h",
                "coat2_w",
                "coat2_h",
            ],
            num_req_mats=5,
        )
        self.set_allowed_parameters(
            [
                "lc_bkg",
                "lc_refine_1",
                "lc_refine_2",
                "lc_refine_3",
                "lc_refine_4",
                "lc_refine_5",
            ]
        )
        self.set_parameter_help(
            {
                "rib_w": "width of raised rib region",
                "rib_h": "height of raised rib region",
                "slab_w": "width of slab substrate region",
                "slab_h": "height of slab substrate region",
                "coat_w": "horizontal thickness of inner coating layer",
                "coat_h": "vertical thickness of inner coating layer",
                "coat2_w": "horizontal thickness of outer coating layer",
                "coat2_h": "vertical thickness of outer coating layer",
                "material_bkg": "background material",
                "material_a": "material of raised rib core",
                "material_b": "material of substrate slab",
                "material_c": "material of inner coating layer",
                "material_d": "material of outer coating layer",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
            ("a1 = 20;", "a1 = %f;", "rib_w"),
            ("a1y = 10;", "a1y = %f;", "rib_h"),
            ("slabx = 80;", "slabx = %f;", "slab_w"),
            ("slaby = 10;", "slaby = %f;", "slab_h"),
            ("coatx = 2;", "coatx = %f;", "coat_w"),
            ("coaty = 2;", "coaty = %f;", "coat_h"),
            ("slab2y = 5;", "slab2y = %f;", "slab_b_y"),
            ("coat2x = 4;", "coat2x = %f;", "coat2_w"),
            ("coat2y = 4;", "coat2y = %f;", "coat2_h"),
            ("lc_bkg = 0.05;", "lc = %f;", "lc_bkg"),
            ("lc_refine_1 = lc_bkg/2.0;", "lc_refine_1 = lc_bkg/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc_bkg/2.0;", "lc_refine_2 = lc_bkg/%f;", "lc_refine_2"),
            ("lc_refine_3 = lc_bkg/2.0;", "lc_refine_3 = lc_bkg/%f;", "lc_refine_3"),
            ("lc_refine_4 = lc_bkg/2.0;", "lc_refine_4 = lc_bkg/%f;", "lc_refine_4"),
            ("lc_refine_5 = lc_bkg/2.0;", "lc_refine_5 = lc_bkg/%f;", "lc_refine_5"),
        ]

        return subs


# TODO: more efficient way to do this logic?!
def _process_rib_mk_2(params):
    slab2_active = "slab2_h" in params and "material_c" in params
    slab3_active = "slab3_h" in params and "material_d" in params
    slab4_active = "slab4_h" in params and "material_e" in params
    slab5_active = "slab5_h" in params and "material_f" in params

    # warn if one is set but not the other
    slab2_incompl = ("slab2_h" in params or "material_c" in params) and not slab2_active
    slab3_incompl = ("slab3_h" in params or "material_d" in params) and not slab3_active
    slab4_incompl = ("slab4_h" in params or "material_e" in params) and not slab4_active
    slab5_incompl = ("slab5_h" in params or "material_f" in params) and not slab5_active

    if slab2_incompl:
        reporting.report("Incomplete parameters for Rib MkII slab 2. Disabling this slab.")
    if slab3_incompl:
        reporting.report("Incomplete parameters for Rib MkII slab 3. Disabling this slab.")
    if slab4_incompl:
        reporting.report("Incomplete parameters for Rib MkII slab 4. Disabling this slab.")
    if slab5_incompl:
        reporting.report("Incomplete parameters for Rib MkII slab 5. Disabling this slab.")

    if slab3_active and not slab2_active:
        reporting.warning("Can't have slab 3 active without active slab 2. Disabling slabs 3 and higher.")
        slab3_active=False

    if slab4_active and not slab3_active:
        reporting.warning("Can't have slab 4 active without active slabs 2 and 3. Disabling slabs 4 and higher.")
        slab4_active=False

    if slab5_active and not slab4_active:
        reporting.warning("Can't have slab 5 active without active slab 2, 3 and 4. Disabling slabs 5.")
        slab5_active=False

    if slab5_active:
        nslabs=5
        nelts=7
        #gmshfile = "ribmk2_5slab_msh_template.geo"
    elif slab4_active:
        nslabs=4
        nelts=6
        #gmshfile = "ribmk2_4slab_msh_template.geo"
    elif slab3_active:
        nslabs=3
        nelts=5
        #gmshfile = "ribmk2_3slab_msh_template.geo"
    elif slab2_active:
        nslabs=2
        nelts=4
        #gmshfile = "ribmk2_2slab_msh_template.geo"
    else:
        nslabs=1
        nelts=3
        #gmshfile = "ribmk2_1slab_msh_template.geo"


    return nslabs, nelts

class RibMkII(UserGeometryBase):
    """NumBAT geometry template for a rib mk 2 waveguide."""

    def init_geometry(self):
        self.set_properties("ribmk2")


        self.set_required_parameters( ["rib_w", "rib_h", "slab1_h"], num_req_mats=3)

        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", 
                                     "slab2_h", "slab3_h", "slab4_h", "slab5_h", "nslabs"], num_allowed_mats=7) # How specify mats_?

        nslabs, nelts = _process_rib_mk_2(self._d_params)

        self._d_params['nslabs'] = nslabs

        self.set_parameter_help(
            {
                "rib_w": "width of raised rib region",
                "rib_h": "height of raised rib region",
                "slab1_h": "height of top slab substrate region",
                "slab2_h": "height of second slab substrate region [optional]",
                "slab3_h": "height of third slab substrate region [optional]",
                "slab4_h": "height of fourth slab substrate region [optional]",
                "material_bkg": "background material",
                "material_a": "material of rib",
                "material_b": "material of top slab substrate region",
                "material_c": "material of second slab substrate region [optional]",
                "material_d": "material of third slab substrate region [optional]",
                "material_e": "material of fourth slab substrate region [optional]",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 1500.0;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 1500.0;", "dy_in_nm = %f;", "domain_y"),
            ("un_rib_w = 300.0;", "un_rib_w = %f;", "rib_w"),
            ("un_rib_h = 200.0;", "un_rib_h = %f;", "rib_h"),
            ("un_slab1_h = 100.0;", "un_slab1_h = %f;", "slab1_h"),
            ("un_slab2_h = 100.0;", "un_slab2_h = %f;", "slab2_h"),
            ("un_slab3_h = 100.0;", "un_slab3_h = %f;", "slab3_h"),
            ("un_slab4_h = 100.0;", "un_slab4_h = %f;", "slab4_h"),
            ("un_slab5_h = 100.0;", "un_slab5_h = %f;", "slab5_h"),
            ("lc_bkg = 0.05;", "lc_bkg = %f;", "lc_bkg"),
            ("lc_refine_1 = lc_bkg/2.0;", "lc_refine_1 = lc_bkg/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc_bkg/2.0;", "lc_refine_2 = lc_bkg/%f;", "lc_refine_2"),
            ("nslabs = 4;", "nslabs = %d;", "nslabs"),
        ]

        return subs

    def draw_mpl_frame(self, ax, styles):
        EDGE_COLOR = styles['edgecolor']
        LINEWIDTH = styles['linewidth']


        rib_w = self.get_param("rib_w") * nmtoum
        rib_h = self.get_param("rib_h") * nmtoum
        slab1_h = self.get_param("slab_w", 0) * nmtoum
        slab2_h = self.get_param("slab_h", 0) * nmtoum
        slab3_h = self.get_param("slab_h", 0) * nmtoum
        slab4_h = self.get_param("slab_h", 0) * nmtoum
        slab5_h = self.get_param("slab_h", 0) * nmtoum

        domx = self.get_param("domain_x") * nmtoum
        domy = self.get_param("domain_x") * nmtoum


        vertices = np.array(
            [
                [-domx / 2, -slab1_h],
                [-domx / 2, 0],
                [-rib_w / 2, 0],
                [-rib_w / 2, rib_h],
                [rib_w / 2, rib_h],
                [rib_w / 2, 0],
                [domx / 2, 0],
                #[slab_w / 2, -slab_h],
                #[-slab_w / 2, -slab_h],
            ]
        )

        ax.add_patch(
            mplpatches.Polygon(
                vertices, facecolor=None, fill=False, edgecolor=EDGE_COLOR, linewidth=LINEWIDTH
            )
        )



class Slot(UserGeometryBase):
    """NumBAT geometry template for a slot waveguide."""

    def init_geometry(self):

        self.set_properties("slot")

        self.set_required_parameters(
            ["rib_w", "rib_h", "slab_w", "slab_h", "slot_w"], num_req_mats=4
        )
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])
        self.set_parameter_help(
            {
                "rib_w": "width of raised ribs",
                "rib_h": "height of raised ribs",
                "slot_w": "width of slot between ribs",
                "slab_w": "width of slab substrate region",
                "slab_h": "height of slab substrate region",
                "material_bkg": "background material",
                "material_a": "material of slot",
                "material_b": "material of substrate slab",
                "material_c": "material of ribs",
                "lc_refine_1": "refine factor for slot and ribs",
                "lc_refine_2": "refine factor for slab and ribs",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 1000;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 500;", "dy_in_nm = %f;", "domain_y"),
            ("slot_w = 200;", "slot_w = %f;", "slot_w"),
            ("rib_h = 100;", "rib_h = %f;", "rib_h"),
            ("rib_w = 200;", "rib_w = %f;", "rib_w"),
            ("slab_w = 800;", "slab_w = %f;", "slab_w"),
            ("slab_h = 100;", "slab_h = %f;", "slab_h"),
            ("lc_bkg = 0.1;", "lc_bkg = %f;", "lc_bkg"),
            ("lc_refine_1 = lc_bkg/1;", "lc_refine_1 = lc_bkg/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc_bkg/1;", "lc_refine_2 = lc_bkg/%f;", "lc_refine_2"),
        ]

        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        slot_w = self.get_param("slot_w")
        rib_w = self.get_param("rib_w")
        rib_h = self.get_param("rib_h")
        slab_w = self.get_param("slab_w")
        slab_h = self.get_param("slab_h")

        msg = ""

        if slab_w >= dom_x:
            msg += "Slab width (slab_w) is larger than domain width (domain_x).\n"
        if slab_h >= dom_y / 2:
            msg += "Slab height (slab_h) is larger than half the domain height (domain_y).\n"
        if rib_h >= dom_y / 2:
            msg += (
                "Rib height (rib_h) is larger than half the domain height (domain_y).\n"
            )

        if slot_w + 2 * rib_w >= dom_x:
            msg += "Slot and ribs are together wider than domain width (domain_x).\n"

        dims_ok = not len(msg)

        return dims_ok, msg


class SlotCoated(UserGeometryBase):
    """NumBAT geometry template for a slot waveguide."""

    def init_geometry(self):

        self.set_properties("slot_coated")

        self.set_required_parameters(
            ["rib_w", "rib_h", "slab_w", "slab_h", "slot_w", "coat_t"], num_req_mats=5
        )
        self.set_allowed_parameters(
            ["lc_bkg", "lc_refine_1", "lc_refine_2", "lc_refine_3"]
        )
        self.set_parameter_help(
            {
                "rib_w": "width of raised ribs",
                "rib_h": "height of raised ribs",
                "slot_w": "width of slot between ribs",
                "slab_w": "width of slab substrate region",
                "slab_h": "height of slab substrate region",
                "coat_t": "thickness of coating layer",
                "material_bkg": "background material",
                "material_a": "material of slot",
                "material_b": "material of substrate slab",
                "material_c": "material of ribs",
                "material_d": "material of coating",
                "lc_refine_1": "refine factor for slot and ribs",
                "lc_refine_2": "refine factor for slab",
                "lc_refine_3": "refine factor for coating",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 1000;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 500;", "dy_in_nm = %f;", "domain_y"),
            ("slot_w = 200;", "slot_w = %f;", "slot_w"),
            ("rib_h = 100;", "rib_h = %f;", "rib_h"),
            ("rib_w = 200;", "rib_w = %f;", "rib_w"),
            ("slab_w = 800;", "slab_w = %f;", "slab_w"),
            ("slab_h = 100;", "slab_h = %f;", "slab_h"),
            ("coat_y = 100;", "coat_y = %f;", "coat_t"),
            ("lc_bkg = 0.1;", "lc_bkg = %f;", "lc_bkg"),
            ("lc_refine_1 = lc_bkg/1;", "lc_refine_1 = lc_bkg/%f;", "lc_refine_1"),
            ("lc_refine_2 = lc_bkg/1;", "lc_refine_2 = lc_bkg/%f;", "lc_refine_2"),
            ("lc_refine_3 = lc_bkg/1;", "lc_refine_3 = lc_bkg/%f;", "lc_refine_3"),
        ]

        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        slot_w = self.get_param("slot_w")
        rib_w = self.get_param("rib_w")
        rib_h = self.get_param("rib_h")
        slab_w = self.get_param("slab_w")
        slab_h = self.get_param("slab_h")
        coat_t = self.get_param("coat_t")

        msg = ""

        if slab_w >= dom_x:
            msg += "Slab width (slab_w) is larger than domain width (domain_x).\n"
        if slab_h >= dom_y / 2:
            msg += "Slab height (slab_h) is larger than half the domain height (domain_y).\n"
        if rib_h >= dom_y / 2:
            msg += "Rib height (rib_h) is larger than half the domain height (domain_y).\n"

        if rib_h + coat_t >= dom_y / 2:
            msg += "Rib height (rib_h) + coat thickness (coat_t) are together larger than half the domain height (domain_y).\n"

        if slot_w + 2 * rib_w >= dom_x:
            msg += "Slot and ribs are together wider than domain width (domain_x).\n"

        dims_ok = not len(msg)

        return dims_ok, msg
