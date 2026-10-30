import matplotlib.patches as mplpatches

from numbat_core.mesh.geometry import UserGeometryBase

from ._common import _process_one_and_two_incls, _process_one_and_two_incls_subs, nmtoum


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
