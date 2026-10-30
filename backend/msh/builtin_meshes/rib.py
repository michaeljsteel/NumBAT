import matplotlib.patches as mplpatches
import numpy as np

from numbat_core import reporting
from usermesh import UserGeometryBase

from ._common import nmtoum


class Pedestal(UserGeometryBase):
    """NumBAT geometry template for a pedestal-type waveguide."""

    def init_geometry(self):
        self.set_properties("pedestal")

        self.set_required_parameters(
            [
                "ped_top_w",
                "ped_base_w",
                "ped_h",
                "pillar_x",
                "pillar_y",
                "slab_a_x",
                "slab_a_y",

            ],
            num_req_mats=4,
        )
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])

        self.set_parameter_help(
            {
                "ped_top_w": "width of the top of the pedestal",
                "ped_base_w": "width of the base of the pedestal",
                "ped_h": "height of the pedestal",
                "pillar_x": "width of the supporting pillar",
                "pillar_y": "height of the supporting pillar",
                "slab_a_x": "width of the slab",
                "slab_a_y": "height of the slab",

                "material_bkg": "material of background",
                "material_a": "material of rib",
                "material_b": "material of pillar",
                "material_c": "material of slab",
                "lc": "grid points arounds boundary as fraction of domain_x",
                "lc_refine1": "refined density along upper rib",
                "lc_refine2": "refined density along buried rib",
            }
        )

    def apply_parameters(self):

        subs = [
            ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
            ("a1top = 15;", "a1top = %f;", "ped_top_w"),
            ("a1 = 20;", "a1 = %f;", "ped_base_w"),
            ("a1y = 10;", "a1y = %f;", "ped_h"),
            ("slabx = 80;", "slabx = %f;", "slab_a_x"),
            ("slaby = 10;", "slaby = %f;", "slab_a_y"),
#            ("slabxtop = 60;", "slabxtop = %f;", "slab_b_x"),
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
                "slab_b_h",
                "coat_w",
                "coat_h",
                "coat2_w",
                "coat2_h",
            ],
            num_req_mats=6,
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
                "slab_w": "width of upper slab substrate region",
                "slab_h": "height of upper slab substrate region",
                "slab_b_h": "height of lower slab substrate region",
                "coat_w": "horizontal thickness of inner coating layer",
                "coat_h": "vertical thickness of inner coating layer",
                "coat2_w": "horizontal thickness of outer coating layer",
                "coat2_h": "vertical thickness of outer coating layer",
                "material_bkg": "background material",
                "material_a": "material of raised rib core",
                "material_b": "material of primary waveguide substrate layer (often the same as material_a)" ,
                "material_c": "material of inner coating layer",
                "material_d": "material of lower substrate layer",
                "material_e": "material of outer coating layer",
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
            ("slab2y = 5;", "slab2y = %f;", "slab_b_h"),
            ("coat2x = 4;", "coat2x = %f;", "coat2_w"),
            ("coat2y = 4;", "coat2y = %f;", "coat2_h"),
            ("lc_bkg = 0.05;", "lc_bkg = %f;", "lc_bkg"),
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
