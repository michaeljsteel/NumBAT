from numbat_core.mesh.geometry import UserGeometryBase


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
