# Deprecated implementation example for file-based ``wg_impl`` waveguides.
# Prefer examples/user_waveguides and numbat.register_waveguide_directory().

from usermesh import UserGeometryBase

class MyGuide(UserGeometryBase):

    def init_geometry(self):
        self.set_properties("myguide")
        self.set_required_parameters([], 2)
        self._gmsh_template_filename = "oneincl"

    def apply_parameters(self):
        return [
            ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
            ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
        ]

    def draw_mpl_frame(self, ax, styles):
        pass
