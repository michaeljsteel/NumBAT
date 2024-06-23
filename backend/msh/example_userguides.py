# Implementation file for user waveguides
# The filename 'example_guides.py' matches the value of 'wg_impl' in user_waveguides.json.

from usermesh import UserGeometryBase

# Name of the derived class matches the value of 'wg_class' in user_waveguides.json.
class MyGuide(UserGeometryBase):

    # Set up any needed adjustable dimensions
    def init_geometry(self):
        desc = '''An example NumBAT geometry template for a user defined waveguide.  '''


        # First argument of this call matches the value of 'inc_shape' in user_waveguides.json.
        # Second argument is the number of distinct regions that exist in the template.
        # Corresponds to the maximum number of distinct materials that can be created for this geometry.
        self.set_properties('myguide', 2, False, desc)


    # Use information from the python call to make_structure to adjust any dimensions defined
    # in the .geo template file.

    def apply_parameters(self):

        # At a minimum we usually need to adjust the unit cell size
        subs = [('dx_in_nm = 100;', 'dx_in_nm = %f;', self.get_param('domain_x'))]
        subs.append(('dy_in_nm = 50;', 'dy_in_nm = %f;', self.get_param('domain_y')))

        return subs

    # This method is optional.
    # It consists of commands to draw the outline of the waveguide onto matplotlib plots of the mode profiles.

    def draw_mpl_frame(self, ax):
        pass

