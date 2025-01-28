# Copyright (C) 2017-2025  Michael Steel.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


import uuid
import os
from pathlib import Path

import reporting
from numbattools import is_real_number

class UserGeometryBase():

    def __init__(self, params, d_materials):
        self._shape_name = ''
        self._num_materials=0          # actual number of materials defined
        self._num_req_materials = 0    # number of materials required by the template
        self.is_curvilinear=False
        self._d_materials = d_materials
        self._d_params = params
        self._descrip = 'Unimplemented user template'
        self._gmsh_template_filename = ''  # internal meshes use this. User meshes should not

        self._req_params =[]
        self._allowed_params =[]
        self.d_param_help = {}

    def set_properties(self, nm, is_curvi=False):
        self._shape_name = nm
        self.is_curvilinear = is_curvi
        self._descrip = self.__doc__

    def set_required_parameters(self, l_nms, num_req_mats):
        self._req_params.extend(['domain_x','domain_y', 'lc_bkg'])
        self._req_params.extend(l_nms)
        self._num_req_materials = num_req_mats
        self._num_materials = max(self._num_materials, num_req_mats)

    def set_allowed_parameters(self, l_nms, num_allowed_mats=0):
        self._allowed_params.extend(l_nms)

        num_allowed_mats = max(num_allowed_mats, self._num_req_materials)
        self._num_materials = max(self._num_materials, num_allowed_mats)


        for im in range(num_allowed_mats): # need mat_a, mat_b, mat_c, etc
            self._allowed_params.append('material_'+'abcdefghijklmnopqrtstuvwxyz'[im])

    def set_parameter_help(self, d_help):
        self.d_param_help = {
            'domain_x': "length of simulation domain along x",
            'domain_y': "length of simulation domain along y"}
        self.d_param_help.update(d_help)

    def get_parameter_help_summary(self):
        head = f'Waveguide parameters for shape {self._shape_name}:\n'
        #for k,v in self.d_param_help.items():
        #    msg += f'{k:>20} : {v}\n'
        body = '\n'.join([f'{k:>20} : {v}' for k,v in self.d_param_help.items()])

        msg = head + body
        return msg

    def check_parameters(self, user_params):
        if not self._req_params: # not yet defined for this template
            return

        reqkws = self._req_params
        reqkws.append('material_bkg')

        for im in range(self._num_req_materials-1): # need mat_a, mat_b, mat_c, etc, -1 because one is mat_bkg
            reqkws.append('material_'+'abcdefghijklmnopqrtstuvwxyz'[im])

        for key in reqkws:
            if key not in user_params:
                msg =(f"Waveguide type '{self._shape_name}' requires a value for the parameter '{key}' in the call to make_structure()."
                      '\n\nNote that some waveguide types have changed their required parameters to adopt more intuitive names.'
                      '\n\nFor this waveguide type, the following guidelines apply:\n\n')

                msg+=self.get_parameter_help_summary()

                reporting.report_and_exit(msg)

        # report unexpected keys
        goodkeys = reqkws + self._allowed_params
        goodkeys.append('lc') # remove special case once everything is moved to lc_bkg
        for key in user_params.keys():
            if key not in goodkeys:
                reporting.report(
                    f"Waveguide '{self._shape_name}' will ignore the parameter '{key}' in make_structure().")


    def check_dimensions(self):  # override to implement check for each mesh design
        dims_ok = True
        msg=''
        return dims_ok, msg

    def validate_dimensions(self):
        '''Checks that the combination of user parameters defines a well-defined consistent geometry.'''

        dims_ok, msg = self.check_dimensions()

        if not dims_ok: reporting.report_and_exit(f'There is a problem with the waveguide structure:\n{msg}')

    def get_param(self, k, dflt=None):
        return self._d_params.get(k, dflt)

    def gmsh_template_filename(self):
        if self._gmsh_template_filename:
            return self._gmsh_template_filename
        else:
            return self._shape_name

    def __str__(self):
        return self._descrip

    def make_geometry(self, p_dir_templates):
        subs = self.apply_parameters()
        #if subs is None:
        #    subs = self._param_subs

        geo = open(Path(p_dir_templates,
                   f'{self.gmsh_template_filename()}_msh_template.geo'), 'r').read()

        for (olds, news, sval) in subs:
            # sval can be either an actual float or a string rep of a float
            # or the name of a parameter
            try:
                val = float(sval)
            except ValueError:
                val = self.get_param(sval)

            if val is not None:
                assert is_real_number(val), f'Parameter {sval} is not a number'
                geo = geo.replace(olds, news % val)

        return geo


    def get_instance_filename(self):   # , l_dims):
        '''Make name for the concrete instantiation of a given mesh .geo template'''
        msh_fname = self._shape_name

        # made crazy long names, not helping
        # for v in l_dims:
        #    if is_real_number(v): msh_name += '_%s' % dec_float_str(v)

        msh_fname += f'--pid-{os.getpid()}'

        # need to make name unique to support parallel processing
        msh_fname += '--'+str(uuid.uuid4())

        return msh_fname



    # Methods to be overridden by derived classes

    def init_geometry(self):
        print('IMPLEMENT ME', __file__, 'make_geometry')
        raise NotImplementedError('UserGeometryBase.init_geometry() must be overriden')
        return ''

    def apply_parameters(self):
        print('IMPLEMENT ME', __file__, 'make_geometry')
        raise NotImplementedError('UserGeometryBase.apply_parameters() must be overriden')
        return ''

    def draw_mpl_frame(self, ax, styles): # optional method
        '''Add elements to a matplotlib axis to draw outline of the structure'''
        pass
        # print('base class drawmplframe on ', ax)
