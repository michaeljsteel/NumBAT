import uuid
import os
from pathlib import Path

def is_real_number(x):
    # return isinstance(x, float) or isinstance(x, int)
    # return isinstance(x, numbers.Number)
    return isinstance(x, (int, float))
    #try:
    #    xx = float(x)
    #    return True
    #except Exception:
    #    return False

class UserGeometryBase(object):

    def __init__(self, params, d_materials):
        self._geom_name = ''
        self._num_type_elements=0
        self._is_curvilinear=False
        self._d_materials = d_materials
        self._d_params = params
        self._descrip = 'Unimplemented user template'
        self._gmsh_template_filename = ''  # internal meshes use this. User meshes should not

    def get_param(self, k):
        return self._d_params.get(k, None)

    def set_num_type_elements(self, n):
        self._num_type_elements = n

    def set_is_curvilinear(self, b):
        self._is_curvilinear = b

    def set_name(self, nm):
        self._geom_name = nm

    def set_description(self, desc):
        self._descrip = desc

    def set_properties(self, nm, n_elements, is_curvi, desc):
        self.set_name(nm)
        self.set_num_type_elements(n_elements)
        self.set_is_curvilinear(is_curvi)
        self.set_description(desc)

    def geom_name(self):
        return self._geom_name

    def gmsh_template_filename(self):
        if self._gmsh_template_filename:
            return self._gmsh_template_filename
        else:
            return self.geom_name()

    def num_type_elements(self):
        return self._num_type_elements

    def is_curvilinear(self):
        return self._is_curvilinear

    def __str__(self):
        return self._descrip

    def apply_parameters(self):
        print('IMPLEMENT ME', __file__, 'make_geometry')
        return ''

    def make_geometry(self, p_dir_templates):
        subs = self.apply_parameters()

  #$      msh_template = self.wg_geom.gmsh_template_filename()
#$
    #        geo = self._load_mesh_template(msh_template)

        geo = open(Path(p_dir_templates,
                   f'{self.gmsh_template_filename()}_msh_template.geo'), 'r').read()

        for (olds, news, val) in subs:
            if val is None:  # unset value not overridden or dropped
                continue
            elif is_real_number(val):
                geo = geo.replace(olds, news % val)
            else:
                geo = geo.replace(olds, news)

        return geo


    def get_instance_filename(self):   # , l_dims):
        '''Make name for the concrete instantiation of a given mesh .geo template'''
        msh_fname = self._geom_name

        # made crazy long names, not helping
        # for v in l_dims:
        #    if is_real_number(v): msh_name += '_%s' % dec_float_str(v)

        msh_fname += '--pid-%d' % os.getpid()

        # need to make name unique to support parallel processing
        msh_fname += '--'+str(uuid.uuid4())

        return msh_fname

    def draw_mpl_frame(self, ax):
        '''Add elements to a matplotlib axis to draw outline of the structure'''
        print('base class drawmplframe')
        pass