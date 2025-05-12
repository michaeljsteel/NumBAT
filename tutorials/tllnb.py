import sys
import numpy as np

sys.path.append("../backend")
import materials


mat_LiNbO3 = materials.make_material("LiNbO3aniso_2021_Steel")
print(mat_LiNbO3.full_str())

#mat_LiNbO3.plot_photoelastic_IJ("tt", ("11","12", "13", "14", "15"))

mat_LiNbO3.plot_bulk_dispersion('tt')

mat_LiNbO3.plot_bulk_dispersion_3D('tt')

materials.compare_bulk_dispersion(mat_LiNbO3, mat_LiNbO3, 'tt')
