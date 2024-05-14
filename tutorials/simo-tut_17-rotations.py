
import sys
import numpy as np

sys.path.append("../backend/")
import numbat
import materials
import mode_calcs
import integration
import plotting
from math import pi
import copy

import starter

# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 5*lambda_nm
unitcell_y = unitcell_x
inc_a_x = 550
inc_a_y = inc_a_x
inc_shape = 'circular'

mat_a=materials.make_material("SiO2_2016_Smith")        # isotropic
mat_b=materials.make_material("Si_2012_Rakich")  # cubic
mat_c=materials.make_material("LiNbO3aniso_2021_Steel") # tetragonal

pref = 'tut_17-'

mat_as = [copy.copy(mat_a) for i in range(4)]
mat_bs = [copy.copy(mat_b) for i in range(4)]
mat_cs = [copy.copy(mat_c) for i in range(4)]
print(mat_as)

mat_as[0].rotate(pi/2, 'x-axis')
mat_as[1].rotate(-pi/2, 'y-axis')
mat_as[2].rotate(pi/4, (1,1,1))
mat_as[3].rotate(pi/4, np.array([1,0,1]))

mat_bs[0].rotate(pi/2, 'x-axis')
mat_bs[1].rotate(-pi/2, 'y-axis')
mat_bs[2].rotate(pi/4, (1,1,1))
mat_bs[3].rotate(pi/4, np.array([1,1,1]))

mat_cs[0].rotate(pi/2, 'x-axis')
mat_cs[1].rotate(-pi/2, 'y-axis')
mat_cs[2].rotate(pi/4, (1,1,1))
mat_cs[3].rotate(pi/4, np.array([1,2,1]))


#for mats in (mat_as, mat_bs, mat_cs):
for mats in (mat_as,):
    for im,m in enumerate(mats): m.make_crystal_plot(pref+f'_{im:d}')
    #for m in mats: m.plot_bulk_dispersion(pref)


