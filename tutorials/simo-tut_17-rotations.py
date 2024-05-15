
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

np.set_printoptions( precision=4)
# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 5*lambda_nm
unitcell_y = unitcell_x
inc_a_x = 550
inc_a_y = inc_a_x
inc_shape = 'circular'

mat_a=materials.make_material("SiO2_2021_Poulton")        # isotropic
#mat_b=materials.make_material("Si_2012_Rakich")  # cubic
mat_b=materials.make_material("Ge_cubic_2014_Wolff")  # cubic
mat_c=materials.make_material("LiNbO3aniso_2021_Steel") # tetragonal

mat_d=materials.make_material("LiNbO3aniso_2021_Steel") # tetragonal

pref = 'tut_17'

mat_as = [copy.deepcopy(mat_a) for i in range(4)]
mat_bs = [copy.deepcopy(mat_b) for i in range(4)]
mat_cs = [copy.deepcopy(mat_c) for i in range(4)]
mat_ds = [copy.deepcopy(mat_d) for i in range(4)]

#print('x pi/4')
mat_as[0].rotate(pi/5, 'x-axis')
#print('\n\ny -pi/2')
#mat_as[1].rotate(-pi/2, 'y-axis')
#print('\n\n111, pi/4')
#mat_as[2].rotate(pi/4, (1,1,1))
#print('\n\n101, pi/4')
#mat_as[3].rotate(pi/4, np.array([1,0,1]))
#

mat_bs[0].rotate(pi/2, 'x-axis')
mat_bs[1].rotate(-pi/2, 'y-axis')
mat_bs[2].rotate(pi/4, (1,1,1))
mat_bs[3].rotate(pi/3, np.array([1,1,1]))

#mat_cs[0].rotate(pi/2, 'x-axis')
mat_cs[1].rotate(-pi/2, 'y-axis')
mat_cs[2].rotate(pi/4, (1,1,1))
mat_cs[3].rotate(pi/4, np.array([1,2,1]))

#mat_ds[0].set_orientation('x-cut')
#mat_ds[1].set_orientation('y-cut')
#mat_ds[2].set_orientation('y-cut')
#mat_ds[2].set_orientation('z-cut')
#mat_ds[3].set_orientation('y-cut')
#mat_ds[3].set_orientation('z-cut')
#mat_ds[3].set_orientation('x-cut')




#for mats in (mat_as, mat_bs, mat_cs, mat_ds):
##for mats in (mat_as,):
#    for im,m in enumerate(mats): 
#        m.make_crystal_plot(pref+f'{im:d}')
#    #for m in mats: m.plot_bulk_dispersion(pref)

print(mat_cs[0].c_tensor.mat)
for am,mats in enumerate((mat_cs,)):
    for im,m in enumerate((mats[0],)): 
        m.plot_bulk_dispersion(pref+f'-{am:d}-{im:d}')


