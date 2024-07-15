
import sys

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
import materials

import starter


# Geometric Parameters - all in nm.
#lambda_nm = 1550
lambda_nm = 1000
domain_x = lambda_nm * 1.5
domain_y = domain_x
tri_wid = 800      
tri_hgt = 500  
tri_xoff = 300   

inc_shape = 'triangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(19, sys.argv, sub='b')

nbapp = numbat.NumBATApp(prefix)

refine_fac = 1
lc_bkg = .1 * refine_fac
lc_norm = 30
lc_corner = 30

wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, 
                              base_width=tri_wid, peak_height=tri_hgt, peak_xoff=tri_xoff, 
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("Si_2021_Poulton"),
                        lc_bkg=lc_bkg, lc_refine_1=lc_norm, lc_refine_2=lc_corner)



wguide.plot_refractive_index_profile(prefix, as_epsilon=True)

#wguide.plot_phase_velocity_z_profile(prefix)
