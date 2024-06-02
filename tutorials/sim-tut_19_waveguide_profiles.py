
import sys

sys.path.append("../backend/")
import numbat
import materials

import starter


# Geometric Parameters - all in nm.
#lambda_nm = 1550
lambda_nm = 1000
unitcell_x = lambda_nm
unitcell_y = unitcell_x
inc_a_x = 600      # base length (always horizontal)
inc_a_y = inc_a_x  # not used
inc_b_x = 500      # displacement of peak from left end of base
inc_b_y = 500      # height of peak from base

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
lc_norm = 3
lc_corner = 6

lc_norm = 1
lc_corner = 1

wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                              inc_b_x = inc_b_x, inc_b_y = inc_b_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("Si_2021_Poulton"),
                        lc_bkg=lc_bkg, lc_refine_1=lc_norm, lc_refine_2=lc_corner)



wguide.plot_refractive_index_profile(prefix, as_epsilon=True)

#wguide.plot_phase_velocity_z_profile(prefix)
