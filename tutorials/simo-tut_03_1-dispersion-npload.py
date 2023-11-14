""" Calculate a dispersion diagram of the acoustic modes
    from q_AC ~ 0 (forward SBS) to q_AC = 2*k_EM (backward SBS).
    Load EM mode data from simo_tut_02.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../backend/")
import numbat
import integration
import mode_calcs
import objects
import materials

import starter

# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 2.5*lambda_nm
unitcell_y = unitcell_x
inc_a_x = 300
inc_a_y = 280
inc_shape = 'rectangular'
# Choose modes to include.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 25
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(3, sys.argv, 'a')

nbapp = numbat.NumBATApp()

wguide = objects.Structure(unitcell_x, inc_a_x, unitcell_y, inc_a_y, inc_shape,
                           material_bkg=materials.make_material("Vacuum"),
                           material_a=materials.make_material("Si_2016_Smith"),
                           lc_bkg=.1, lc_refine_1=4.0*refine_fac, lc_refine_2=4.0*refine_fac)

# wguide.check_mesh()
# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Assuming this calculation is run directly after simo-tut_02
# we don't need to recalculate EM modes, but can load them in.
sim_EM_pump = mode_calcs.load_simulation('tut02_wguide_data')
sim_EM_Stokes = mode_calcs.load_simulation('tut02_wguide_data2')

# Will scan from forward to backward SBS so need to know q_AC of backward SBS.
q_AC = np.real(sim_EM_pump.kz_EM(0) - sim_EM_Stokes.kz_EM(0))

# Number of wavevector steps.
nu_ks = 20

fig, ax = plt.subplots()
for i_ac, q_ac in enumerate(np.linspace(0.0, q_AC, nu_ks)):
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_ac, EM_sim=sim_EM_pump)
    prop_AC_modes = np.array(
        [np.real(x) for x in sim_AC.nu_AC_all() if abs(np.real(x)) > abs(np.imag(x))])
    sym_list = integration.symmetries(sim_AC)

    for i in range(len(prop_AC_modes)):
        Om = prop_AC_modes[i]*1e-9
        if sym_list[i][0] == 1 and sym_list[i][1] == 1 and sym_list[i][2] == 1:
            sym_A, = ax.plot(np.real(q_ac/q_AC), Om, 'or')
        if sym_list[i][0] == -1 and sym_list[i][1] == 1 and sym_list[i][2] == -1:
            sym_B1, = ax.plot(np.real(q_ac/q_AC), Om, 'vc')
        if sym_list[i][0] == 1 and sym_list[i][1] == -1 and sym_list[i][2] == -1:
            sym_B2, = ax.plot(np.real(q_ac/q_AC), Om, 'sb')
        if sym_list[i][0] == -1 and sym_list[i][1] == -1 and sym_list[i][2] == 1:
            sym_B3, = ax.plot(np.real(q_ac/q_AC), Om, '^g')

    print("Wavevector loop", i_ac+1, "/", nu_ks)

ax.set_ylim(0, 25)
ax.set_xlim(0, 1)
ax.legend([sym_A, sym_B1, sym_B2, sym_B3], [
          'A', r'B$_1$', r'B$_2$', r'B$_3$'], loc='lower right')
ax.set_xlabel(r'Normalised axial wavevector $q/(2\beta)$')
ax.set_ylabel(r'Frequency $\Omega/(2\pi)$ [GHz]')
fig.savefig('tut_03a-dispersion_symmetrised.png', bbox_inches='tight')

print(nbapp.final_report())
