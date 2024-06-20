"""
Calculate dispersion diagram of the acoustic modes in a rectangular Si waveguide
"""

# Import the necessary packages



import sys
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))
import numbat
import materials
import mode_calcs
import integration


import starter



# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 3.01*wl_nm
unitcell_y = unitcell_x
inc_a_x = 450 # Waveguide widths.
inc_a_y = 200
inc_shape = 'rectangular'
# Choose modes to include.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 100
EM_ival_pump = 0
EM_ival_Stokes = EM_ival_pump
AC_ival = 'All'

prefix, refine_fac = starter.read_args(2, sys.argv, sub='b')

nbapp = numbat.NumBATApp(prefix)

# Use all specified parameters to create a waveguide object
wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("Si_2021_Poulton"),
                        lc_bkg=0.05, # mesh coarseness in background, larger lc_bkg = coarser along horizontal outer edge
                        lc_refine_1=4*refine_fac, # mesh refinement factor near the interface of waveguide, larger = finer along horizontal interface
                        lc_refine_2=5*refine_fac) # mesh refinement factor near the origin/centre of waveguide

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1
# Calculate Electromagnetic modes.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

# Print EM mode info
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print(f'{i:3d}  {np.real(kz):.4e}')

n_eff_sim = np.real(sim_EM_pump.neff(EM_ival_pump))
print("\n Fundamental optical mode ")
print(" n_eff = ", np.round(n_eff_sim, 4))

#   q_AC of backward SBS.
q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))
# Number of wavevectors steps.
nu_ks = 50

plt.clf()
plt.figure(figsize=(10,6))
ax = plt.subplot(1,1,1)
symmetries_working = False
for i_ac, q_ac in enumerate(np.linspace(0.0,q_AC,nu_ks)):
    print("Wavevector loop", i_ac+1, "/", nu_ks)

    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_ac, EM_sim=sim_EM_pump)
    prop_AC_modes = np.array([np.real(x) for x in sim_AC.eigs_nu if abs(np.real(x)) > abs(np.imag(x))])

    if symmetries_working:
        sym_list = integration.symmetries(sim_AC)   # FIXME Symmetries
        for i in range(len(prop_AC_modes)):
            Om = prop_AC_modes[i]*1e-9
            if sym_list[i][0] == 1 and sym_list[i][1] == 1 and sym_list[i][2] == 1:
                sym_A, = plt.plot(np.real(q_ac/q_AC), Om, 'or')
            if sym_list[i][0] == -1 and sym_list[i][1] == 1 and sym_list[i][2] == -1:
                sym_B1, = plt.plot(np.real(q_ac/q_AC), Om, 'vc')
            if sym_list[i][0] == 1 and sym_list[i][1] == -1 and sym_list[i][2] == -1:
                sym_B2, = plt.plot(np.real(q_ac/q_AC), Om, 'sb')
            if sym_list[i][0] == -1 and sym_list[i][1] == -1 and sym_list[i][2] == 1:
                sym_B3, = plt.plot(np.real(q_ac/q_AC), Om, '^g')


ax.set_ylim(0,15)
ax.set_xlim(0,1)

if symmetries_working:
    plt.legend([sym_A, sym_B1, sym_B2, sym_B3],['A',r'B$_1$',r'B$_2$',r'B$_3$'], loc='lower right')

plt.xlabel(r'Axial wavevector (normalised)')
plt.ylabel(r'Frequency (GHz)')
plt.savefig(prefix+'-disp-qnu.png', bbox_inches='tight')
plt.close()

print(nbapp.final_report())
