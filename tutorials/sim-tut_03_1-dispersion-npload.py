""" Calculate a dispersion diagram of the acoustic modes
    from q_AC ~ 0 (forward SBS) to q_AC = 2*k_EM (backward SBS).
    Load EM mode data from simo_tut_02.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
import integration

import materials
from nbtypes import SI_GHz

import starter

# Geometric Parameters - all in nm.
lambda_nm = 1550
domain_x = 2.5*lambda_nm
domain_y = domain_x
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

nbapp = numbat.NumBATApp(prefix)

wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                           material_bkg=materials.make_material("Vacuum"),
                           material_a=materials.make_material("Si_2016_Smith"),
                           lc_bkg=.1, lc_refine_1=4.0*refine_fac, lc_refine_2=4.0*refine_fac)

# wguide.check_mesh()
# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Assuming this calculation is run directly after sim-tut_02
# we don't need to reuse_old_fieldsulate EM modes, but can load them in.
reuse_old_fields = False
if reuse_old_fields:
    simres_EM_pump = numbat.load_simulation('tut02_em_pump')
    simres_EM_Stokes = numbat.load_simulation('tut02_em_stokes')
else:
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    simres_EM_Stokes = simres_EM_pump.bkwd_Stokes_modes()



# Will scan from forward to backward SBS so need to know q_AC of backward SBS.
q_AC = np.real(simres_EM_pump.kz_EM(0) - simres_EM_Stokes.kz_EM(0))

# Number of wavevector steps.
nu_qs = 50

fig, ax = plt.subplots()
symmetries_working = False
for i_ac, q_ac in enumerate(np.linspace(0.0, q_AC, nu_qs)):
    print("Wavevector loop", i_ac+1, "/", nu_qs)
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_ac, EM_sim=simres_EM_pump)
    v_Nu = np.array(
        [np.real(x) for x in sim_AC.nu_AC_all() if abs(np.real(x)) > abs(np.imag(x))])

    if symmetries_working:
        sym_list = integration.symmetries(sim_AC)
        for i in range(len(prop_AC_modes)):
            Nu = v_Nu[i]/SI_GHz
            if sym_list[i][0] == 1 and sym_list[i][1] == 1 and sym_list[i][2] == 1:
                sym_A, = ax.plot(np.real(q_ac/q_AC), Nu, 'or', markersize=2)
            if sym_list[i][0] == -1 and sym_list[i][1] == 1 and sym_list[i][2] == -1:
                sym_B1, = ax.plot(np.real(q_ac/q_AC), Nu, 'vc', markersize=2)
            if sym_list[i][0] == 1 and sym_list[i][1] == -1 and sym_list[i][2] == -1:
                sym_B2, = ax.plot(np.real(q_ac/q_AC), Nu, 'sb', markersize=2)
            if sym_list[i][0] == -1 and sym_list[i][1] == -1 and sym_list[i][2] == 1:
                sym_B3, = ax.plot(np.real(q_ac/q_AC), Nu, '^g', markersize=2)

    else:
        ax.plot(np.real(q_ac/q_AC)*np.ones(len(v_Nu)), v_Nu/SI_GHz, 'o', markersize=2)


ax.set_ylim(0, 25)
ax.set_xlim(0, 1)
if symmetries_working:
    ax.legend([sym_A, sym_B1, sym_B2, sym_B3], [
          'A', r'B$_1$', r'B$_2$', r'B$_3$'], loc='lower right')
ax.set_xlabel(r'Normalised axial wavevector $q/(2\beta)$')
ax.set_ylabel(r'Frequency $\Omega/(2\pi)$ [GHz]')
fig.savefig('tut_03a-dispersion_symmetrised.png', bbox_inches='tight')

print(nbapp.final_report())
