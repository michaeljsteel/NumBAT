""" Calculate a dispersion diagram of the acoustic modes
    from q_AC ~ 0 (forward SBS) to q_AC = 2*k_EM (backward SBS).
    Load EM mode data from simo_tut_02.
"""

import sys
import math

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
sys.path.append(str(Path('../backend')))
import numbat
import materials
import mode_calcs
import integration


import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 2000
unitcell_y = 1000
inc_a_x = 450
inc_a_y = 200
inc_shape = "rectangular"
# Choose modes to include.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 25
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = "All"

prefix, refine_fac = starter.read_args(11, sys.argv)

nbapp = numbat.NumBATApp(prefix)

mat_core = materials.make_material("Si_2021_Poulton")
wguide = nbapp.make_structure(
    inc_shape,
    unitcell_x,
    unitcell_y,
    inc_a_x,
    inc_a_y,
    material_bkg=materials.make_material("Vacuum"),
    material_a=mat_core,
    lc_bkg=0.05,
    lc_refine_1=2.0 * refine_fac,
    lc_refine_2=2.0 * refine_fac,
)

# wguide.check_mesh()

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material("a").refindex_n - 0.1

sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

# Will scan from forward to backward SBS so need to know q_AC of backward SBS.
q_AC = np.real(sim_EM_pump.kz_EM(0) - sim_EM_Stokes.kz_EM(0))

# Number of wavevector steps.
nu_ks = 40

fig, ax = plt.subplots()
for i_ac, q_ac in enumerate(np.linspace(0.0, q_AC, nu_ks)):
    shift_Hz = q_ac * mat_core.Vac_shear() / (2 * math.pi)
    sim_AC = wguide.calc_AC_modes(
        num_modes_AC, q_ac, EM_sim=sim_EM_pump, shift_Hz=shift_Hz
    )
    prop_AC_modes = np.array(
        [np.real(x) for x in sim_AC.nu_AC_all() if abs(np.real(x)) > abs(np.imag(x))]
    )
    sym_list = integration.symmetries(sim_AC)

    for i in range(len(prop_AC_modes)):
        Om = prop_AC_modes[i] * 1e-9
        if sym_list[i][0] == 1 and sym_list[i][1] == 1 and sym_list[i][2] == 1:
            (sym_A,) = ax.plot(np.real(q_ac / q_AC), Om, "or")
        if sym_list[i][0] == -1 and sym_list[i][1] == 1 and sym_list[i][2] == -1:
            (sym_B1,) = ax.plot(np.real(q_ac / q_AC), Om, "vc")
        if sym_list[i][0] == 1 and sym_list[i][1] == -1 and sym_list[i][2] == -1:
            (sym_B2,) = ax.plot(np.real(q_ac / q_AC), Om, "sb")
        if sym_list[i][0] == -1 and sym_list[i][1] == -1 and sym_list[i][2] == 1:
            (sym_B3,) = ax.plot(np.real(q_ac / q_AC), Om, "^g")

    print("Wavevector loop", i_ac + 1, "/", nu_ks)

ax.set_ylim(0, 15)
ax.set_xlim(0, 1)
ax.legend(
    [sym_A, sym_B1, sym_B2, sym_B3],
    ["A", r"B$_1$", r"B$_2$", r"B$_3$"],
    loc="lower right",
)
ax.set_xlabel(r"Normalised axial wavevector $q/(2\beta)$")
ax.set_ylabel(r"Frequency $\nu = \Omega/(2\pi)$ [GHz]")
fig.savefig(prefix + "-disp-qnu.png", bbox_inches="tight")

print(nbapp.final_report())
