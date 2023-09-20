""" Replicating the results of
    Brillouin light scattering from surface acoustic
    waves in a subwavelength-diameter optical fibre
    Beugnot et al.
    http://dx.doi.org/10.1038/ncomms6242
"""


import os
import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
import matplotlib
import matplotlib.pyplot as plt

sys.path.append("../backend/")
import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT

import starter


start = time.time()

# Select the number of CPUs to use in simulation.
num_cores = 1

# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 4000
unitcell_y = unitcell_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 80
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

# Expected effective index of fundamental guided mode.
n_eff = 1.18

freq_min = 4e9
freq_max = 12e9

#width_min = 600
#width_max = 1200
#num_widths = 301
width_min = 900
width_max = 1300
num_widths = 21
num_widths = 11

v_widths = np.linspace(width_min, width_max, num_widths)
num_interp_pts = 2000

prefix, refine_fac = starter.read_args(3, sys.argv)

def modes_n_gain(inc_a_x):
    inc_a_y = inc_a_x
    # Use all specified parameters to create a waveguide object.
    wguide = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                            material_bkg=materials.make_material("Vacuum"),
                            material_a=materials.make_material("SiO2_2016_Smith"),
                            lc_bkg=.1, lc_refine_1=4.0, lc_refine_2=4.0)
    #wguide.plot_mesh(prefix+'_%3d'%int(inc_a_x))

    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)

    #plotting.plot_mode_fields(sim_EM_pump, ivals=range(5), EM_AC='EM_E', prefix=prefix)

    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
    q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))
    shift_Hz = 4e9
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

    #plotting.plot_mode_fields(sim_AC, ivals=range(5), prefix=prefix)


    set_q_factor = 600.
    SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)#, fixed_Q=set_q_factor)
    interp_values = plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
                                               EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min, freq_max,
                                               num_interp_pts=num_interp_pts, semilogy=True, dB=True,
                                               prefix = prefix, suffix='_w%i' %int(inc_a_x))

    return interp_values

# Run widths in parallel across num_cores CPUs using multiprocessing package.
pool = Pool(num_cores)
width_objs = pool.map(modes_n_gain, v_widths)
# Note pool.map() doesn't pass errors back from fortran routines very well.
# It's good practise to run the extrema of your simulation range through map()
# before launcing full multicore simulation.


m_gain= np.zeros((num_interp_pts, num_widths))
m_gaindB = np.zeros((num_interp_pts, num_widths))
for w, width_interp in enumerate(width_objs):
    m_gain[:,w] = width_interp[::1]

L=0.08  # 8 cm nanowire
Pp=0.001
m_gaindB = 10*np.log10(m_gaindB*L*Pp + 1e-30)

fig, ax1 = plt.subplots()
im = ax1.imshow(m_gaindB.T, aspect='auto', interpolation='none',
                #vmin=-60, vmax=0,
                extent=[freq_min, freq_max, v_widths[0], v_widths[-1]],
                cmap='jet', origin='lower')

num_xticks = 5
num_yticks = 5
#ax1.xaxis.set_ticks_position('bottom')
#ax1.set_yticks(np.linspace(0,(num_widths-1),num_xticks))
#ax1.set_xticks(np.linspace((num_interp_pts-1),0,num_yticks))
#ax1.set_yticklabels(["%4.0f" % (w*0.001) for w in np.linspace(width_min,width_max,num_xticks)])
#ax1.set_xticklabels(["%4.0f" % (nu*1e-9) for nu in np.linspace(freq_min,freq_max,num_yticks)])
ax1.set_xlim(5,11)

ax1.set_ylabel(r'Microwire diameter (um)')
ax1.set_xlabel('Acoustic frequency (GHz)')
fig.savefig(prefix+'-gain_width_scan.png')

end = time.time()
print("\n Simulation time (sec.)", (end - start))

print("--------------------------------------------------------------------\n\n\n")
