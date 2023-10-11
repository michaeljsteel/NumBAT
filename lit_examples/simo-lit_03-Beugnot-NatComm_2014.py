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

def make_gain_plot(m_gain, extents, pref):
    L=0.08  # 8 cm nanowire
    Pp=0.001
    m_gaindB = 10*np.log10(np.abs(m_gain)*L*Pp + 1e-300)


    fig, ax1 = plt.subplots()
    im = ax1.imshow(m_gaindB.T, aspect='auto', interpolation='none',
                extent=extents, cmap='jet', origin='lower')
    for axis in ['top', 'bottom','left','right']:
     ax1.spines[axis].set_linewidth(1)

    cb=fig.colorbar(im, ax=ax1)
    cb.outline.set_visible(False)

    ax1.set_ylabel('Microwire diameter (μm)')
    ax1.set_xlabel('Acoustic frequency (GHz)')
    fig.savefig(pref+'-diam_scan.png')



start = time.time()

# Select the number of CPUs to use in simulation.
num_cores = 4

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

freq_min = 5e9
freq_max = 11e9

diam_0  =  1000
diam_min = 950
diam_max = 1300
diam_steps = 50
# find a d_diam that will give v_diams hitting 1000 exactly, with about 20 diam_steps total across (diam_max-diam_min)
d_diam = (diam_0-diam_min)/round((diam_0-diam_min)/((diam_max-diam_min)/diam_steps))
v_diams = np.arange(diam_min, diam_max, d_diam).round().astype(int) # make sure that 1000 nm is included.
num_diams = len(v_diams)
num_interp_pts = 2000

prefix, refine_fac = starter.read_args(3, sys.argv)

def modes_n_gain(diam):
    print('Handling diam', diam)
    inc_a_x = diam
    inc_a_y = diam
    # Use all specified parameters to create a waveguide object.
    wguide = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                            material_bkg=materials.make_material("Vacuum"),
                            material_a=materials.make_material("SiO2_2016_Smith"),
                            #lc_bkg=.1, lc_refine_1=4.0, lc_refine_2=4.0)
                            lc_bkg=.1, lc_refine_1=2.0, lc_refine_2=2.0)
    #wguide.plot_mesh(prefix+'_%3d'%int(inc_a_x))

    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)

    #plotting.plot_mode_fields(sim_EM_pump, ivals=range(5), EM_AC='EM_E', prefix=prefix)

    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
    q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))
    shift_Hz = 4e9
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)



    set_q_factor = 600.
    gainbox = integration.get_gains_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)#, fixed_Q=set_q_factor)

    gainbox.set_allowed_EM_pumps(EM_ival_pump)
    gainbox.set_allowed_EM_Stokes(EM_ival_Stokes)
    (nu_gain_tot, nu_gain_PE, nu_gain_MB) = gainbox.plot_spectra(freq_min, freq_max, 
                                                                 num_interp_pts=num_interp_pts, semilogy=True, dB=True, 
                                                                 prefix = prefix, suffix='_w%i' %int(inc_a_x))

    if abs(diam - 1000)<1e-10:  # print fields for 1 micron guide
        plotting.plot_mode_fields(sim_EM_pump, EM_AC = 'EM_E', ivals=range(5), prefix=prefix+'-diam-1000')
        plotting.plot_mode_fields(sim_AC, ivals=range(40), prefix=prefix+'-diam-1000')
        for m in range(num_modes_AC):
            print(f'{m}, {sim_AC.nu_AC(m)*1e-9:.4f}, {gainbox.gain_total(m):.3e}, ',
                 f'{gainbox.gain_PE(m):.4e}, {gainbox.gain_MB(m):.4e}')


    return (nu_gain_tot, nu_gain_PE, nu_gain_MB) 

# Run diams in parallel across num_cores CPUs using multiprocessing package.
pool = Pool(num_cores)
gain_specs = pool.map(modes_n_gain, v_diams)


m_gain_tot= np.zeros((num_interp_pts, num_diams))
m_gain_PE= np.zeros((num_interp_pts, num_diams))
m_gain_MB= np.zeros((num_interp_pts, num_diams))
m_gaindB = np.zeros((num_interp_pts, num_diams))

for idiam, gains in enumerate(gain_specs):
    (nu_gain_tot, nu_gain_PE, nu_gain_MB)  = gains
    m_gain_tot[:,idiam] = nu_gain_tot
    m_gain_PE[:,idiam] = nu_gain_PE
    m_gain_MB[:,idiam] = nu_gain_MB 

# label axes in GHz and microns
extents=(freq_min*1e-9, freq_max*1e-9, v_diams[0]*1e-3, v_diams[-1]*1e-3)
make_gain_plot(m_gain_tot, extents, prefix+'-gain_tot')
make_gain_plot(m_gain_PE, extents, prefix+'-gain_PE')
make_gain_plot(m_gain_MB, extents, prefix+'-gain_MB')

end = time.time()
print("\n Simulation time (sec.)", (end - start))

print("--------------------------------------------------------------------\n\n\n")
