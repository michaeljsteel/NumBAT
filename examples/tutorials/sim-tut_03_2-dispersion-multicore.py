""" Calculate a dispersion diagram of the acoustic modes
    from q_AC ~ 0 (forward SBS) to q_AC = 2*k_EM (backward SBS).
    Use python's (embarrassing parallel) multiprocessing package.
"""

import sys
import os
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import copy

from pathlib import Path
sys.path.append(str(Path('../../backend')))


import numbat

import materials

import starter





# Geometric Parameters - all in nm.
lambda_nm = 1550
#domain_x = 3.0*lambda_nm
#domain_y = domain_x
domain_x = 2.0*lambda_nm
domain_y = domain_x

inc_a_x = 800.
inc_a_y = 220.
inc_shape = 'rectangular'
# Choose modes to include.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 60
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(3, sys.argv, 'b')

nbapp = numbat.NumBATApp(prefix)

# Note that this mesh is quite fine, may not be required if purely using dispersive sims
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                           material_bkg=materials.make_material("Vacuum"),
                           material_a=materials.make_material("Si_2016_Smith"),
                           #lc_bkg=.1, lc_refine_1=5.0*refine_fac, lc_refine_2=5.0*refine_fac)
                           lc_bkg=.05, lc_refine_1=1.0*refine_fac, lc_refine_2=1.0*refine_fac)

#wguide.plot_mesh(prefix)

# Estimated effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic modes.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()

# Will scan from forward to backward SBS so need to know q_AC of backward SBS.
q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) -
               sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))

# Rather than calculating with a loop we can use pool to do a multi core sim


def solve_ac_mode_freqs(qset):
    ik, nk, q_AC, sim_EM = qset
    print('\nPID: %d  commencing mode calculation for q_AC %d/%d = %f /m' % (
        os.getpid(), ik+1, nk, q_AC))

    # Calculate the modes, grab the output frequencies only and convert to GHz
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM)

    print('PID: %d got ac modes for q_AC = %f' % (os.getpid(), q_AC))

    mode_freqs = np.real(sim_AC.nu_AC_all()) * 1.e-9  # convert to GHz

    # Clear memory
    sim_AC = None

    print('PID: %d completed mode calculation for width a_x = %f' % (
        os.getpid(), q_AC))

    # Return the frequencies and simulated q_AC value in a list
    return mode_freqs


# Now we utilise multi-core calculations to perform parallel simulations and speed up the simulation
n_qs = 100  # start with a low number of q_AC values to get an idea of the scale of the problem
acoustic_qs = np.linspace(5., q_AC*1.1, n_qs)

# Output the normalisation k value for reference
print(f"The acoustic wavevector 2*kp = {q_AC:f}")

multiproc = nbapp.is_linux()  # multiproc not yet working for windows and macos
#multiproc = False

# make jobs list with entries of form  (iq, n_qs, q)  (qstep, total qs, qval)
emsim_copies = [copy.deepcopy(sim_EM_pump) for i in range(n_qs)]
qsets = zip(np.arange(n_qs), np.arange(n_qs)*0+n_qs,
            acoustic_qs, emsim_copies)

if multiproc:
    num_cores = os.cpu_count()  # Let OS decide how many processes to run
    num_cores = 4
    with  Pool(num_cores) as pool:
        pooled_mode_freqs = pool.map(solve_ac_mode_freqs, qsets)
else:
    pooled_mode_freqs = []
    for qset in qsets:
        pooled_mode_freqs.append(solve_ac_mode_freqs(qset))


# We will pack the above values into a single array for plotting purposes, initialise first
m_nu_ac = np.zeros((n_qs, num_modes_AC+1))
for i_w, sim_freqs in enumerate(pooled_mode_freqs):
    # Set the value to the values in the frequency array
    m_nu_ac[i_w, 0] = acoustic_qs[i_w]
    m_nu_ac[i_w, 1:] = sim_freqs

# Now that we have packed will save to a numpy file for better plotting and reference
file_name = prefix+'-disp'
# mat=np.transpose([acoustic_qs, m_nu_ac])  # arrange the data as two columns of wavenumber and frequency
np.savetxt(file_name, m_nu_ac)

# Also plot a figure for reference
plot_range = num_modes_AC
fig, ax = plt.subplots()
for idx in range(plot_range):
    # slicing in the row direction for plotting purposes
    nu_slice = m_nu_ac[:, 1+idx]
    ax.plot(acoustic_qs/q_AC, nu_slice, '.', markersize=2)

# Set the limits and plot axis labels
ax.set_ylim(0, 20)
ax.set_xlim(0, 1.1)
ax.set_xlabel(r'Normalised acoustic wavenumber $q/(2\beta)$')
ax.set_ylabel(r'Frequency $\Omega/(2\pi)$ [GHz]')
fig.savefig(prefix+'-dispersion_multicore.png', bbox_inches='tight')


print(nbapp.final_report())
