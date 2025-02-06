""" Replicating the results of
    Brillouin light scattering from surface acoustic
    waves in a subwavelength-diameter optical fibre
    Beugnot et al.
    http://dx.doi.org/10.1038/ncomms6242
"""


import time

import sys
from multiprocessing import Pool

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
sys.path.append(str(Path('../backend')))
import numbat
import materials
import modecalcs
import integration

import starter

def make_gain_plot(m_gain, extents, pref, swap_axes):
    L=0.08  # 8 cm nanowire
    Pp=0.001
    m_gaindB = 10*np.log10(np.abs(m_gain)*L*Pp + 1e-300)


    fig, ax1 = plt.subplots()
    if swap_axes:
        exts = [extents[2], extents[3], extents[0], extents[1]]
        im = ax1.imshow(m_gaindB, aspect='auto', interpolation='none',
                extent=exts, cmap='jet', origin='lower')
        ax1.set_xlabel('Microwire diameter (μm)')
        ax1.set_ylabel('Acoustic frequency (GHz)')
    else:
        im = ax1.imshow(m_gaindB.T, aspect='auto', interpolation='none',
                extent=extents, cmap='jet', origin='lower')
        ax1.set_ylabel('Microwire diameter (μm)')
        ax1.set_xlabel('Acoustic frequency (GHz)')

    for axis in ['top', 'bottom','left','right']:
        ax1.spines[axis].set_linewidth(1)

    cb=fig.colorbar(im, ax=ax1)
    cb.outline.set_visible(False)

    fig.savefig(pref+'-diam_scan.png')



# Select the number of CPUs to use in simulation.
num_cores = 1

# Geometric Parameters - all in nm.
wl_nm = 1550
inc_shape = 'circular'

num_modes_EM_pump = 40
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 200
EM_mode_index_pump = 1    # exploit degeneracy to escape weird mode m=0 fails
EM_mode_index_Stokes = 1  # exploit degeneracy to escape weird mode m=0 fails
AC_mode_index = 'All'

# Expected effective index of fundamental guided mode.
n_eff = 1.44 # just a little under the core index

prefix, refine_fac = starter.read_args(3, sys.argv)

nbapp = numbat.NumBATApp(prefix, prefix+'-out')

prefix = 'lit_03_m4'

freq_min = 5e9
freq_max = 11e9

widerange=False
if 'widerange' in sys.argv[1:]:
    diam_0  =  1050
    diam_min = 600
    diam_max = 3500
    #diam_steps = 301
    diam_steps = 51
    widerange = True

    postfix = '-dwide'
else:
    diam_0  =  1050
    diam_min = 950
    diam_max = 1300
    diam_steps = 101
    postfix = ''


# find a d_diam that will give v_diams hitting diam_0 exactly, with about 20 diam_steps total across (diam_max-diam_min)
d_diam = (diam_0-diam_min)/round((diam_0-diam_min)/((diam_max-diam_min)/diam_steps))
v_diams = np.arange(diam_min, diam_max+d_diam/2, d_diam).round().astype(int) # make sure that diam_0 nm is included.
num_diams = len(v_diams)
print('Using', num_diams, 'diameters:', v_diams)
num_interp_pts = 2000


prefix+=postfix

def modes_n_gain(diam):
    print('Handling diam', diam)

    mat_bkg=materials.make_material("Vacuum")
    mat_a=materials.make_material("SiO2_2023_Steel")
    print(mat_a.elastic_properties())

    inc_a_x = diam
    inc_a_y = diam
    # Use all specified parameters to create a waveguide object.
    domain_x = 3*diam
    domain_y = domain_x
    wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                            material_bkg=mat_bkg, material_a=mat_a,
                            #lc_bkg=.1, lc_refine_1=4.0, lc_refine_2=4.0)
                            lc_bkg=.2/refine_fac, lc_refine_1=refine_fac, lc_refine_2=refine_fac)
    #wguide.plot_mesh(prefix+'_%3d'%int(inc_a_x))
    #wguide.check_mesh()

    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)

    #sim_EM_pump.plot_modes(mode_indices=range(5), field_type='EM_E', )

    sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()

    # find correct fundamental EM mode in case there are spurious ones found as ground state
    EM_mode_index_pump = 0
    for m_p in range(num_modes_EM_pump):
        if sim_EM_pump.neff(m_p)<mat_a.refindex_n: # first mode in physical region
            EM_mode_index_pump = m_p
            break

    EM_mode_index_Stokes = EM_mode_index_pump

    q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) - sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))
    shift_Hz = 4e9
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)



#    set_q_factor = 600.
    gain_box = integration.get_gains_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)#, fixed_Q=set_q_factor)

    #gain_box.set_allowed_EM_pumps(EM_mode_index_pump)
    #gain_box.set_allowed_EM_Stokes(EM_mode_index_Stokes)
    #gain_box.set_EM_modes(EM_mode_index_pump, EM_mode_index_Stokes)

    (nu_gain_tot, nu_gain_PE, nu_gain_MB) = gain_box.plot_spectra(
            freq_min, freq_max, num_interp_pts=num_interp_pts, logy=True,
            prefix = prefix, suffix=f'_w{int(inc_a_x):04d}')

    if abs(diam - diam_0)<1e-10:  # print fields for 1 micron guide
        #pass

        sim_EM_pump.plot_modes(mode_indices=range(10),
                                  prefix=prefix+f'-diam-{round(diam):04d}')
        sim_AC.plot_modes(mode_indices=range(10), prefix=prefix+f'-diam-{round(diam):04d}')
        for m in range(num_modes_AC):
            print(f'{m}, {sim_AC.nu_AC(m)*1e-9:.4f}, {gain_box.gain_total(m):.3e}, ',
                 f'{gain_box.gain_PE(m):.4e}, {gain_box.gain_MB(m):.4e}')


    return (nu_gain_tot, nu_gain_PE, nu_gain_MB)

# Run diams in parallel across num_cores CPUs using multiprocessing package.
if num_cores>1:
    with Pool(num_cores) as pool:
        gain_specs = pool.map(modes_n_gain, v_diams)
else:
    gain_specs = map(modes_n_gain, v_diams)


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
make_gain_plot(m_gain_tot, extents, prefix+'-gain_tot', widerange)
make_gain_plot(m_gain_PE, extents, prefix+'-gain_PE', widerange)
make_gain_plot(m_gain_MB, extents, prefix+'-gain_MB', widerange)

end = time.time()

print(nbapp.final_report())
