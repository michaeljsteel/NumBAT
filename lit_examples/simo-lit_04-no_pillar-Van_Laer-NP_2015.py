""" Replicating the results of
    Interaction between light and highly confined
    hypersound in a silicon photonic nanowire
    Van Laer et al.
    http://dx.doi.org/10.1038/nphoton.2015.11

    Making simplification of ignoring the pedestal.
"""

import sys
import copy
import time
import numpy as np

sys.path.append("../backend/")
import numbat
import materials
import mode_calcs
import integration
import plotting
from plotting import Decorator

import starter


# use this class to add or alter features to the final plots
class EMDecorator(Decorator):
  def __init__(self):
    super().__init__()

  def extra_axes_commands(self, ax):
    ax.tick_params(axis='x',color='gray', which='both')
    ax.tick_params(axis='y',color='gray', which='both')
    if self.is_single_plot():
      ax.tick_params(axis='x',length=20)
      ax.tick_params(axis='y',length=20)
      ax.tick_params(axis='x',width=2)
      ax.tick_params(axis='y',width=2)

      ax.tick_params(axis='x',length=10,which='minor')
      ax.tick_params(axis='y',length=10,which='minor')
      ax.tick_params(axis='x',width=2, which='minor')
      ax.tick_params(axis='y',width=2, which='minor')


emdecorate=EMDecorator()
#acdecorate=ACDecorator()


# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 5*wl_nm
unitcell_y = 0.5*unitcell_x
inc_a_x = 485
inc_a_y = 230
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 60
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(4, sys.argv, sub='a')

nbapp = numbat.NumBATApp(prefix)
# Rotate crystal axis of Si from <100> to <110>, starting with same Si_2016_Smith data.
Si_110 = copy.deepcopy(materials.make_material("Si_2016_Smith"))
# Si_110 = copy.deepcopy(materials.materials_dict["Si_2015_Van_Laer"])
Si_110.rotate_axis(np.pi/4,'z-axis', save_rotated_tensors=True)
# Use all specified parameters to create a waveguide object.
wguide = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=Si_110, symmetry_flag=False,
                        lc_bkg=.1, lc_refine_1=20.0*refine_fac, lc_refine_2=20.0*refine_fac)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1


doem=True
doac=True
new_calcs=True

if doem:
  # Calculate Electromagnetic Modes
  if new_calcs:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
    np.savez(prefix+'-wguide_data', sim_EM_pump=sim_EM_pump)
  else:
    npzfile = np.load(prefix+'-wguide_data.npz', allow_pickle=True)
    sim_EM_pump = npzfile['sim_EM_pump'].tolist()

  sim_EM_Stokes = mode_calcs.fwd_Stokes_modes(sim_EM_pump)
  np.savez(prefix+'-wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
  #npzfile = np.load(prefix+'-wguide_data2.npz', allow_pickle=True)
  #sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

  plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.43, xlim_max=0.43, ivals=[EM_ival_pump],
                           ylim_min=0.43, ylim_max=0.43, EM_AC='EM_E',
                           n_points=2000, quiver_points=10,  decorator=emdecorate)

  # Print the wavevectors of EM modes.
  print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.kz_EM_all()), 4))

  # Calculate the EM effective index of the waveguide.
  print("n_eff = ", np.round(sim_EM_pump.neff_all(), 4))


if doac:
  q_AC = 5 # close but not quite zero

  # Calculate Acoustic Modes
  if new_calcs:
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)
    #np.savez(prefix+'-wguide_data_AC', sim_AC=sim_AC)  # Pickle problem
  else:
    npzfile = np.load(prefix+'-wguide_data_AC.npz', allow_pickle=True)
    sim_AC = npzfile['sim_AC'].tolist()

  plotting.plot_mode_fields(sim_AC, ivals=range(20), )

set_q_factor = 306

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC, EM_ival_pump=EM_ival_pump,
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Print the frequencies and gains of AC modes.
nus = sim_AC.nu_AC_all()
gains = SBS_gain[EM_ival_pump, EM_ival_Stokes, :]

print('Acoustic modes')
print('m      Freq (GHz)      Total gain (1/(Wm))')
for m in range(len(nus)):
    print('{0:3d}    {1:12.6f}   {2:8.4f}'.format(m, np.real(nus[m])*1e-9, gains[m]))


# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)

print("\n Displaying results with negligible components masked out")
print("SBS_gain [1/(Wm)] PE contribution \n", masked_PE)
print("SBS_gain [1/(Wm)] MB contribution \n", masked_MB)
print("SBS_gain [1/(Wm)] total \n", masked)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min =5e9 # 9.1e9 # Hz
freq_max = 50e9 # 9.3e9 # Hz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    )

print(nbapp.final_report())
