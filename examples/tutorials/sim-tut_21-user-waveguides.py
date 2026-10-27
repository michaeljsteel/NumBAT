"""Calculate optical and elastic modes with an externally registered waveguide."""

import sys
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT / "backend"))

import materials
import numbat


print("\nCommencing NumBAT tutorial 21")

# Register a source-tree collection of user waveguides before creating NumBATApp.
numbat.register_waveguide_directory(PROJECT_ROOT / "examples" / "user_waveguides")

prefix = "tut_21"
lambda_nm = 1550

# All geometric dimensions are in nm.
domain_x = 4000
domain_y = 3000
inc_a_x = 1200
inc_a_y = 800

num_modes_em = 20
num_modes_ac = 20

background = materials.make_material("SiO2_2016_Smith")
core = materials.make_material("As2S3_2023_Steel")

nbapp = numbat.NumBATApp(prefix)

wguide = nbapp.make_structure(
    inc_shape="elliptical_rod",
    domain_x=domain_x,
    domain_y=domain_y,
    inc_a_x=inc_a_x,
    inc_a_y=inc_a_y,
    material_bkg=background,
    material_a=core,
    lc_bkg=0.05,
    lc_refine_1=5,
    lc_refine_2=5,
)

wguide.plot_mesh(prefix)
wguide.plot_refractive_index_profile(prefix)

# Solve the optical modes of the externally supplied EllipticalRod geometry.
sim_em = wguide.calc_EM_modes(num_modes_em, lambda_nm, n_eff=1.4)
sim_em.plot_modes(mode_indices=[0])

print("\nOptical effective indices:")
for mode_index, effective_index in enumerate(sim_em.neff_all()):
    print(f"{mode_index:3d}  {np.real(effective_index):.6f}")

# A nonzero acoustic wavevector gives a representative elastic-mode calculation.
# It is the usual backward-scattering value, but this tutorial does not calculate SBS gain.
q_ac = 2 * np.real(sim_em.kz_EM(0))
sim_ac = wguide.calc_AC_modes(num_modes_ac, q_ac, EM_sim=sim_em, shift_Hz=4e9)
sim_ac.plot_modes(mode_indices=range(3))

print("\nElastic mode frequencies (GHz):")
for mode_index, frequency in enumerate(sim_ac.nu_AC_all()):
    print(f"{mode_index:3d}  {np.real(frequency) * 1e-9:.6f}")

print(nbapp.final_report())
