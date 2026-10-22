"""Calculate EM and acoustic modes of a rounded silicon bowtie waveguide."""

import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[2] / "backend"))

import materials  # noqa: I001
import numbat

def output_geometry(wguide, prefix):
    #wguide.check_mesh()
    wguide.plot_mesh(prefix)

    index_plotter = wguide.get_structure_plotter_refractive_index()
    index_plotter.make_plot_2D(prefix)
    index_plotter.make_plot_1D_cut(prefix, "y", 0)

    # Needs fixing
    #c11_plotter = wguide.get_structure_plotter_stiffness(1,1)
    #c11_plotter.make_plot_2D(prefix)
    #c11_plotter.make_plot_1D_cut(prefix, "y", 0)

    vp_plotter = wguide.get_structure_plotter_acoustic_velocity()
    vp_plotter.make_plot_2D(prefix)
    vp_plotter.make_plot_1D_cut(prefix, "y", 0)

def bowtie_field_intensity_ratio(sim_EM, bowlen_l, mode_index, prefix):
    """Return |E(bowlen_l, 0)|^2 / |E(0, 0)|^2 from a saved horizontal field slice."""
    sim_EM.write_modes_1D("y", 0, mode_indices=[mode_index], n_points=5001, prefix=prefix)
    field_path = numbat.NumBATApp().outdir_fields_path(prefix) / f"EM_mode_{mode_index:02d}_ycut.txt"
    field_data = np.loadtxt(field_path)
    x_um = field_data[:, 0]
    ref_x_um = -bowlen_l * 0.001  # convert nm to um

    if not x_um[0] <= ref_x_um <= x_um[-1]:
        raise ValueError(f"bowlen_l={bowlen_l} nm is outside the saved horizontal field slice.")

    electric_field = field_data[:, 1:7:2] + 1j * field_data[:, 2:7:2]

    def field_intensity_at(target_x_um):
        field_at_target = np.array(
            [
                np.interp(target_x_um, x_um, electric_field[:, component].real)
                + 1j * np.interp(target_x_um, x_um, electric_field[:, component].imag)
                for component in range(electric_field.shape[1])
            ]
        )
        return np.sum(np.abs(field_at_target) ** 2)

    intensity_at_origin = field_intensity_at(0)
    del_x = .01 # come in a little since the field may not be smooth at the bowtie edge
    target_x_um = ref_x_um + del_x
    intensity_at_edge = field_intensity_at(target_x_um)

    if intensity_at_edge == 0:
        raise ZeroDivisionError("The electric-field intensity at the bowtie edge is zero.")

    return  intensity_at_origin / intensity_at_edge




# Geometric parameters in nm. alpha_tip is a half-angle in radians.
lambda_nm = 1550
domain_x = 3000
domain_y = 3000
rad_tip = 20
alpha_tip = np.deg2rad(20)
gap_g = 30
bowlen_l = 220
rad_corner = 20
rad_wg = 500

num_modes_EM = 50
num_modes_AC = 100
em_mode_index = 0
ac_mode_index = 0
q_AC = 4e7  # 1/m

prefix = "tut_20_rounded_bowtie"
refine_fac = 1
nbapp = numbat.NumBATApp(prefix)

vacuum = materials.make_material("Vacuum")
silica = materials.make_material("SiO2_2023_Steel")
silicon = materials.make_material("Si_2016_Smith")
chalc = materials.make_material("As2S3_2023_Steel")

one_wedge=False
mat_bg = silica
mat_a = chalc
mat_b = mat_bg if one_wedge else mat_a
wguide = nbapp.make_structure(
    "rounded_bowtie",
    domain_x,
    domain_y,
    rad_tip=rad_tip,
    rad_wg=rad_wg,
    alpha_tip=alpha_tip,
    gap_g=gap_g,
    bowlen_l=bowlen_l,
    rad_corner=rad_corner,
    material_bkg=mat_bg,
    material_a=mat_a,
    material_b=mat_b,
    lc_bkg=0.02,
    lc_refine_1=5.0 * refine_fac,
    lc_refine_2=5.0 * refine_fac,
)

# Plot the 2D refractive-index profile and its horizontal cut through y = 0 um.
output_geometry(wguide, prefix)

# Solve and plot the fundamental EM mode in 2D and along the horizontal centre line.
emout_m=range(5)
acout_m=range(50)
n_eff = silicon.refindex_n - 0.1
sim_EM = wguide.calc_EM_modes(num_modes_EM, lambda_nm, n_eff=n_eff)
sim_EM.plot_modes(mode_indices=emout_m, hide_vector_field=True)
sim_EM.plot_modes(mode_indices=emout_m, hide_vector_field=True,
                  xlim_min=.3, ylim_min=.3, xlim_max=.3, ylim_max=.3, suffix='zoom')
sim_EM.plot_modes_1D("x", 0, mode_indices=emout_m)
sim_EM.plot_modes_1D("y", 0, mode_indices=emout_m)



field_intensity_ratio = bowtie_field_intensity_ratio(sim_EM, bowlen_l, em_mode_index, prefix)

print(
    "\nNormalized field intensity at the outer bowtie boundary: "
    f"(||E(0, 0)/E({bowlen_l} nm, 0)| )^2 = {field_intensity_ratio:.6e}"
)


sys.exit(0)
print("\nEM effective indices:")
for mode_index, mode_neff in enumerate(sim_EM.neff_all()):
    print(f"{mode_index:3d}  {mode_neff:.6f}")

# Solve and plot the requested acoustic-wavevector modes.
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM)
sim_AC.plot_modes(mode_indices=acout_m, quiver_points=20)
#sim_AC.plot_modes_1D("y", 0, mode_indices=acout_m)

print("\nAcoustic mode frequencies (GHz):")
for mode_index, frequency in enumerate(sim_AC.nu_AC_all()):
    print(f"{mode_index:3d}  {np.real(frequency) * 1e-9:.6f}")

print(nbapp.final_report())
