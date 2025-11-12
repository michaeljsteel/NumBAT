from pathlib import Path
from typing import Optional, Tuple

from PIL import Image

import nbgmsh
import numbat
import plotting.plottools as plottools
from numbattools import f2f_with_subs, run_subprocess


def plot_mail_mesh(mesh_mail_fname: Path | str, outpref: Optional[str]) -> None:
    """Visualise the mesh in .mail format.

    Args:
        mesh_mail_fname: Path to the .mail mesh file.
        outpref: Output prefix for the plot files.
    """
    path = Path(numbat.NumBATApp().outpath())
    mail_data = nbgmsh.MailData(mesh_mail_fname)
    mail_data.plot_mesh(path)


def plot_gmesh_wireframe(
    msh_loc_in: Path | str, msh_loc_out: Path | str, msh_name: str, outpref: str
) -> str:
    """Create a wireframe visualization of the geometry.

    Args:
        msh_loc_in: Input directory containing gmsh script templates.
        msh_loc_out: Output directory for generated files.
        msh_name: Base name of the mesh file.
        outpref: Output prefix for the plot files.

    Returns:
        Path to the generated wireframe PNG file.
    """
    nbapp = numbat.NumBATApp()
    gmsh_exe = str(nbapp.path_gmsh())

    tmpoutpref = Path.cwd() / outpref

    # Make the wire frame image
    scr_in = Path(msh_loc_in) / "geo2png.scr"
    scr_out = Path(msh_loc_out) / f"{msh_name}_geo2png.scr"

    fn_out = Path(f"{tmpoutpref}-wireframe")
    f2f_with_subs(scr_in, scr_out, {"tmp": str(fn_out)})
    cmd = [gmsh_exe, f"{msh_name}.geo", scr_out.name]

    run_subprocess(cmd, "Gmsh", cwd=str(msh_loc_out))

    # tidy it up
    fin_out = fn_out.with_suffix(".png")
    im = Image.open(fin_out)

    im = plottools.fig_trim_border(im, clip=(0, 0, 1, 0), trimwhite=True, border=20)
    im.save(fin_out)

    return str(fin_out)


def plot_gmesh_mesh(
    msh_loc_in: Path | str, msh_loc_out: Path | str, msh_name: str, outpref: str
) -> str:
    """Create a visualization of the mesh nodes.

    Args:
        msh_loc_in: Input directory containing gmsh script templates.
        msh_loc_out: Output directory for generated files.
        msh_name: Base name of the mesh file.
        outpref: Output prefix for the plot files.

    Returns:
        Path to the generated mesh PNG file.
    """
    nbapp = numbat.NumBATApp()
    gmsh_exe = str(nbapp.path_gmsh())

    tmpoutpref = Path.cwd() / outpref

    # Make the mesh image
    scr_in = Path(msh_loc_in) / "msh2png.scr"
    scr_out = Path(msh_loc_out) / f"{msh_name}_msh2png.scr"
    fn_out = Path(f"{tmpoutpref}-mesh_nodes")
    f2f_with_subs(scr_in, scr_out, {"tmp": str(fn_out)})

    cmd = [gmsh_exe, f"{msh_name}.msh", scr_out.name]
    run_subprocess(cmd, "Gmsh", cwd=str(msh_loc_out))

    # tidy
    fin_out = fn_out.with_suffix(".png")
    im = Image.open(fin_out)
    im = plottools.fig_trim_border(im, clip=(0, 0, 1, 0), trimwhite=True, border=20)
    im.save(fin_out)

    return str(fin_out)


def plot_mesh(
    msh_loc_in: Path | str,
    msh_loc_out: Path | str,
    msh_name: str,
    outpref: str,
    combo_plot: bool = True,
) -> Tuple[str, str, str]:
    """Visualise mesh with gmsh and save to a file.

    Args:
        msh_loc_in: Input directory containing gmsh script templates.
        msh_loc_out: Output directory for generated files.
        msh_name: Base name of the mesh file.
        outpref: Output prefix for the plot files.
        combo_plot: If True, create a combined plot of wireframe and mesh.

    Returns:
        Tuple of (wireframe_file, mesh_file, combined_file) paths.
    """
    # Manipulate scripts in backend/fortran/build
    # Writes final png file to user directory

    # Make the individual ones
    fout_wire = plot_gmesh_wireframe(msh_loc_in, msh_loc_out, msh_name, outpref)
    fout_mesh = plot_gmesh_mesh(msh_loc_in, msh_loc_out, msh_name, outpref)

    # join them
    fout_join = ""
    if combo_plot:
        nbapp = numbat.NumBATApp()

        outprefix = Path(nbapp.outpath(outpref))

        fout_join = str(outprefix.with_name(f"{outprefix.name}-mesh.png"))

        plottools.join_figs(
            [
                fout_wire,
                fout_mesh,
            ],
            fout_join,
            # clip=(0,0,1,0), # gmsh leaves a vertical black line on right edge
            # trimwhite=True
            delete_inputs=True,
        )
    return fout_wire, fout_mesh, fout_join


def check_mesh(msh_loc_out: Path | str, msh_name: str) -> None:
    """Visualise geometry and mesh with gmsh.

    Opens the geometry and mesh files in the gmsh GUI for interactive inspection.

    Args:
        msh_loc_out: Directory containing the mesh files.
        msh_name: Base name of the mesh file (without extension).
    """
    nbapp = numbat.NumBATApp()
    gmsh_exe = str(nbapp.path_gmsh())

    msh_loc_out = Path(msh_loc_out)

    gmsh_cmd = [gmsh_exe, str(msh_loc_out / f"{msh_name}.geo")]
    run_subprocess(gmsh_cmd, "Gmsh", cwd=str(msh_loc_out))

    gmsh_cmd = [gmsh_exe, str(msh_loc_out / f"{msh_name}.msh")]
    run_subprocess(gmsh_cmd, "Gmsh", cwd=str(msh_loc_out))
    run_subprocess(gmsh_cmd, "Gmsh", cwd=str(msh_loc_out))
