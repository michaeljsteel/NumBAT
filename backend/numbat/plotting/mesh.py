from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from numbat_core import runtime
from numbat_core.mesh.gmsh_io import read_msh2_mesh
from numbat_core.mesh.mail import MailData
from numbat_core.plotting import tools as plottools
from numbat_core.plotting.axes import TidyAxes
from numbat_core.plotting.tools import save_and_close_figure
from numbat_core.utils.numeric import f2f_with_subs, run_subprocess


def plot_mail_mesh(mesh_mail_fname: Path | str, outpref: Optional[str]) -> None:
    """Visualise the mesh in .mail format.

    Args:
        mesh_mail_fname: Path to the .mail mesh file.
        outpref: Output prefix for the plot files.
    """
    path = Path(runtime.outpath())
    plot_mail_data(MailData(mesh_mail_fname), path)


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
    gmsh_exe = str(runtime.path_gmsh())

    tmpoutpref = Path.cwd() / outpref

    # Make the wire frame image
    scr_in = Path(msh_loc_in) / "geo2png.scr"
    scr_out = Path(msh_loc_out) / f"{msh_name}_geo2png.scr"

    fn_out = Path(f"{tmpoutpref}-wireframe")
    f2f_with_subs(scr_in, scr_out, {"tmp": str(fn_out)})
    cmd = [gmsh_exe, f"{msh_name}.geo", scr_out.name]

    cmd_flat = ' '.join(cmd)
    print(f'Calling gmsh cmd: "{cmd_flat}"')
    exit_code = run_subprocess(cmd, "Gmsh", cwd=str(msh_loc_out))
    if exit_code:
        print(f'Gmsh plot command "{cmd_flat}" returned error code {exit_code}.')

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
    gmsh_exe = str(runtime.path_gmsh())

    tmpoutpref = Path.cwd() / outpref

    # Make the mesh image
    scr_in = Path(msh_loc_in) / "msh2png.scr"
    scr_out = Path(msh_loc_out) / f"{msh_name}_msh2png.scr"
    fn_out = Path(f"{tmpoutpref}-mesh_nodes")
    f2f_with_subs(scr_in, scr_out, {"tmp": str(fn_out)})

    cmd = [gmsh_exe, f"{msh_name}.msh", scr_out.name]
    cmd_flat = ' '.join(cmd)
    print(f'Calling gmsh cmd: "{cmd_flat}" in cwd="{str(msh_loc_out)}"')

    exit_code = run_subprocess(cmd, "Gmsh", cwd=str(msh_loc_out))

    if exit_code:
        print(f'Gmsh plot command "{cmd_flat}" returned error code {exit_code}.')

    # tidy
    fin_out = fn_out.with_suffix(".png")
    im = Image.open(fin_out)
    im = plottools.fig_trim_border(im, clip=(0, 0, 1, 0), trimwhite=True, border=20)
    im.save(fin_out)

    return str(fin_out)


def _render_gmsh_png(
    msh_loc_in: Path | str,
    msh_loc_out: Path | str,
    msh_name: str,
    outpref: str,
    template_name: str,
    input_suffix: str,
    output_suffix: str,
) -> Path:
    """Render a Gmsh model through its GUI context and return the PNG path."""
    gmsh_exe = str(runtime.path_gmsh())
    output_prefix = Path.cwd() / outpref
    output_stem = output_prefix.with_name(f"{output_prefix.name}-{output_suffix}")
    output_png = output_stem.with_suffix(".png")

    template_path = Path(msh_loc_in) / template_name
    render_script = Path(msh_loc_out) / f"{msh_name}_{template_name.removesuffix('.scr')}.geo"
    f2f_with_subs(template_path, render_script, {"tmp": str(output_stem)})

    cmd = [gmsh_exe, "-nopopup", f"{msh_name}{input_suffix}", render_script.name]
    cmd_flat = " ".join(cmd)
    print(f'Calling Gmsh render command: "{cmd_flat}" in cwd="{msh_loc_out}"')
    exit_code = run_subprocess(cmd, "Gmsh", cwd=str(msh_loc_out), exit_on_fail=False)
    if exit_code:
        raise RuntimeError(f'Gmsh render command failed with exit code {exit_code}: "{cmd_flat}"')
    if not output_png.is_file():
        raise RuntimeError(
            f'Gmsh render command completed without creating "{output_png}": "{cmd_flat}"'
        )

    return output_png


def plot_gmsh_wireframe(
    msh_loc_in: Path | str, msh_loc_out: Path | str, msh_name: str, outpref: str
) -> str:
    """Create a wireframe PNG using Gmsh's GUI renderer without popup dialogs."""
    output_png = _render_gmsh_png(
        msh_loc_in,
        msh_loc_out,
        msh_name,
        outpref,
        "geo2png.scr",
        ".geo",
        "wireframe",
    )
    image = Image.open(output_png)
    image = plottools.fig_trim_border(image, clip=(0, 0, 1, 0), trimwhite=True, border=20)
    image.save(output_png)
    return str(output_png)


def plot_gmsh_mesh(
    msh_loc_in: Path | str, msh_loc_out: Path | str, msh_name: str, outpref: str
) -> str:
    """Create a mesh-node PNG using Gmsh's GUI renderer without popup dialogs."""
    output_png = _render_gmsh_png(
        msh_loc_in,
        msh_loc_out,
        msh_name,
        outpref,
        "msh2png.scr",
        ".msh",
        "mesh_nodes",
    )
    image = Image.open(output_png)
    image = plottools.fig_trim_border(image, clip=(0, 0, 1, 0), trimwhite=True, border=20)
    image.save(output_png)
    return str(output_png)


def _wireframe_edges(
    triangles: list[tuple[tuple[int, int, int], int]], line_elements: set[tuple[int, int]]
) -> tuple[set[tuple[int, int]], set[tuple[int, int]]]:
    """Return exterior and interior edges that define the physical geometry."""
    edge_materials = {}
    edge_counts = {}
    for node_ids, physical_id in triangles:
        for edge in (
            tuple(sorted((node_ids[0], node_ids[1]))),
            tuple(sorted((node_ids[1], node_ids[2]))),
            tuple(sorted((node_ids[2], node_ids[0]))),
        ):
            edge_materials.setdefault(edge, set()).add(physical_id)
            edge_counts[edge] = edge_counts.get(edge, 0) + 1

    exterior_edges = {
        edge
        for edge, physical_ids in edge_materials.items()
        if edge_counts[edge] == 1
    }
    interior_edges = {
        edge for edge, physical_ids in edge_materials.items() if len(physical_ids) > 1
    }

    for edge in line_elements:
        (exterior_edges if edge in exterior_edges else interior_edges).add(edge)

    return exterior_edges, interior_edges


def _draw_edges(
    axes: plt.Axes,
    coordinates: np.ndarray,
    edges: set[tuple[int, int]],
    color: str,
    linewidth: float,
    alpha: float = 1.0,
) -> None:
    """Draw a collection of mesh edges."""
    for start_index, end_index in edges:
        axes.plot(
            coordinates[[start_index, end_index], 0],
            coordinates[[start_index, end_index], 1],
            color=color,
            linewidth=linewidth,
            alpha=alpha,
        )


def _save_msh_plot(
    output_png: Path,
    coordinates: np.ndarray,
    exterior_edges: set[tuple[int, int]],
    interior_edges: set[tuple[int, int]],
    triangles: list[tuple[tuple[int, int, int], int]] | None = None,
) -> None:
    """Save a Matplotlib rendering of MSH 2.2 boundary edges or triangular elements."""
    figure, axes = plt.subplots()
    if triangles is None:
        _draw_edges(axes, coordinates, exterior_edges, "#2e7d32", 1.2)
        _draw_edges(axes, coordinates, interior_edges, "#1565c0", 0.8)
    else:
        triangle_nodes = np.array([node_ids for node_ids, _ in triangles])
        axes.triplot(
            coordinates[:, 0],
            coordinates[:, 1],
            triangle_nodes,
            color="#5b9bd5",
            linewidth=0.25,
        )
        _draw_edges(
            axes,
            coordinates,
            exterior_edges | interior_edges,
            "#a5543f",
            0.7,
            alpha=0.65,
        )

    axes.set_aspect("equal")
    axes.set_xlabel(r"$x$ [μm]")
    axes.set_ylabel(r"$y$ [μm]")
    axes.tick_params(colors="#666666")
    axes.xaxis.label.set_color("#666666")
    axes.yaxis.label.set_color("#666666")
    for spine in axes.spines.values():
        spine.set_color("#666666")
    figure.savefig(output_png, dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(figure)


def plot_msh_wireframe(
    msh_loc_out: Path | str, msh_name: str, outpref: str
) -> str:
    """Create a deterministic wireframe PNG from a Gmsh 2.2 mesh file."""
    output_png = (Path.cwd() / outpref).with_name(f"{outpref}-wireframe.png")
    coordinates, triangles, line_elements = read_msh2_mesh(
        Path(msh_loc_out) / f"{msh_name}.msh"
    )
    exterior_edges, interior_edges = _wireframe_edges(triangles, line_elements)
    _save_msh_plot(
        output_png,
        coordinates,
        exterior_edges,
        interior_edges,
    )
    return str(output_png)


def plot_msh_mesh(msh_loc_out: Path | str, msh_name: str, outpref: str) -> str:
    """Create a deterministic triangular-mesh PNG from a Gmsh 2.2 mesh file."""
    output_png = (Path.cwd() / outpref).with_name(f"{outpref}-mesh_nodes.png")
    coordinates, triangles, line_elements = read_msh2_mesh(
        Path(msh_loc_out) / f"{msh_name}.msh"
    )
    exterior_edges, interior_edges = _wireframe_edges(triangles, line_elements)
    _save_msh_plot(output_png, coordinates, exterior_edges, interior_edges, triangles)
    return str(output_png)


def plot_mesh(
    msh_loc_in: Path | str,
    msh_loc_out: Path | str,
    msh_name: str,
    outpref: str,
    combo_plot: bool = True,
    keep_singles: bool = False,
) -> Tuple[str, str, str]:
    """Visualise mesh with gmsh and save to a file.

    Args:
        msh_loc_in: Input directory containing gmsh script templates.
        msh_loc_out: Output directory for generated files.
        msh_name: Base name of the mesh file.
        outpref: Output prefix for the plot files.
        combo_plot: If True, create a combined plot of wireframe and mesh.
        keep_singles: If True, keep the single plots when a combo is made.

    Returns:
        Tuple of (wireframe_file, mesh_file, combined_file) paths.
    """
    # Manipulate scripts in backend/fortran/build
    # Writes final png file to user directory

    # Render from the mesh in Python to avoid FLTK/OpenGL-dependent Gmsh screenshots.
    fout_wire = plot_msh_wireframe(msh_loc_out, msh_name, outpref)
    fout_mesh = plot_msh_mesh(msh_loc_out, msh_name, outpref)

    # join them
    fout_join = ""
    if combo_plot:
        outprefix = Path(runtime.outpath(outpref))

        fout_join = str(outprefix.with_name(f"{outprefix.name}-mesh.png"))

        plottools.join_figs(
            [
                fout_wire,
                fout_mesh,
            ],
            fout_join,
            # clip=(0,0,1,0), # gmsh leaves a vertical black line on right edge
            # trimwhite=True
            delete_inputs= not keep_singles,
        )
    if combo_plot:
        if keep_singles:
            return fout_join, fout_wire, fout_mesh
        else:
            return fout_join
    else:
        return fout_wire, fout_mesh


def check_mesh(msh_loc_out: Path | str, msh_name: str) -> None:
    """Visualise geometry and mesh with gmsh.

    Opens the geometry and mesh files in the gmsh GUI for interactive inspection.

    Args:
        msh_loc_out: Directory containing the mesh files.
        msh_name: Base name of the mesh file (without extension).
    """
    gmsh_exe = str(runtime.path_gmsh())

    msh_loc_out = Path(msh_loc_out)

    gmsh_cmd = [gmsh_exe, str(msh_loc_out / f"{msh_name}.geo")]
    run_subprocess(gmsh_cmd, "Gmsh", cwd=str(msh_loc_out))

    gmsh_cmd = [gmsh_exe, str(msh_loc_out / f"{msh_name}.msh")]
    run_subprocess(gmsh_cmd, "Gmsh", cwd=str(msh_loc_out))
    run_subprocess(gmsh_cmd, "Gmsh", cwd=str(msh_loc_out))
def plot_mail_data(mail_data, path):
    fig, axs = plt.subplots(2,2, dpi=300)
    fig.subplots_adjust(hspace=0.3, wspace=.1)


    axs = axs.flat
    cf=axs[0].scatter(mail_data.v_x, mail_data.v_y, s=4, c=mail_data.v_node_physindex, cmap='Dark2')

    axs[0].set_title('Gmsh surface and curve types', fontsize=8)

    cbar = fig.colorbar(cf, shrink=.5, pad=.025, location='right')
    cbar.set_label(label='Phys curve/surface number')


    elts = mail_data.v_elts-1  # built for fortran indexing
    axs[1].scatter(mail_data.v_x[elts[:,0]], mail_data.v_y[elts[:,0]], s=6, c='red')
    axs[1].scatter(mail_data.v_x[elts[:,1]], mail_data.v_y[elts[:,1]], s=6, c='red')
    axs[1].scatter(mail_data.v_x[elts[:,2]], mail_data.v_y[elts[:,2]], s=6, c='red')
    axs[1].scatter(mail_data.v_x[elts[:,3]], mail_data.v_y[elts[:,3]], s=4, c='blue')
    axs[1].scatter(mail_data.v_x[elts[:,4]], mail_data.v_y[elts[:,4]], s=4, c='blue')
    axs[1].scatter(mail_data.v_x[elts[:,5]], mail_data.v_y[elts[:,5]], s=4, c='blue')

    axs[1].set_title('Vertex and edge nodes', fontsize=8)

    # axs[2].scatter(mail_data.v_x[elts[:,1]], mail_data.v_y[elts[:,1]], s=4, c='red')
    # axs[2].scatter(mail_data.v_x[elts[:,2]], mail_data.v_y[elts[:,2]], s=4, c='green')
    # axs[2].scatter(mail_data.v_x[elts[:,3]], mail_data.v_y[elts[:,3]], s=4, c='blue')
    # axs[3].scatter(mail_data.v_x[elts[:,4]], mail_data.v_y[elts[:,4]], s=4, c='cyan')
    # axs[3].scatter(mail_data.v_x[elts[:,5]], mail_data.v_y[elts[:,5]], s=4, c='magenta')
    # axs[3].scatter(mail_data.v_x[elts[:,6]], mail_data.v_y[elts[:,6]], s=4, c='yellow')

    # Draw elements (solid) and sub-element (dotted) triangles
    cmap = plt.get_cmap('tab10')
    for i in range(mail_data.n_msh_elts):
        icor = [elts[i,0],elts[i,1],elts[i,2],elts[i,0]]
        imid = [elts[i,3],elts[i,4],elts[i,5],elts[i,3]]
        ic = (2*i)%10
        axs[2].plot(mail_data.v_x[icor], mail_data.v_y[icor], color=cmap(ic), lw=1)
        axs[2].plot(mail_data.v_x[imid], mail_data.v_y[imid], ':',color=cmap(ic+1), lw=.5)

    for off in range(mail_data.n_msh_elts):  # look for top left element
        if 10 in elts[off,:6]:
            break
    for i in range(6):
        axs[2].plot(mail_data.v_x[elts[off, i]],mail_data.v_y[elts[off, i]], 'o',color=cmap(i), markersize=3, label=f'Node {i}')
    axs[2].legend(loc='upper right', fontsize=8)
    axs[2].set_xlim(min(mail_data.v_x), max(mail_data.v_x))
    axs[2].set_ylim(min(mail_data.v_y), max(mail_data.v_y))
    axs[2].set_title('NumBAT FEM sub-elements', fontsize=8)

    for ax in axs:
        ax.set_xlabel(r'$x$ [μm]')
        ax.set_ylabel(r'$y$ [μm]')

    tidyaxes = TidyAxes(nax=4, props={'ax_label_fs':8, 'ax_ticklabel_fs':6,
                                 'ax_tickwidth':0.25, 'ax_linewidth':1, 'aspect':1,
                                 'cb_linewidth':1, 'cb_label_fs':10,
                                 'cb_ticklabel_fs':8, 'cb_tickwidth':0.25})


    tidyaxes.apply_to_axes(axs)
    tidyaxes.apply_to_cbars(cbar)

    save_and_close_figure(fig, str(path) + '-mailmesh.png')
