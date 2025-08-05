from pathlib import Path
from PIL import Image

import numbat
import nbgmsh

import plotting.plottools as plottools
from numbattools import f2f_with_subs, run_subprocess

def plot_mail_mesh(mesh_mail_fname, outpref):
    """Visualise the mesh in .mail format."""

    path = numbat.NumBATApp().outpath()
    mail_data = nbgmsh.MailData(mesh_mail_fname)
    mail_data.plot_mesh(path)


def plot_gmesh_wireframe(msh_loc_in, msh_loc_out, msh_name, outpref):
    nbapp = numbat.NumBATApp()
    gmsh_exe = nbapp.path_gmsh()

    #tdir = tempfile.TemporaryDirectory()
    #tmpoutpref = str(Path(tdir.name, outpref))
    tmpoutpref = Path.cwd() / outpref

    # Make the wire frame image
    scr_in = Path(msh_loc_in) / 'geo2png.scr'
    scr_out = Path(msh_loc_out) / (msh_name + '_geo2png.scr')

    fn_out = str(tmpoutpref) + '-wireframe'
    f2f_with_subs(scr_in, scr_out, {'tmp': fn_out})
    cmd = [gmsh_exe, msh_name + '.geo', scr_out.name]

    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    # tidy it up
    fin_out = fn_out + '.png'
    im = Image.open(fin_out)

    im = plottools.fig_trim_border(im, clip=(0,0,1,0), trimwhite=True, border = 20)
    im.save(fin_out)

    return fin_out

def plot_gmesh_mesh(msh_loc_in, msh_loc_out, msh_name, outpref):

    nbapp = numbat.NumBATApp()
    gmsh_exe = nbapp.path_gmsh()

    #tdir = tempfile.TemporaryDirectory()
    #tmpoutpref = str(Path(tdir.name, outpref))
    #tmpoutpref = nbapp.outpath(outpref)
    tmpoutpref = Path.cwd() / outpref


    # Make the mesh image
    scr_in = Path(msh_loc_in) / 'msh2png.scr'
    scr_out = Path(msh_loc_out) / (msh_name + '_msh2png.scr')
    fn_out = str(tmpoutpref) + '-mesh_nodes'
    f2f_with_subs(scr_in, scr_out, {'tmp': fn_out})

    cmd = [gmsh_exe, msh_name + '.msh', scr_out.name]

    if numbat.NumBATApp().is_macos(): # macOS requires the -a flag to open applications
        cmd = ['open', '-a'] + cmd


    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    # tidy
    fin_out = fn_out + '.png'
    im = Image.open(fin_out)
    im = plottools.fig_trim_border(im, clip=(0,0,1,0), trimwhite=True, border = 20)
    im.save(fin_out)

    return fin_out


def plot_mesh(msh_loc_in, msh_loc_out, msh_name, outpref,
              combo_plot=True):
    """Visualise mesh with gmsh and save to a file."""

    # Manipulate scripts in backend/fortran/build
    # Writes final png file to user directory

    # Make the individual ones
    fout_wire = plot_gmesh_wireframe(msh_loc_in, msh_loc_out, msh_name, outpref)
    fout_mesh = plot_gmesh_mesh(msh_loc_in, msh_loc_out, msh_name, outpref)

    # join them
    fout_join = ''
    if combo_plot:
        nbapp = numbat.NumBATApp()

        outprefix = nbapp.outpath(outpref)

        fout_join = str(outprefix)+'-mesh.png'

        plottools.join_figs([fout_wire, fout_mesh,],
                            fout_join,
                            #clip=(0,0,1,0), # gmsh leaves a vertical black line on right edge
                            #trimwhite=True
                            )
    return fout_wire, fout_mesh, fout_join


def check_mesh(msh_loc_out, msh_name):
    """Visualise geometry and mesh with gmsh."""

    nbapp = numbat.NumBATApp()
    gmsh_exe = str(nbapp.path_gmsh())

    gmsh_cmd1 = [gmsh_exe, f'{msh_loc_out}/{msh_name}.geo']

    gmsh_cmd2 = [gmsh_exe, f'{msh_loc_out}/{msh_name}.msh']

    if numbat.NumBATApp().is_macos(): # macOS requires the -a flag to open applications
        gmsh_cmd1 = ['open', '-a'] + gmsh_cmd1
        gmsh_cmd2 = ['open', '-a'] + gmsh_cmd2

    run_subprocess(gmsh_cmd1, 'Gmsh', cwd=msh_loc_out)
    run_subprocess(gmsh_cmd2, 'Gmsh', cwd=msh_loc_out)
