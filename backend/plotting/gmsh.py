import numbat
import nbgmsh

import plottools
import tempfile
from pathlib import Path
from numbattools import f2f_with_subs, run_subprocess

def plot_mail_mesh(mesh_mail_fname, outpref):
    """Visualise the mesh in .mail format."""

    path = numbat.NumBATApp().outpath()
    mail_data = nbgmsh.MailData(mesh_mail_fname)
    mail_data.plot_mesh(path)


def plot_gmesh_wireframe(msh_loc_in, msh_loc_out, msh_name, outpref):
    nbapp = numbat.NumBATApp()
    gmsh_exe = nbapp.path_gmsh()

    tdir = tempfile.TemporaryDirectory()
    tmpoutpref = str(Path(tdir.name, outpref))

    # Make the wire frame image
    scr_in = Path(msh_loc_in) / 'geo2png.scr'
    scr_out = Path(msh_loc_out) / (msh_name + '_geo2png.scr')

    fn_out = str(tmpoutpref) + '-entities'
    f2f_with_subs(scr_in, scr_out, {'tmp': fn_out})
    cmd = [gmsh_exe, msh_name + '.geo', scr_out.name]
    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    return fn_out

def plot_gmesh_mesh(msh_loc_in, msh_loc_out, msh_name, outpref):

    nbapp = numbat.NumBATApp()
    gmsh_exe = nbapp.path_gmsh()

    tdir = tempfile.TemporaryDirectory()
    tmpoutpref = str(Path(tdir.name, outpref))

    # Make the mesh image
    scr_in = Path(msh_loc_in) / 'msh2png.scr'
    scr_out = Path(msh_loc_out) / (msh_name + '_msh2png.scr')
    fn_out = str(tmpoutpref) + '-mesh_nodes'
    f2f_with_subs(scr_in, scr_out, {'tmp': fn_out})

    cmd = [gmsh_exe, msh_name + '.msh', scr_out.name]
    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    return fn_out


def plot_mesh(msh_loc_in, msh_loc_out, msh_name, outpref):
    """Visualise mesh with gmsh and save to a file."""

    # Manipulate scripts in backend/fortran/build
    # Writes final png file to user directory

    fout_wire = plot_gmesh_wireframe(msh_loc_in, msh_loc_out, msh_name, outpref)
    fout_mesh = plot_gmesh_mesh(msh_loc_in, msh_loc_out, msh_name, outpref)

    nbapp = numbat.NumBATApp()
    #gmsh_exe = nbapp.path_gmsh()

    outprefix = nbapp.outpath(outpref)

    # tdir = tempfile.TemporaryDirectory()
    # tmpoutpref = str(Path(tdir.name, outpref))

    # # Make the wire frame image
    # scr_in = Path(msh_loc_in) / 'geo2png.scr'
    # scr_out = Path(msh_loc_out) / (msh_name + '_geo2png.scr')
    # f2f_with_subs(scr_in, scr_out, {'tmp': str(tmpoutpref) + '-entities'})
    # cmd = [gmsh_exe, msh_name + '.geo', scr_out.name]
    # run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    # # Make the mesh image
    # scr_in = Path(msh_loc_in) / 'msh2png.scr'
    # scr_out = Path(msh_loc_out) / (msh_name + '_msh2png.scr')
    # f2f_with_subs(scr_in, scr_out, {'tmp': str(tmpoutpref) + '-mesh_nodes'})

    # cmd = [gmsh_exe, msh_name + '.msh', scr_out.name]
    # run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    # Join the two images
    plottools.join_figs([fout_wire + '.png',
                        fout_mesh +'.png',],
                        str(outprefix)+'-mesh.png',
                        #clip=(60,50,60,50)
                        trimwhite=True
                        )


def check_mesh(msh_loc_out, msh_name):
    """Visualise geometry and mesh with gmsh."""

    nbapp = numbat.NumBATApp()
    gmsh_exe = str(nbapp.path_gmsh())

    gmsh_cmd = [gmsh_exe, f'{msh_loc_out}/{msh_name}.geo']
    run_subprocess(gmsh_cmd, 'Gmsh', cwd=msh_loc_out)

    gmsh_cmd = [gmsh_exe, f'{msh_loc_out}/{msh_name}.msh']
    run_subprocess(gmsh_cmd, 'Gmsh', cwd=msh_loc_out)
