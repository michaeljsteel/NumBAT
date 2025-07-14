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

    outprefix = nbapp.outpath(outpref)

    tdir = tempfile.TemporaryDirectory()
    tmpoutpref = str(Path(tdir.name, outpref))

    #tmpoutpref = str(Path(tempfile.TemporaryFile(), outpref))

    # Make the wire frame image
    fn_in = Path(msh_loc_in) / 'geo2png.scr'
    fn_out = Path(msh_loc_out) / (msh_name + '_geo2png.scr')
    f2f_with_subs(fn_in, fn_out, {'tmp': str(tmpoutpref) + '-entities'})
    cmd = [gmsh_exe, msh_name + '.geo', fn_out.name]
    print(outprefix, tmpoutpref, fn_out)
    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

def plot_gmesh_mesh(msh_loc_in, msh_loc_out, msh_name, outpref):

    nbapp = numbat.NumBATApp()
    gmsh_exe = nbapp.path_gmsh()

    outprefix = nbapp.outpath(outpref)

    tdir = tempfile.TemporaryDirectory()
    tmpoutpref = str(Path(tdir.name, outpref))

    # Make the mesh image
    fn_in = Path(msh_loc_in) / 'msh2png.scr'
    fn_out = Path(msh_loc_out) / (msh_name + '_msh2png.scr')
    f2f_with_subs(fn_in, fn_out, {'tmp': str(tmpoutpref) + '-mesh_nodes'})

    cmd = [gmsh_exe, msh_name + '.msh', fn_out.name]
    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)



def plot_mesh(msh_loc_in, msh_loc_out, msh_name, outpref):
    """Visualise mesh with gmsh and save to a file."""

    # Manipulate scripts in backend/fortran/build
    # Writes final png file to user directory


    #plot_gmesh_wireframe(msh_loc_in, msh_loc_out, msh_name, outpref)
    #plot_gmesh_mesh(msh_loc_in, msh_loc_out, msh_name, outpref)

    nbapp = numbat.NumBATApp()
    gmsh_exe = nbapp.path_gmsh()

    outprefix = nbapp.outpath(outpref)

    tdir = tempfile.TemporaryDirectory()
    tmpoutpref = str(Path(tdir.name, outpref))

    # Make the wire frame image
    fn_in = Path(msh_loc_in) / 'geo2png.scr'
    fn_out = Path(msh_loc_out) / (msh_name + '_geo2png.scr')
    f2f_with_subs(fn_in, fn_out, {'tmp': str(tmpoutpref) + '-entities'})
    cmd = [gmsh_exe, msh_name + '.geo', fn_out.name]
    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    # Make the mesh image
    fn_in = Path(msh_loc_in) / 'msh2png.scr'
    fn_out = Path(msh_loc_out) / (msh_name + '_msh2png.scr')
    f2f_with_subs(fn_in, fn_out, {'tmp': str(tmpoutpref) + '-mesh_nodes'})

    cmd = [gmsh_exe, msh_name + '.msh', fn_out.name]
    run_subprocess(cmd, 'Gmsh', cwd=msh_loc_out)

    # Join the two images
    plottools.join_figs([tmpoutpref+'-entities.png',
                        tmpoutpref+'-mesh_nodes.png',],
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
