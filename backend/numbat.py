# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import platform
import os
import shutil
import time
import datetime
from pathlib import Path

import reporting
import objects
from modecalcs import Simulation
import meshing.templates as mshtemplates

from fortran import nb_fortran
import nbversion

_envvar_gmsh_binary = 'NUMBAT_PATH_GMSH'
_envvar_numbat_root = 'NUMBAT_ROOT_DIR'

def _confirm_file_exists(nm, fn, envvar=''):
    if not Path(fn).exists():
        s=f"Can't find the {nm} executable at {fn}."
        if envvar:
            s += f'You may need to set the environment variable {envvar}.'
        reporting.report_and_exit(s)

class _NumBATApp:
    my_num_instances = 0

    __instance = None

    #def __new__(cls):
    #    print('seeking nba')
    #    if not hasattr(cls, 'instance'):
    #        cls.instance = super(NumBATApp, cls).__new__(cls)
    #        print('newing nba with insta', cls.instance)
    #    return cls.instance

    def __init__(self, outprefix='nbtmp', outdir='.'):

        _NumBATApp.__instance = self

        _NumBATApp.my_num_instances += 1

        self._outprefix=outprefix
        self._outdir=outdir
        self._paths={}
        self._start_time=time.time()
        self._codedir = Path(__file__).parents[0]

        # location of top level numbat tree containing other libraries etc. Mainly for windows
        # this seems a bit flaky depending on installations
        self._nbrootdir =  os.getenv(_envvar_numbat_root,
                                          default=Path(__file__).resolve().parents[3])

        self._plot_extension = '.png'
        #self._plot_extension = '.pdf'

        self._check_versions()
        self._setup_paths()
        reporting.init_logger()
        mshtemplates.initialise_waveguide_templates(self)

    @staticmethod
    def get_instance(outprefix='', outdir=''):
        if _NumBATApp.__instance is None:
            _NumBATApp(outprefix, outdir)  # instance gets attached inside __init__
        return _NumBATApp.__instance

    def plotfile_ext(self):
        return self._plot_extension

    def is_linux(self):
        return platform.system()=='Linux'

    def is_macos(self):
        return platform.system()=='Darwin'

    def is_windows(self):
        return platform.system()=='Windows'

    def can_multiprocess(self):
        return self.is_linux()

    # Only needed to change prefix at some point
    def set_outprefix(self, s): # change to getter/setter
        self._outprefix=s

    def outprefix(self):
        """Return the current global output prefix as a string."""
        return self._outprefix

    def outdir(self):
        """Return the current global output directory as a string."""
        return self._outdir

    def set_outdir(self, outdir):
        """Set the current global output directory."""
        self._outdir = outdir
        path=Path(outdir)

        if not path.exists():
            os.makedirs(path)
        elif not path.is_dir():
            reporting.report_and_exit(f'Output path {outdir} exists and is not a directory.')

    def outpath(self, prefix=''):
        '''Returns a string composed of the current output directory and the given prefix.'''

        outpref = prefix if prefix else self._outprefix
        return str(Path(self._outdir, outpref))

    def outdir_fields_path(self, prefix=''):
        '''Returns a string composed of the current output directory for modal fields and the given prefix.

        Creates the directory if it does not exist.'''
        outpref = prefix if prefix else self._outprefix
        pout = Path(self._outdir, outpref+'-fields')

        if not pout.exists():
            pout.mkdir()

        return pout


    def path_gmsh(self):
        return self._paths['gmsh']

    def path_mesh_templates(self):
        return Path(self._codedir, 'msh')

    def final_report(self, outprefix=''):
        dt=time.time()-self._start_time
        s_dt = datetime.timedelta(seconds=round(dt))

        s=f"\nNumBAT calculations concluded. \n Simulation time: {s_dt}.\n"
        s+=reporting.report_warnings()+'\n'

        if outprefix:
            with open(outprefix+'_warnings.txt', 'w') as fout:
                fout.write(s)

        return s

    def make_structure(self, *args, **kwargs):
        return objects.Structure(*args, direct_call=False, **kwargs)

    def wg_structure_help(self, inc_shape):
        mshtemplates.print_waveguide_help(inc_shape)

    def _setup_paths(self):
        # numbat paths
        if not Path(self._outdir).is_dir():
            try:
                Path(self._outdir).mkdir()
            except OSError as ex:
                reporting.report_and_exit(f"Can't open output directory {self._outdir}: {str(ex)}")

        # paths to other tools
        if self.is_linux():
            gmpath = shutil.which('gmsh')
            self._paths['gmsh'] = Path(os.getenv(_envvar_gmsh_binary, default=gmpath))

        if self.is_windows():
            gmpath = Path(self._nbrootdir, 'usr_local/packages/gmsh/gmsh.exe')
            self._paths['gmsh'] = Path(os.getenv(_envvar_gmsh_binary, default=gmpath))

        elif self.is_macos():
            self._paths['gmsh'] = Path(os.getenv(_envvar_gmsh_binary, default=
                                            '/Applications/Gmsh.app/Contents/MacOS/gmsh'))
        else:
            pass

        _confirm_file_exists('Gmsh', self._paths['gmsh'], _envvar_gmsh_binary)


    def _check_versions(self):
        pyver = platform.python_version_tuple()

        if int(pyver[0])<3 or int(pyver[1])<10:
            reporting.report_and_exit('NumBAT must be run with a Python version of 3.11 or later.')

        s_fortver = str(nb_fortran.get_fortran_compiled_nb_version(), "utf-8")
        version_match = s_fortver == nbversion.NUMBAT_VERSION_STR_MMM
        if not version_match:
            reporting.report_and_exit(
                f'NumBAT has detected different versions for the Python (ver={nbversion.NUMBAT_VERSION_STR_MMM}) and Fortran components (ver={s_fortver}).\n Most likely, you need to rebuild the Fortran module.'
                )

    def version(self):
        return nbversion.NUMBAT_VERSION_STR_MMM


# always returns the singleton NumBATApp object
def NumBATApp(outprefix='', outdir='.'):
    '''Returns the same singleton NumBATApp object on every call.'''

    nba = _NumBATApp.get_instance(outprefix, outdir)
    return nba


def assert_numbat_object_created():
    if _NumBATApp.my_num_instances != 1:
        reporting.report_and_exit('In NumBAT 2.0, you must now create a NumBAT object before calling any other NumBAT functions.  See the tutorials for examples.')



#TODO: move this to plottools.py and
# make the NumBATPlotPrefs call return a reference to object in plottools.py
class _NumBATPlotPrefs:
    """Selection of color scales and line properties for different field plots."""

    # Add color combinations here
    # See https://matplotlib.org/stable/users/explain/colors/colormaps.html

                  #  signed, unsigned, vector arrow
    color_tup_1 = ('seismic', 'OrRd', 'dimgray')
    color_tup_2 = ('PRGn', 'GnBu', 'brown')
    color_tup_3 = ('BrBG', 'YlOrBr', 'black')
    color_tup_4 = ('PuOr', 'YlOrBr', 'dimgray')
    color_tup_5 = ('coolwarm', 'Reds', 'dimgray')
    color_tup_6 = ('PiYG', 'GnBu', 'dimgray')
    color_tup_7 = ('RdYlGn', 'GnBu', 'dimgray')




    def __init__(self):

        # Select color combinations here
        coltup_em = self.color_tup_1
        coltup_ac = self.color_tup_7

        # electromagnetic plots
        (self.cmap_em_field_signed,        # Ex, Ey, Ez, Hx, Hy, Hz
            self.cmap_em_field_unsigned,   # |E|^2, |H|^2,
            self.vecplot_arrow_color_em
            ) = coltup_em


        # acoustic plots
        (self.cmap_ac_field_signed,        # ux, uy, uz
            self.cmap_ac_field_unsigned,   # |u|^2
            self.vecplot_arrow_color_ac
            ) = coltup_ac


        self.vecplot_arrow_scale = 10  # larger makes smaller arrows
        self.vecplot_arrow_linewidth = 0.005       # width of the shaft in fractions of plot width
        self.vecplot_arrow_headwidth = 2             # multiple of shaft width



        #self.vecplot_arrow_head = 2
        #self.vecplot_arrow_len = 2



        # colormap for refractive index plots
        self.cmap_ref_index = 'cool'


        # properties of waveguide boundary frames

        # EDGE_COLOR="gray"
        # EDGE_COLOR="dimgray"
        # EDGE_COLOR="black"
        # EDGE_COLOR="brown"
        EDGE_COLOR="darkred"

        self.WG_FRAME_EDGE_COLOR = EDGE_COLOR

        self.WG_FRAME_LINEWIDTH_WHOLEFIG = 0.75   # if field plot is whole figure
        self.WG_FRAME_LINEWIDTH_SUBFIG = 0.25      # if field plot is a part figure



    def cmap_field_signed(self, ftag):
        if ftag.is_EM():
            return self.cmap_em_field_signed
        else:
            return self.cmap_ac_field_signed

    def cmap_field_unsigned(self, ftag):
        if ftag.is_EM():
            return self.cmap_em_field_unsigned
        else:
            return self.cmap_ac_field_unsigned

    def vector_field_arrow_color(self, ftag):
        if ftag.is_EM():
            return self.vecplot_arrow_color_em
        else:
            return self.vecplot_arrow_color_ac


    def vector_field_arrow_scale(self):
        return self.vecplot_arrow_scale

    def vector_field_arrow_linewidth(self):
        return self.vecplot_arrow_linewidth


def NumBATPlotPrefs():
    return _NumBATPlotPrefs()

# def NumBAT_color_styles():


#     #EDGE_COLOR="mediumblue"
#     #EDGE_COLOR="darkblue"

#     styles = {'WG_FRAME_EDGE_COLOR': WG_FRAME_EDGE_COLOR,
#               'WG_FRAME_LINEWIDTH': WG_FRAME_LINEWIDTH,
#               }


#     return styles

def version():
        return nbversion.NUMBAT_VERSION_STR_MMM

def load_simulation(prefix):
    return Simulation.load_simulation(prefix)
