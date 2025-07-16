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
import sys
import shutil
import time
import datetime
from pathlib import Path

try:
    from fortran import nb_fortran
except ImportError:
    print("""
          NumBAT can't load the backend Fortran module that implements the finite element solvers!
          Are you sure you have compiled the module?
          See the installation chapter of the manual for details.
          """)
    sys.exit(1)

import reporting
import objects
from modecalcs import Simulation
import meshing.templates as mshtemplates
from plotting.plotprefs import PlotPrefs


import nbversion

_envvar_gmsh_binary = 'NUMBAT_PATH_GMSH'
_envvar_numbat_root = 'NUMBAT_ROOT_DIR'

g_numbat_plot_prefs = None


def _confirm_file_exists(nm, filetype, fn, envvar=''):
    if not Path(fn).exists():
        s=f"NumBAT can't find the {nm} {filetype} at {fn}."
        if envvar:
            s += f'You may need to set the environment variable {envvar}.'
        reporting.report_and_exit(s)

class _NumBATApp:
    _num_instances = 0

    __instance = None

    @staticmethod
    def get_instance(outprefix='', outdir='',
                     user_settings_file='',
                     ignore_user_settings=False):

        # instance gets attached inside __init__
        if _NumBATApp.__instance is None:
            _NumBATApp(outprefix, outdir, user_settings_file, ignore_user_settings)

        return _NumBATApp.__instance

    def __init__(self, outprefix='nbtmp', outdir='.',
                 user_settings_file='',
                 ignore_user_settings=False
                 ):


        global g_numbat_plot_prefs

        _NumBATApp.__instance = self
        _NumBATApp._num_instances += 1

        reporting.init_logger()
        self._start_time=time.time()
        self._check_versions()


        self._paths={}
        self._codedir = Path(__file__).parents[0]


        self._setup_paths()


        self._outprefix=outprefix
        self.set_outdir(outdir)

        self._user_settings={}


        #self._plot_extension = '.png'
        #self._plot_extension = '.pdf'

        g_numbat_plot_prefs = PlotPrefs(user_settings_file, ignore_user_settings)




        mshtemplates.initialise_waveguide_templates(self)



    def plotfile_ext(self):
        return g_numbat_plot_prefs._plot_extension


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
        # if not Path(self._outdir).is_dir():
        #     try:
        #         Path(self._outdir).mkdir()
        #     except OSError as ex:
        #         reporting.report_and_exit(f"Can't open output directory {self._outdir}: {str(ex)}")

         # location of top level numbat tree containing other libraries etc. Mainly for windows
        # this seems a bit flaky depending on installations
        nbrootdir =  os.getenv(_envvar_numbat_root,
                               default=Path(__file__).resolve().parents[3])

        # paths to other tools
        if self.is_linux():
            gmpath = shutil.which('gmsh')
            self._paths['gmsh'] = Path(os.getenv(_envvar_gmsh_binary, default=gmpath))

        if self.is_windows():
            gmpath = Path(nbrootdir, 'usr_local/packages/gmsh/gmsh.exe')

        elif self.is_macos():
            gmpath = '/Applications/Gmsh.app/Contents/MacOS/gmsh'

        self._paths['gmsh'] = Path(os.getenv(_envvar_gmsh_binary, default=gmpath))

        _confirm_file_exists('Gmsh', 'executable', self._paths['gmsh'], _envvar_gmsh_binary)


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
def NumBATApp(outprefix='', outdir='.',
              user_settings_file='',
              ignore_user_settings=False):
    '''Returns the same singleton NumBATApp object on every call.'''

    nba = _NumBATApp.get_instance(outprefix, outdir, user_settings_file, ignore_user_settings)
    return nba


def assert_numbat_object_created():
    if _NumBATApp._num_instances != 1:
        reporting.report_and_exit('You must create a NumBAT object before calling any other NumBAT functions.  See the tutorials for examples.')








def NumBATPlotPrefs():
    return g_numbat_plot_prefs


def version():
        return nbversion.NUMBAT_VERSION_STR_MMM

def load_simulation(prefix):
    return Simulation.load_simulation(prefix)
