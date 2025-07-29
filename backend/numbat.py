# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

"""
NumBAT main application module.

This module provides the main application singleton, handles environment setup, and provides interfaces for running and managing NumBAT simulations, including output management, version checking, and integration with the Fortran backend.

Key classes and functions:

- :class:`NumBATApp` -- Singleton class for managing global NumBAT state and configuration.
- :func:`NumBATPlotPrefs` -- Returns the global plot preferences object.
- :func:`load_simulation` -- Loads a saved simulation from disk.
- :func:`version` -- Returns the NumBAT version string.
"""

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
import structure
from modecalcs import Simulation
import meshing.templates as mshtemplates
from plotting.plotprefs import PlotPrefs


import nbversion

_envvar_gmsh_binary = 'NUMBAT_PATH_GMSH'
_envvar_numbat_root = 'NUMBAT_ROOT_DIR'

g_numbat_plot_prefs = None


def _confirm_file_exists(nm, filetype, fn, envvar=''):
    """
    Check that a file exists at the given path, and exit with an error if not.

    Args:
        nm (str): Name of the resource (for error messages).
        filetype (str): Type of file (for error messages).
        fn (str): Path to the file.
        envvar (str): Environment variable to suggest if missing.
    """
    if not Path(fn).exists():
        s=f"NumBAT can't find the {nm} {filetype} at {fn}."
        if envvar:
            s += f'You may need to set the environment variable {envvar}.'
        reporting.report_and_exit(s)

class _SingletonMeta(type):
    """
    Metaclass for singleton pattern. Ensures only one instance exists and provides _is_instantiated().
    """
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args, **kwargs)
            cls._instantiated = True
        return cls._instances[cls]

class NumBATApp(metaclass=_SingletonMeta):
    """
    Singleton class for managing the NumBAT application state, configuration, and output paths.
    Provides methods for output management, environment setup, and version checking.
    Use NumBATApp() to get the singleton instance.
    """
    _instantiated = False

    @classmethod
    def _is_instantiated(cls):
        """Return True if the singleton has been instantiated."""
        return cls._instantiated

    def __init__(self, outprefix='nbtmp', outdir='.',
                 user_settings_file='',
                 ignore_user_settings=False
                 ):
        """
        Initialize the NumBAT application singleton.

        Args:
            outprefix (str): Output file prefix.
            outdir (str): Output directory.
            user_settings_file (str): Path to user settings file.
            ignore_user_settings (bool): If True, ignore user settings file.
        """


        global g_numbat_plot_prefs

        reporting.init_logger()
        self._start_time=time.time()
        self._check_versions()


        self._paths={}
        self._codedir = Path(__file__).parents[0]


        self._setup_paths()


        self._outprefix=outprefix
        self.set_outdir(outdir)

        #self._user_settings={}


        g_numbat_plot_prefs = PlotPrefs(user_settings_file, ignore_user_settings)

        mshtemplates.initialise_waveguide_templates(self)

    def make_structure(self, *args, **kwargs):
        """
        Create and return a Structure object for the simulation.
        """
        return structure.Structure(*args, direct_call=False, **kwargs)

    def is_linux(self):
        """
        Return True if running on Linux.
        """
        return platform.system()=='Linux'

    def is_macos(self):
        """
        Return True if running on macOS.
        """
        return platform.system()=='Darwin'

    def is_windows(self):
        """
        Return True if running on Windows.
        """
        return platform.system()=='Windows'

    def can_multiprocess(self):
        """
        Return True if the platform supports multiprocessing (Linux only).
        """
        return self.is_linux()

    def plotfile_ext(self):
        """
        Return the current plot file extension as a string.
        """
        return g_numbat_plot_prefs._plot_extension

    # Only needed to change prefix at some point
    def set_outprefix(self, s): # change to getter/setter
        """
        Set the global output file prefix.

        Args:
            s (str): Output prefix.
        """
        self._outprefix=s

    def outprefix(self):
        """
        Return the current global output prefix as a string.
        """
        return self._outprefix


    def set_outdir(self, outdir):
        """
        Set the current global output directory, creating it if necessary.

        Args:
            outdir (str): Output directory path.
        """
        self._outdir = outdir
        path=Path(outdir)

        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        elif not path.is_dir():
            reporting.report_and_exit(f'Output path {outdir} exists and is not a directory.')

    def outdir(self):
        """
        Return the current global output directory as a string.
        """
        return self._outdir

    def outpath(self, prefix=''):
        """
        Return the full output path for a given prefix.

        Args:
            prefix (str): Output prefix (optional).
        Returns:
            str: Full output path.
        """
        outpref = prefix if prefix else self._outprefix
        return str(Path(self._outdir, outpref))

    def outdir_fields_path(self, prefix=''):
        """
        Return the output directory path for modal fields, creating it if necessary.

        Args:
            prefix (str): Output prefix (optional).
        Returns:
            Path: Path object for the fields directory.
        """
        outpref = prefix if prefix else self._outprefix
        pout = Path(self._outdir, outpref+'-fields')

        if not pout.exists():
            pout.mkdir(parents=True, exist_ok=True)

        return pout


    def path_gmsh(self):
        """
        Return the path to the gmsh executable.
        """
        return self._paths['gmsh']

    def path_mesh_templates(self):
        """
        Return the path to the mesh templates directory.
        """
        return Path(self._codedir, 'msh')

    def final_report(self, outprefix=''):
        """
        Return a string summarizing the simulation run and any warnings.
        Optionally writes warnings to a file.

        Args:
            outprefix (str): Output prefix for warnings file (optional).
        Returns:
            str: Summary report string.
        """
        dt=time.time()-self._start_time
        s_dt = datetime.timedelta(seconds=round(dt))

        s=f"\nNumBAT calculations concluded. \n Simulation time: {s_dt}.\n"
        s+=reporting.report_warnings()+'\n'

        if outprefix:
            with open(outprefix+'_warnings.txt', 'w') as fout:
                fout.write(s)

        return s


    def _wg_structure_help(self, inc_shape):
        """
        Print help for waveguide structure templates.
        """
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
        nbrootdir = os.getenv(_envvar_numbat_root, str(Path(__file__).resolve().parents[3]))

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
        """
        Check that the Python and Fortran versions match the required NumBAT version.
        Exits with an error if not.
        """
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
        """
        Return the NumBAT version string.
        """
        return nbversion.NUMBAT_VERSION_STR_MMM


# # always returns the singleton NumBATApp object
# def NumBATApp(outprefix='', outdir='.',
#               user_settings_file='',
#               ignore_user_settings=False):
#     """
#     Returns the singleton NumBATApp object on every call.

#     Args:
#         outprefix (str): Output file prefix.
#         outdir (str): Output directory.
#         user_settings_file (str): Path to user settings file.
#         ignore_user_settings (bool): If True, ignore user settings file.
#     Returns:
#         _NumBATApp: The singleton application instance.
#     """
#     nba = _NumBATApp(outprefix, outdir, user_settings_file, ignore_user_settings)
#     return nba


def _is_NumBATApp_created():
    """
    Returns True if the NumBATApp singleton has been instantiated, False otherwise.
    """
    return NumBATApp._is_instantiated()

def _assert_numbat_object_created():
    """
    Ensure that the NumBAT singleton object has been created before using other functions.
    Exits with an error if not.
    """
    if not _is_NumBATApp_created():
        reporting.report_and_exit('You must create a NumBAT object before calling any other NumBAT functions.  See the tutorials for examples.')

def NumBATPlotPrefs():
    """
    Returns the global NumBAT plot preferences object.
    """
    return g_numbat_plot_prefs


def version():
    """
    Returns the NumBAT version string.
    """
    return nbversion.NUMBAT_VERSION_STR_MMM

def load_simulation(prefix):
    """
    Loads a saved Simulation object from disk using the given prefix.

    Args:
        prefix (str): Prefix for the simulation file.
    Returns:
        Simulation: The loaded simulation object.
    """
    return Simulation.load_simulation(prefix)
