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
from mode_calcs import Simulation


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
        #if _NumBATApp.my_num_instances:
        #    reporting.report_and_exit('You may only create a single NumBAT object.')

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
        objects.initialise_waveguide_templates(self)

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

    def outpath_fields(self, prefix=''):
        '''Returns a string composed of the current output directory for modal fields and the given prefix.'''
        outpref = prefix if prefix else self._outprefix
        return str(Path(self._outdir, outpref+'-fields'))


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
        objects.print_waveguide_help(inc_shape)

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


# always returns the singleton NumBATApp object
def NumBATApp(outprefix='', outdir='.'):
    '''Returns the same singleton NumBATApp object on every call.'''

    nba = _NumBATApp.get_instance(outprefix, outdir)
    return nba


def assert_numbat_object_created():
    if _NumBATApp.my_num_instances != 1:
        reporting.report_and_exit('In NumBAT 2.0, you must now create a NumBAT object before calling any other NumBAT functions.  See the tutorials for examples.')



#TODO: move this to plotting.py and
# make the NumBATPlotPrefs call return a reference to object in plotting.py
class _NumBATPlotPrefs:
    def __init__(self):
        self.cmap_field_signed = 'seismic'
        self.cmap_field_unsigned = 'OrRd'
        self.cmap_ref_index = 'cool'

def NumBATPlotPrefs():
    return _NumBATPlotPrefs()

def load_simulation(prefix):
    return Simulation.load_simulation(prefix)
