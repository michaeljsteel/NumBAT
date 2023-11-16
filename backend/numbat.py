import platform
import reporting
import os
import shutil
import time
import datetime

import objects


_evar_gmsh_path = 'NUMBAT_PATH_GMSH'

def _confirm_file_exists(nm, path, evar=''):
    if not os.path.exists(path):
        s=f"Can't find the {nm} executable at {path}."
        if evar:
            s += f'You may need to set the environment variable {evar}.'
        reporting.report_and_exit(s)

class _NumBATApp(object):
    my_num_instances = 0

    __instance = None

    #def __new__(cls):
    #    print('seeking nba')
    #    if not hasattr(cls, 'instance'):
    #        cls.instance = super(NumBATApp, cls).__new__(cls)
    #        print('newing nba with insta', cls.instance)
    #    return cls.instance

    def __init__(self, prefix='nbtmp'):
        #if _NumBATApp.my_num_instances:
        #    reporting.report_and_exit('You may only create a single NumBAT object.')

        _NumBATApp.__instance = self

        _NumBATApp.my_num_instances += 1

        self._prefix=prefix
        self._paths={}
        self._start_time=time.time()

        self._check_versions()
        self._setup_paths()
        reporting.init_logger()

    @staticmethod
    def get_instance(prefix=''):
        if _NumBATApp.__instance is None:
            _NumBATApp(prefix)  # instance gets attached inside __init__
        return _NumBATApp.__instance

    def is_linux(self):
        return platform.system()=='Linux'

    def is_macos(self):
        return platform.system()=='Darwin'

    def is_windows(self):
        return platform.system()=='Windows'

    def set_prefix(self, s): # change to getter/setter
        self._prefix=s

    def prefix(self):
        return self._prefix

    def path_gmsh(self):
        return self._paths['gmsh']

    def final_report(self, prefix=''):
        dt=time.time()-self._start_time
        s_dt = datetime.timedelta(seconds=round(dt))

        s=f"\nNumBAT calculations concluded. \n Simulation time: {s_dt}.\n"
        s+=reporting.report_warnings()+'\n'

        if prefix:
            with open(prefix+'_warnings.txt', 'w') as fout:
                fout.write(s)

        return s

    def make_structure(self, *args, **kwargs):
        return objects.Structure(*args, direct_call=False, **kwargs)


    def _setup_paths(self):
        if self.is_linux():
            path = shutil.which('gmsh')
            self._paths['gmsh'] = os.environ.get(_evar_gmsh_path, path)

        elif self.is_macos():
            self._paths['gmsh'] = os.environ.get(_evar_gmsh_path,
                                                 '/Applications/Gmsh.app/Contents/MacOS/gmsh')
        else:
            pass

        _confirm_file_exists('Gmsh', self._paths['gmsh'], _evar_gmsh_path)


    def _check_versions(self):
        pyver = platform.python_version_tuple()

        if int(pyver[0])<3 or int(pyver[1])<7:
            reporting.report_and_exit('NumBAT must be run with a Python version of 3.6 or later.')


def NumBATApp(prefix=''):
    nba = _NumBATApp.get_instance(prefix)  # always returns the singleton NumBATApp object
    return nba


def assert_numbat_object_created():
    if _NumBATApp.my_num_instances != 1:
        reporting.report_and_exit('In NumBAT 2.0, you must now create a NumBAT object before calling any other NumBAT functions.  See the tutorials for examples.')
