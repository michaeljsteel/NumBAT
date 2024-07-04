import os 
import sys 
import importlib
from pathlib import Path
import subprocess
import platform

def is_linux():
    return platform.system()=='Linux'

def is_macos():
    return platform.system()=='Darwin'

def is_windows():
    return platform.system()=='Windows'


def nb_platform():
    if is_linux(): return 'Linux'
    elif is_macos(): return 'MacOS Darwin'
    elif is_windows(): return 'Windows'
    else:
        return 'Unknown platform'

def check_gmsh(nbapp, run_gmsh):

    passed=False
    print(f'\nChecking gmsh')

    envar = 'NUMBAT_PATH_GMSH'
    if envar in os.environ:
        print('  Path variable NUMBAT_PATH_GMSH is set.')
    else:
        print('  Path variable NUMBAT_PATH_GMSH is not set, using default location.')

    path_gmsh = nbapp.path_gmsh()
    print(f'  Path expected: {path_gmsh}')

    s_exists = 'yes' if Path(path_gmsh).exists() else 'no'
    print(f'  File exists: {s_exists}')

    if s_exists == 'yes' and run_gmsh:
        print('\n  Attempting to open a Gmsh windows.')
        print('  Close the window if it appears...')
        procret = subprocess.run(path_gmsh)
        if procret.returncode == 0:
            print('  Executes ok: yes')
            passed=True
        else:
            print('  Executes ok: no')

    if passed:
        print('Gmsh tests completed: Pass')
    else:
        print('Gmsh tests completed: Fail')

    return passed

def check_asymptote(nbapp):
    print('\nChecking asymptote')
    print('   no tests yet')

    passed=True
    return passed

def check_nb_fortran():
    passed=False

    print('\nChecking NumBAT core Fortran module')

    s_lib_nbfort = 'fortran.nb_fortran'
    s_lib_nbfort_full = 'fortran/nb_fortranb.pyd' if is_windows() else 'fortran/nb_fortranb.so' 

    mod = None
    try:
        mod = importlib.import_module(s_lib_nbfort)
    except ModuleNotFoundError:
        print(f"  Can't find module file {s_lib_nbfort_full}.")
    except ImportError as err:
        print(f"Module failed to load: {err=}.\n One or more shared libraries can't be located.")
    else:
        print(f"  Successfully loaded.")
        passed = True

    print('\n\n')
    
    return passed, mod

def check_helper_apps(nbapp):

    pass_gmsh = check_gmsh(nbapp, False)
    pass_asy = check_asymptote(nbapp)

    return pass_gmsh and pass_asy and pass_nb_fortran



def do_main(argv):

    print('\n\n')
    print('NumBAT install tester')
    print('---------------------')

    print(f'\n\nPerforming tests for platform: {nb_platform()}.')
    print('\n\n')

    # Can't create NumBATApp without the main module

    pass_nbfort, mod_nb = check_nb_fortran()

    if not pass_nbfort: 
        sys.exit(1)
    
    import numbat
    nbapp = numbat.NumBATApp()



    pass_apps = check_helper_apps(nbapp)

    if pass_apps and pass_nbfor:
        print('\n\nAll tests passed. You are ready to NumBAT!')
    else:
        print('\n\nPlease attend to the failed tests before starting to work with NumBAT.')

    print('\n\n')

if __name__=='__main__':
    do_main(sys.argv)
