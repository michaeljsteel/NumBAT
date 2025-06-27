project(
  'NumBAT Fortran backend',
  ['fortran', 'c'],
  version: '2.1.2', # this has to match the version in nbversion.py
  meson_version: '>1.4.0',
)

version_output = run_command(find_program('python3'), files('../nbversion.py'), check: true)
project_version = version_output.stdout().strip()

c_comp = meson.get_compiler('c')
f_comp = meson.get_compiler('fortran')

#message('@0@'.format(c_comp.get_id()))
#message('@0@'.format(f_comp.get_id()))

message('C compiler:', c_comp.get_id())
message('F compiler:', f_comp.get_id())

message('host:', host_machine.system())

libs_extra = []

plat_pref = host_machine.system() + '_'

dir_usrlocal = get_option(plat_pref + 'dir_nb_usrlocal')
dir_usrlocal_lib = dir_usrlocal / 'lib'

dir_usrlocal_inc = dir_usrlocal / 'include'

if host_machine.system() == 'linux'

  if c_comp.get_id() == 'gcc'
    # Get the compiler version as a version object
    gcc_version = c_comp.version()

    # Define the minimum required version
    min_gcc_version = '13.0'

    # Compare the versions
    if gcc_version.version_compare('< ' + min_gcc_version)
      error('GCC version ' + gcc_version + ' is too old. NumBAT requires a minimum version of ' + min_gcc_version + '.')
    endif
  endif

  use_mkl_gcc = false
  inc_nblocal = include_directories([dir_usrlocal_inc,])
  dir_oneapi = get_option(plat_pref + 'dir_oneapi')
  dir_mkl = dir_oneapi / 'mkl/latest/lib'

  if use_mkl_gcc
    dep_blas = c_comp.find_library('mkl_rt', dirs: dir_mkl, header_include_directories: inc_nblocal)
  else
    dep_blas = dependency('blas', required: true)
  endif

elif host_machine.system() == 'darwin'
  inc_nblocal = include_directories([dir_usrlocal_inc])
  dep_blas = [
    c_comp.find_library('lapack', dirs: '/opt/local/lib/lapack'),
    c_comp.find_library('blas', dirs: '/opt/local/lib/lapack'),
  ]

else
  # Windows

  dir_oneapi = get_option(plat_pref + 'dir_oneapi')
  dir_mkl = dir_oneapi / 'mkl/latest/lib'
  inc_nblocal = include_directories([dir_usrlocal_inc])
  dep_blas = c_comp.find_library('mkl_rt', dirs: dir_mkl, header_include_directories: inc_nblocal)
endif

if c_comp.get_id() == 'intel-llvm' # using intel compiler
  lib_oneapi = c_comp.find_library('ifport', dirs: dir_mkl)
  libs_extra = [lib_oneapi]
endif

lib_arpack = c_comp.find_library('arpack', dirs: dir_usrlocal_lib)
lib_umfpack = c_comp.find_library('umfpack', dirs: dir_usrlocal_lib)

# enforce implicit none
if host_machine.system() == 'windows'
  add_global_arguments(
    ['/warn:declarations'],
    language: 'fortran',
  )
  add_global_arguments(
    ['/warn:no-format', '/warn:no-misleading-indentation'],
    language: 'c',
  )

elif c_comp.get_id() == 'intel-llvm'
  # intel on linux
  add_global_arguments(
    ['-implicitnone'],
    language: 'fortran',
  )
  add_global_arguments(
    ['-Wno-format', '-Wno-misleading-indentation'],
    language: 'c',
  ) # numpy f2py has bad format strings that do no harm

else
  # gcc on linux and mac
  #add_global_arguments(['-fimplicit-none', '-std=f2023', '-fall-intrinsics'], language: 'fortran',)
  add_global_arguments(
    ['-fimplicit-none'],
    language: 'fortran',
  )
  add_global_arguments(
    ['-Wno-format'],
    language: 'c',
  ) # numpy f2py has bad format strings that do no harm
endif

if host_machine.system() == 'windows'
  add_global_arguments(
    ['-fpp'],
    language: 'fortran',
  )
  add_global_arguments(
    ['-D UPPERCASE_FORTRAN', '-D NO_APPEND_FORTRAN', '-D __INTEL_COMPILER'],
    language: 'c',
  )
  add_global_arguments(['-D __INTEL_COMPILER'], language: 'cpp')
else
  add_global_arguments(
    ['-cpp'],
    language: 'fortran',
  )
endif

add_global_arguments(
  [
    '-DZLONG',
  ],
  language: 'c',
)

comp_f2py = find_program('f2py')

subdir('sswrap')
subdir('src')

pymod = import('python')
py_installation = pymod.find_installation()
#py_exe = find_program('python', py_installation.path() / 'bin')
py_exe = py_installation.full_path()
dep_py = py_installation.dependency()

# F2PY details
r = run_command(py_exe, '-c', 'import numpy; print(numpy.get_include())', check: true)
dir_numpy = r.stdout().strip()
r = run_command(
  py_exe,
  '-c', 'import numpy.f2py; print(numpy.f2py.get_include())',
  check: true,
)
dir_f2py = r.stdout().strip()
src_f2py_fortobj = dir_f2py / 'fortranobject.c'

incs_f2py = include_directories([dir_numpy, dir_f2py])

message([dir_numpy, dir_f2py])

incs = [inc_nblocal, incs_f2py]

srcs = [
  src_f2py_fortobj,
  nb_pyf,
  nb_fortmodc,
  nb_pyinterface_fsrcs,
  nb_module_fsrcs,
  nb_impl_fsrcs,
  src_sswrap,
]

nblib = shared_module(
  'nb_fortran',
  srcs,
  include_directories: incs,
  name_prefix: '',
  override_options: ['b_lundef=false'],
  dependencies: [lib_arpack, lib_umfpack, dep_blas, libs_extra, dep_py],
)
