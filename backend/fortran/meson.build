project(
  'NumBAT Fortran backend',
  ['fortran', 'c'],
  version: '2.0.1',
)

ccomp = meson.get_compiler('c')
fcomp = meson.get_compiler('fortran')

message('@0@'.format(ccomp.get_id()))
message('@0@'.format(fcomp.get_id()))


libs_extra = []

plat_pref = host_machine.system() + '_'

dir_usrlocal = get_option(plat_pref + 'dir_nb_usrlocal')
dir_usrlocal_lib = dir_usrlocal / 'lib'

dir_usrlocal_inc = dir_usrlocal / 'include'
if host_machine.system() == 'linux'
dir_oneapi = get_option(plat_pref + 'dir_oneapi')
dir_mkl = dir_oneapi / 'mkl/latest/lib'

  inc_nblocal = include_directories([dir_usrlocal_inc,])
  lib_blas = ccomp.find_library('mkl_rt',
    dirs: dir_mkl,
    header_include_directories: inc_nblocal)

elif host_machine.system() == 'darwin'
  inc_nblocal = include_directories([dir_usrlocal_inc])
  lib_blas = [
    ccomp.find_library('lapack', dirs: '/opt/local/lib/lapack'),
    ccomp.find_library('blas', dirs: '/opt/local/lib/lapack'),
  ]

else # Windows

dir_oneapi = get_option(plat_pref + 'dir_oneapi')
dir_mkl = dir_oneapi / 'mkl/latest/lib'
  inc_nblocal = include_directories([dir_usrlocal_inc])
  lib_blas = ccomp.find_library('mkl_rt', dirs: dir_mkl,
  header_include_directories:inc_nblocal)
  endif

if ccomp.get_id() == 'intel-llvm'
  lib_oneapi = ccomp.find_library('ifport',
                                  dirs: dir_mkl
                                  )
  libs_extra = [lib_oneapi]
endif

lib_arpack = ccomp.find_library('arpack', dirs: dir_usrlocal_lib)
lib_umfpack = ccomp.find_library('umfpack', dirs: dir_usrlocal_lib)

if host_machine.system() == 'windows'
  add_global_arguments(['-fpp', ], language: 'fortran',)
  add_global_arguments(['-D UPPERCASE_FORTRAN','-D NO_APPEND_FORTRAN', ], language: 'c',)
else
  add_global_arguments(['-cpp', ], language: 'fortran',)
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
message('pypath', py_installation.path())
py_exe = find_program('python', py_installation.path() / 'bin')
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
  dependencies: [lib_arpack, lib_umfpack, lib_blas, libs_extra, dep_py],
)
