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

inc_nblocal = ''
if host_machine.system() == 'linux'
  dir_mkl = dir_oneapi / 'mkl/latest/lib'
  dir_nblocal = ['/usr/local/lib']

  inc_nblocal = include_directories('/usr/local/include')
  lib_blas = ccomp.find_library('mkl_rt',
    dirs: dir_mkl,
    header_include_directories: inc_nblocal)

elif host_machine.system() == 'darwin'
  lib_blas = [
    ccomp.find_library('lapack', dirs: '/opt/local/lib'),
    ccomp.find_library('blas', dirs: '/opt/local/lib'),
  ]

else

  #MKL_LIB=$(MKLROOT)\lib\mkl_intel_lp64.lib $(MKLROOT)\lib\mkl_rt.lib \
  #   $(MKLROOT)\lib\mkl_core_dll.lib  $(MKLROOT)\lib\mkl_intel_thread_dll.lib

  dir_nblocal = 'C:/Users/mq98002769/numbat/usr_local/lib'
  dir_mkl = 'c:/Program Files (x86)/Intel/oneAPI/mkl/2024.1/lib'
  inc_nblocal = include_directories('C:/Users/mq98002769/numbat/usr_local/include')

  lib_blas = ccomp.find_library('mkl_rt', dirs: dir_mkl,
  header_include_directories:inc_nblocal)
endif

if ccomp.get_id() == 'intel-llvm'
  lib_oneapi = ccomp.find_library('ifport',
                                  dirs: dir_mkl
                                  )
  libs_extra = [lib_oneapi]
endif

lib_arpack = ccomp.find_library('arpack', dirs: dir_nblocal)
lib_umfpack = ccomp.find_library('umfpack', dirs: dir_nblocal)

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

#add_global_arguments(['-DNO_APPEND_FORTRAN',  '-DNO_UPPERCASE_FORTRAN'], language: 'c')

comp_f2py = find_program('f2py')

subdir('sswrap')
subdir('src')

message(
  '\n\n\n\nmeson settings',
  meson.current_build_dir(),
  meson.build_root(),
  meson.current_source_dir(),
  meson.source_root(),
)

pymod = import('python')
py_installation = pymod.find_installation()
py_exe = find_program('python', py_installation.path() / 'bin')
dep_py = py_installation.dependency()

#message('Using python installation at: ', py_installation.path())

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

# TODO: unneeded
incdir = include_directories(
  [],
)
incs_f2py = include_directories([dir_numpy, dir_f2py])

message([dir_numpy, dir_f2py])

incs = [inc_nblocal, incdir, incs_f2py]

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
