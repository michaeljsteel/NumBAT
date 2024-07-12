project('mes_test', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

# These arguments are only used to build the shared library
# not the executables that use the library.
lib_args = ['-DBUILDING_MES_TEST']

shlib = shared_library('mes_test', 'mes_test.f90',
  install : true,
  fortran_args : lib_args,
  gnu_symbol_visibility : 'hidden',
)

test_exe = executable('mes_test_test', 'mes_test_test.f90',
  link_with : shlib)
test('mes_test', test_exe)

# Make this library usable as a Meson subproject.
mes_test_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : shlib)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name : 'mes_test',
  filebase : 'mes_test',
  description : 'Meson sample project.',
  subdirs : 'mes_test',
  libraries : shlib,
  version : '0.1',
)
