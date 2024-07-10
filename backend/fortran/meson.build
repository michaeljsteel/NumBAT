
project('NumBAT Fortran backend', ['fortran','c'],
  version: '2.0.1',)

  #fortran_args : lib_args,
  #gnu_symbol_visibility : 'hidden',

srcs = files('src/csr_make_arrays.f90')

incdir = include_directories(
  [
'/usr/local/include',
'/usr/include',
])

nblib = library('NumBAT', srcs,
  incdir)

# Set build configuration based on environment variables
build_target(
  'gcc',
  compiler : 'gcc',
  cflags : [
    '-Wall',
    '-W',
    '-g',
    '-O2',
    '-fPIC',
    '-Ofast',
    '-funroll-loops'
  ],
  fortran_flags : [
    '-Wall',
    '-Wno-surprising',
    '-W',
    '-g',
    '-O2',
    '-fPIC',
    '-Ofast',
    '-funroll-loops',
    '-fimplicit-none',
    '-fno-second_underscore',
    '-O3',
    '-cpp',
    '@OSF77FLAGS@',
    '-J ./mods',
    '-fall-intrinsics',
    '-pedantic'
  ],
  links : '-lgfortran -lm',
  sources : [
    'src/*.c',
    'src/*.f',
    'src/*.f90'
  ],
  depends : [ 'nb_fortran_module' ]  # Dependency on custom target
)

build_target(
  'intel',
  compiler : 'intel',
  cflags : [
    '-fPIC',
    '-O2',
    '-Wall',
    '-g'
  ],
  fortran_flags : [
    '-nologo',
    '-fpp',
    '-fPIC',
    '-O2',
    '-g',
    '-warn all',
    '-warn noexternals',
    '-module mods',
    '-I mods'
  ],
  links : [
    '-L@ONEAPI@/2024.2/lib',
    '-lifcoremt',
    '-lifport',
    '-limf'
  ],
  include_directories : '@ONEAPI@/include',  # Using environment variable
  sources : [
    'src/*.c',
    'src/*.f',
    'src/*.f90'
  ],
  depends : [ 'nb_fortran_module' ]  # Dependency on custom target
)

# Custom target for Fortran module compilation
custom_target(
  'nb_fortran_module',
  input : 'src/nb_fortranmodule.f90',
  output : 'build/nb_fortranmodule.mod',
  command : [
    executable('fcomp'),  # Assuming a custom Fortran compiler wrapper
    '-c',
    '@F95FLAGS@',
    '@MATHINC@',
    '-o',
    '@OUTPUT@',
    '@INPUT@'
  ]
)

# Define executables or libraries based on your needs
# ... (e.g., using `executable` or `shared_library` functions)

# Define test targets if applicable
# ... (e.g., using `test` function)

