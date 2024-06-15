# Makefile for NumBAT. May require some adjustment to reflect your installation.

NumBAT_MOD=NumBAT

all: $(NumBAT_MOD).so

## Set your python virtualenv platform
MYPYENV = $(HOME)/mypy3-12
PYVERMAJ=3
#PYVERMAJMIN=$(origin NUMBAT_PYVER),3.12)
PYVERMAJMIN=3.12
#PYVERMAJMIN=3.11

# virtual env doesn't contain Python.h, so we grab it from here:
SYSTEMPYINC = /usr/include/python$(PYVERMAJMIN)


## Set your other libraries

MYPLAT="linux_inteloneapi"
#MYPLAT="linux_vanilla"
#MYPLAT="MacOS"
#MYPLAT="myplatform"

COMPILER_VENDOR="gcc"
#COMPILER_VENDOR="intel"


ifeq ($(COMPILER_VENDOR), "gcc")
F_COMP_EXE = f95
C_COMP_EXE = gcc

F77FLAGS=-Wall -Wno-surprising -g0 -fPIC -Ofast -funroll-loops -fimplicit-none -fno-second-underscore -O3 -cpp $(OSF77FLAGS) -J ./build
F95FLAGS=$(F77FLAGS)
CFLAGS=-Wall -g0 -fPIC -Ofast -funroll-loops
COMP_LIB = -lgfortran -lm

else

#TODO: smooth selection of version folders
ONEAPI=/opt/intel/oneapi
ONEAPICOMP=/opt/intel/oneapi/compiler/latest/bin
COMP_LIB = -L $(ONEAPI)/2024.1/lib -limf

F_COMP_EXE = $(ONEAPICOMP)/ifx
C_COMP_EXE = $(ONEAPICOMP)/icx

F77FLAGS= -nologo -fpp -fPIC -O2 -warn all -warn noexternals -module mods
CFLAGS=  -fPIC -O2 -Wall
F95FLAGS=$(F77FLAGS)

endif





PLAT = $(MYPLAT)

# or set it from an environment variable
#ifeq ($(origin NUMBAT_PLATFORM), environment)
#PLAT=$(NUMBAT_PLATFORM)
#endif

ifndef PLAT
$(error 'Error: Must define MYPLAT variable in Makefile or environment variable NUMBAT_PLATFORM')
endif

UNAME_S := $(shell uname -s)



ifeq ($(PLAT), "linux_vanilla")   ##### LINUX VANILLA
                                  ##### Uses BLAS and LAPACK from generic GNU-Linux

ifneq ($(UNAME_S), Linux)  # TODO: make this a make function?
$(error 'We do not seem to be building on Linux. Check the value of MYPLAT in backend/fortran/Makefile')
endif



## Check from here until EndDefs
ARPACK_INC = -I /usr/local/include
ARPACK_LIB  = -L /usr/lib/x86_64-linux-gnu -larpack

UMFPACK_INC = -I/usr/include/suitesparse
UMFPACK_LIB = -L /usr/lib/x86_64-linux-gnu -lumfpack -lcholmod -lamd -lsuitesparseconfig
## EndDefs

PYROOT = $(MYPYENV)
PYBIN = $(PYROOT)/bin
PYBASEINC= $(PYROOT)/include/python$(PYVERMAJMIN)
PYPACKAGES = $(PYROOT)/lib/python$(PYVERMAJMIN)/site-packages
F2PY = $(PYBIN)/f2py

MATHINC = $(ARPACK_INC)
MATHLIB = $(ARPACK_LIB) $(UMFPACK_LIB) -llapack -lblas
GENLIB = -lrt

NUMPYDIR = $(PYPACKAGES)/numpy
PYINC = -I $(PYBASEINC) -I $(NUMPYDIR)/core/include/ -I $(NUMPYDIR)/f2py/src/ -I $(SYSTEMPYINC)


else ifeq ($(PLAT), "linux_inteloneapi")  ##### LINUX INTELONEAPI with gfortran
# Uses BLAS and LAPACK from free Intel OneApi MKL library

ifneq ($(UNAME_S), Linux)
$(error 'We do not seem to be building on Linux. Check the value of MYPLAT in backend/fortran/Makefile')
endif

## Check from here until EndDefs
ARPACK_INC = -I /usr/local/include
ARPACK_LIB  = -L /usr/local/lib -larpack

MKLROOT = /opt/intel/oneapi/mkl/latest

UMFPACK_INC = -I/usr/include/suitesparse
UMFPACK_LIB = -L/usr/local/lib -lumfpack -lcholmod -lamd
## EndDefs

PYROOT = $(MYPYENV)
PYBIN = $(PYROOT)/bin
PYBASEINC= $(PYROOT)/include/python$(PYVERMAJMIN)
PYPACKAGES = $(PYROOT)/lib/python$(PYVERMAJMIN)/site-packages

F2PY = $(PYBIN)/f2py

MKL_INC = -I $(MKLROOT)/include
MKL_LIB  = -L $(MKLROOT)/lib/intel64 -lmkl_rt

MATHINC = $(ARPACK_INC) $(MKL_INC)
MATHLIB = $(ARPACK_LIB) $(UMFPACK_LIB) $(MKL_LIB)
GENLIB = -lrt

NUMPYDIR = $(PYPACKAGES)/numpy
PYINC = -I $(PYBASEINC) -I $(NUMPYDIR)/core/include/ -I $(NUMPYDIR)/f2py/src/ -I $(SYSTEMPYINC)


else ifeq ($(PLAT), "MacOS")
# MacOS with standard LAPACK, BLAS and ARPACK

ifneq ($(UNAME_S), Darwin)
$(error 'We do not seem to be building on Mac OS. Check the value of MYPLAT in backend/fortran/Makefile')
endif

GCCROOT = /opt/local
GCCBIN = $(GCCROOT)/bin
GCCVER = mp-13

F_COMP_EXE = $(GCCBIN)/gfortran-$(GCCVER)
C_COMP_EXE = $(GCCBIN)/gcc-$(GCCVER)


LAPACK_LIB = -L /opt/local/lib -llapack -lblas
ARPACK_INC = -I/opt/local/include
ARPACK_LIB = -L/opt/local/lib -larpack $(LAPACK_LIB)
#MKLROOT =
#MKL_INC =
#MKL_LIB =
UMFPACK_INC = -I /opt/local/include
UMFPACK_LIB = -L/opt/local/lib/ -lumfpack -lcholmod -lamd
MATHINC = /opt/local/include $(MKL_INC) $(ARPACK_INC)
MATHLIB = $(ARPACK_LIB) $(UMFPACK_LIB) $(MKL_LIB) -L /opt/local/lib/gcc13
#GENLIB = -lrt

PYROOT = $(MYPYENV)
PYBIN = $(PYROOT)/bin
PYBASEINC= $(PYROOT)/include/python$(PYVERMAJMIN)
PYPACKAGES = $(PYROOT)/lib/python$(PYVERMAJMIN)/site-packages

PYNAME = python$(PYVERMAJMIN)
#NUMPYDIR = $(PYROOT)/lib/$(PYNAME)/site-packages/numpy

PYFRAMEWORK = /Library/Frameworks/Python.framework/Versions/Current
PYCINC = -I $(PYFRAMEWORK)/Headers

NUMPYDIR = $(PYPACKAGES)/numpy
PYINC = -I $(PYROOT)/include/$(PYNAME) -I $(NUMPYDIR)/core/include -I $(NUMPYDIR)/f2py/src/ $(PYCINC)/ -I /opt/local/include/gcc13/c++

PYLIB = -L $(PYROOT)/lib -l$(PYNAME)
PYLIB = -L $(PYFRAMEWORK)/lib -l$(PYNAME)

OSF77FLAGS = -D __APPLE__
F2PY = $(PYROOT)/bin/f2py

else ifeq ($(PLAT), "myplatform")
# You complete me!
MATHINC =
MATHLIB =
GENLIB = -lrt

else
$(error "Unknown platform: $(NUMBAT_PLATFORM)" )
endif


SRCS_F = $(wildcard src/*.f)
SRCS_F90 = $(wildcard src/*.f90)

#OBJSA = $(SRCS_F:.f=.o) $(SRCS_F90:.f90=.o)
#OBJS = $(subst src,build,$(OBJSA))
OBJS = $(subst src/,build/,$(SRCS_F:.f=.o) $(SRCS_F90:.f90=.o))

#The build order of these matters
SRCS_MODS = src/numbatmod_module.f90 src/nbinterfaces_module.f90

OBJS_MODA = $(SRCS_MODS:.f90=.o)
#OBJS_MOD= $(subst src,build,$(OBJS_MODA))
MODULES_MODS = $(subst src/,build/,$(SRCS_MODS:.f90=.mod))

tt:
	echo $(SRCS_F)
	echo $(SRCS_MODS)
	echo $(OBJS_MOD)
	echo $(MODULES_MODS)

PYINTERFACE_SRCS_IN = pyinterfaces.f90 conv_gmsh.f90 \
	em_mode_energy_int_v2_ez.f em_mode_energy_int_ez.f em_mode_e_energy_int.f \
	h_mode_field_ez.f \
	ac_mode_power_int.f ac_mode_power_int_v2.f ac_mode_power_int_v4.f \
	ac_mode_elastic_energy_int.	f ac_mode_elastic_energy_int_v4.f \
	ac_alpha_int.f ac_alpha_int_v2.f  \
	array_size.f array_material_ac.f \
	photoelastic_int.f photoelastic_int_v2.f moving_boundary.f \
	py_calc_modes.f90 py_calc_modes_ac.f90

PYINTERFACE_SRCS = $(addprefix src/, $(PYINTERFACE_SRCS_IN))

# Local SuiteSparse wrapper
SSWRAP_DIR = sswrap
SS_VER=7.1.0_linux
SSWRAP_OBJ = build/umf4_f77zwrapper64.o
SSWRAP_SRC = $(SSWRAP_DIR)/umf4_f77zwrapper_SS_$(SS_VER).c

$(SSWRAP_OBJ): $(SSWRAP_SRC)
	$(C_COMP_EXE) $(CFLAGS) -I $(SSWRAP_DIR) -D ZLONG $(UMFPACK_INC) -c $< -o $@


mods: $(MODULES_MODS)

build/%.mod : src/%.f90
	$(F_COMP_EXE) -c $(F95FLAGS) $(MATHINC) $< -o $@

# Any module file change will cause rebuilding of all .f files
build/%.o : src/%.f90
	$(F_COMP_EXE) -c $(F95FLAGS) $(MATHINC) $< -o $@

build/%.o : src/%.f $(MODULES_MODS)
	$(F_COMP_EXE) -c $(F77FLAGS) $(MATHINC) $< -o $@

build/%.o : src/%.c
	$(C_COMP_EXE) -c $(CFLAGS) $(PYINC) $< -o $@

build/%.o : build/%.f
	$(F_COMP_EXE) -c $(F77FLAGS) $(MATHINC) $< -o $@

build/%.o : build/%.c
	$(C_COMP_EXE) -c $(CFLAGS) $(PYINC) $< -o $@

build/%.e : src/%.f
	cpp -I src $< -o $@

build/%.e : src/%.c
	$(C_COMP_EXE) -E $(CFLAGS) $(PYINC) $< -o $@

build/%.e : build/%.c
	$(C_COMP_EXE) -E $(CFLAGS) $(PYINC) $< -o $@

F2PY_SRCS = build/$(NumBAT_MOD)module.c build/$(NumBAT_MOD)-f2pywrappers.f build/fortranobject.f
F2PY_OBJSA = $(F2PY_SRCS:.c=.o)
F2PY_OBJS = $(F2PY_OBJSA:.f=.o)

$(NumBAT_MOD).pyf: $(SRCS)
	-rm $(F2PY_SRCS)
	-mkdir build
	$(F2PY) -m $(NumBAT_MOD) -h $(NumBAT_MOD).pyf $(PYINTERFACE_SRCS) --overwrite-signature --lower
	$(F2PY) $(PYINTERFACE_SRCS) -m $(NumBAT_MOD) --build-dir build --lower -include\"../src/myhead.h\"
	cp $(NUMPYDIR)/f2py/src/fortranobject.c build


$(NumBAT_MOD).so: $(NumBAT_MOD).pyf $(MODULES_MODS) $(OBJS) $(F2PY_OBJS) $(SSWRAP_OBJ)
	$(C_COMP_EXE) -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(SSWRAP_OBJ) \
	 $(MATHLIB) $(GENLIB) $(PYLIB) $(COMP_LIB) -lm -o $(NumBAT_MOD).so


clean:
	-rm -f $(SS_WRAP)/*.o
	-rm -f build/*.o build/*.f build/*.c build/*.f90
	-rm ./*.mod src/*.mod build/*.mod
	-rm mods/*.mod mods/*genmod.mod mods/*genmod.f90
	-rm -f $(NumBAT_MOD).pyf $(NumBAT_MOD).so $(NumBAT_MOD).*.so


check:
	echo objs: $(OBJS)
	echo f2pyobjs: $(F2PY_OBJS)
	echo mods: $(SRCS_MODS)
	echo mods: $(MODULES_MODS)




