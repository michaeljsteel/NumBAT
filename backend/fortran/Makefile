# Makefile for NumBAT. May require some adjustment to reflect your installation.

NumBAT_MOD=NumBAT

all: $(NumBAT_MOD).so

## Set your python virtualenv platform
MYPYENV = $(HOME)/mypy3-12
PYVERMAJ=3
#PYVERMAJMIN=$(origin NUMBAT_PYVER),3.12)
PYVERMAJMIN=3.12
#PYVERMAJMIN=3.11


## Set your other libraries

MYPLAT="linux_inteloneapi"
#MYPLAT="linux_intel_ifx"
#MYPLAT="linux_vanilla"
#MYPLAT="MacOS"
#MYPLAT="myplatform"

COMPILER_VENDOR="gcc"
#COMPILER_VENDOR="intel"


ifeq ($(COMPILER_VENDOR), "gcc")
#F_COMP_EXE = g95
F_COMP_EXE = f95
C_COMP_EXE = gcc

F77FLAGS=-Wall -Wno-surprising -g0 -fPIC -Ofast -funroll-loops -fimplicit-none -fno-second-underscore -O3 -cpp $(OSF77FLAGS) -J ./build
F95FLAGS=$(F77FLAGS)
CFLAGS=-Wall -g0 -fPIC -Ofast -funroll-loops
COMP_LIB = -l gfortran

else

INTELONEAPIDIR=/opt/intel/oneapi
INTELCOMPDIR=/opt/intel/oneapi/compiler/latest/bin
#F_COMP_EXE = $(INTELCOMPDIR)/ifx
F_COMP_EXE = $(INTELCOMPDIR)/ifx
C_COMP_EXE = $(INTELCOMPDIR)/icx

F77FLAGS= -nologo -fpp -fPIC -O2 -warn all -warn noexternals -module mods
CFLAGS=  -O2 -Wall -fPIC
F95FLAGS=$(F77FLAGS)
COMP_LIB = -L $(INTELONEAPIDIR)/2024.1/lib -limf


endif



# virtual env doesn't contain Python.h, so we grab it from here:
SYSTEMPYINC = /usr/include/python$(PYVERMAJMIN)


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
ARPACK_INCL = -I /usr/local/include
ARPACK_LIB  = -L /usr/lib/x86_64-linux-gnu -larpack

#MKLROOT = /opt/intel/oneapi/mkl/latest

UMFPACK_INC = -I/usr/include/suitesparse
UMFPACK_LIB = -L /usr/lib/x86_64-linux-gnu -lumfpack -lcholmod -lamd -lsuitesparseconfig
## EndDefs

PYROOT = $(MYPYENV)
PYBIN = $(PYROOT)/bin
PYBASEINC= $(PYROOT)/include/python$(PYVERMAJMIN)
PYPACKAGES = $(PYROOT)/lib/python$(PYVERMAJMIN)/site-packages
F2PY = $(PYBIN)/f2py




MATHINC = $(ARPACK_INCL)
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
ARPACK_INCL = -I /usr/local/include
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

MKL_INCL = -I $(MKLROOT)/include
MKL_LIB  = -L $(MKLROOT)/lib/intel64 -lmkl_rt

MATHINC = $(ARPACK_INCL) $(MKL_INCL)
MATHLIB = $(ARPACK_LIB) $(UMFPACK_LIB) $(MKL_LIB)
GENLIB = -lrt

NUMPYDIR = $(PYPACKAGES)/numpy
PYINC = -I $(PYBASEINC) -I $(NUMPYDIR)/core/include/ -I $(NUMPYDIR)/f2py/src/ -I $(SYSTEMPYINC)

else ifeq ($(PLAT), "MacOS")

else ifeq ($(PLAT), "MacOS")
# MacOS with standard LAPACK, BLAS and ARPACK

ifneq ($(UNAME_S), Darwin)
$(error 'We do not seem to be building on Mac OS. Check the value of MYPLAT in backend/fortran/Makefile')
endif

GCCROOT = /opt/local
GCCBIN = $(GCCROOT)/bin
GCCVER = mp-13
#F_COMP_EXE = $(GCCBIN)/gfortran-$(GCCVER)
F_COMP_EXE = $(GCCBIN)/gfortran-$(GCCVER)
C_COMP_EXE = $(GCCBIN)/gcc-$(GCCVER)


LAPACK_LIB = -L /opt/local/lib -llapack -lblas
ARPACK_INCL = -I/opt/local/include
ARPACK_LIB = -L/opt/local/lib -larpack $(LAPACK_LIB)
#MKLROOT =
#MKL_INCL =
#MKL_LIB =
UMFPACK_INC = -I /opt/local/include
UMFPACK_LIB = -L/opt/local/lib/ -lumfpack -lcholmod -lamd
MATHINC = /opt/local/include $(MKL_INCL) $(ARPACK_INCL)
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


SRCS_MODS = src/numbatmod_module.f90 src/nbinterfaces_module.f90
SRCS_F = $(wildcard src/*.f)
SRCS_F95 = $(wildcard src/*.f90)

OBJSA = $(SRCS_F:.f=.o) $(SRCS_F95:.f90=.o)
OBJS = $(subst src,build,$(OBJSA))

OBJS_MODA = $(SRCS_MODS:.f90=.o)
#OBJS_MOD= $(subst src,build,$(OBJS_MODA))
MODULES_MODS = build/numbatmod_module.mod build/nbinterfaces_module.mod

tt:
	echo $(SRCS_F)
	echo $(SRCS_MODS)
	echo $(OBJS_MOD)
	echo $(MODULES_MODS)

PYINTERFACE_SRCSA = pyinterfaces.f90 py_calc_modes_ac.f90 py_calc_modes.f90 \
	em_mode_energy_int_v2_ez.f conv_gmsh.f90 \
  em_mode_energy_int_ez.f photoelastic_int.f photoelastic_int_v2.f \
  ac_mode_power_int.f ac_mode_power_int_v2.f ac_mode_power_int_v4.f \
  ac_alpha_int.f ac_alpha_int_v2.f moving_boundary.f array_material_ac.f \
  ac_mode_elastic_energy_int.f ac_mode_elastic_energy_int_v4.f em_mode_e_energy_int.f \
  h_mode_field_ez.f array_size.f

PYINTERFACE_SRCS = $(addprefix src/, $(PYINTERFACE_SRCSA))

#SuiteSparse details
SS_WRAP = lib
UMFPACK_OBJ = umf4_f77zwrapper64.o
UMFPACK_WRAPPER = $(SS_WRAP)/$(UMFPACK_OBJ)
#SSP_VER=3.7.1
#SSP_VER=4.2.0
SS_VER=7.1.0_linux
UMFPACK_SRC = $(SS_WRAP)/umf4_f77zwrapper_SS_$(SS_VER).c




mods: $(MODULES_MODS)

build/%.mod : src/%.f90
	$(F_COMP_EXE) -c $(F95FLAGS) $(MATHINCL) $< -o $@

# Any module file change will cause rebuilding of all .f files
build/%.o : src/%.f90
	$(F_COMP_EXE) -c $(F95FLAGS) $(MATHINCL) $< -o $@

build/%.o : src/%.f $(MODULES_MODS)
	$(F_COMP_EXE) -c $(F77FLAGS) $(MATHINCL) $< -o $@

build/%.o : src/%.c
	$(C_COMP_EXE) -c $(CFLAGS) $(PYINC) $< -o $@

build/%.o : build/%.f
	$(F_COMP_EXE) -c $(F77FLAGS) $(MATHINCL) $< -o $@

build/%.o : build/%.c
	$(C_COMP_EXE) -c $(CFLAGS) $(PYINC) $< -o $@

build/%.e : src/%.f
	cpp -I src $< -o $@

build/%.e : src/%.c
	$(C_COMP_EXE) -E $(CFLAGS) $(PYINC) $< -o $@

build/%.e : build/%.c
	$(C_COMP_EXE) -E $(CFLAGS) $(PYINC) $< -o $@

F2PYLINK=Ccomp
#F2PYLINK=Fcomp

ifeq ($(F2PYLINK), Ccomp)
FMETHOD=null
CMETHOD=$(NumBAT_MOD)
else
FMETHOD=$(NumBAT_MOD)
CMETHOD=null
endif


F2PY_SRCS = build/$(NumBAT_MOD)module.c build/$(NumBAT_MOD)-f2pywrappers.f build/fortranobject.f
#build/$(NumBAT_MOD)-f2pywrappers2.f90
F2PY_OBJSA = $(F2PY_SRCS:.c=.o)
F2PY_OBJS = $(F2PY_OBJSA:.f=.o)

$(CMETHOD).pyf: $(SRCS)
	-rm $(F2PY_SRCS)
	-mkdir build
	$(F2PY) -m $(NumBAT_MOD) -h $(NumBAT_MOD).pyf $(PYINTERFACE_SRCS) --overwrite-signature --lower
	$(F2PY) $(PYINTERFACE_SRCS) -m $(NumBAT_MOD) --build-dir build --lower -include\"../src/myhead.h\"
	cp $(NUMPYDIR)/f2py/src/fortranobject.c build


$(CMETHOD).so: $(CMETHOD).pyf $(MODULES_MODS) $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)
	$(C_COMP_EXE) -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)  $(MATHLIB) $(GENLIB) $(PYLIB) $(COMP_LIB) -lm -o $(NumBAT_MOD).so

$(UMFPACK_WRAPPER): $(UMFPACK_SRC)
	$(C_COMP_EXE) -O5 -DZLONG $(UMFPACK_INC) -c $(UMFPACK_SRC) \
		-fPIC -o $(UMFPACK_WRAPPER)



$(CMETHOD).$(SHLIB): $(CMETHOD).pyf $(MODULES_MODS) $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)
	echo "New compilation method Step B"
	$(C_COMP_EXE) -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)  $(MATHLIB) $(GENLIB) $(PYLIB) -lgfortran -lm -o $(NumBAT_MOD).$(SHLIB)


# gcc -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)
# $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) -lgfortran -lm -o $(NumBAT_MOD).$(SHLIB)
# gfortran -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)
# $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) $(GENLIB) -lgfortran -lm -o $(NumBAT_MOD).$(SHLIB)






clean:
	-rm -f $(SS_WRAP)/*.o
	-rm -f build/*.o build/*.f build/*.c build/*.f90
	-rm ./*.mod src/*.mod build/*.mod
	-rm mods/*.mod mods/*genmod.mod mods/*genmod.f90
	-rm -f $(NumBAT_MOD).pyf $(NumBAT_MOD).so $(NumBAT_MOD).*.so


	#cp /Library/Frameworks/Python.framework/Headers/Python.h src
	#mv $(NumBAT_MOD)-f2pywrappers.f $(NumBAT_MOD)module.c build


#gcc -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER) $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) -lgfortran -lm -o $(NumBAT_MOD).so
#gfortran -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER) $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) $(GENLIB) -lgfortran -lm -o $
#TRASH SECTION


#NumBAT.so: NumBAT.pyf $(SRCS)  $(SS_LIB_LOCATION)/$(UMFPACK_NAME)
#	f2py3 -c NumBAT.pyf $(SRCS) --verbose --fcompiler=$(F_COMP_EXE) --f77exec=$(F_COMP_EXE) \
#	$(DEFS)  $(SS_LIB_LOCATION)/$(UMFPACK_NAME)   \
#	$(MATHINCL) $(MATHLIB)  $(UMFPACK_LIB)


#$(FMETHOD).so: NumBAT.pyf $(OBJS)  $(UMFPACK_WRAPPER)
#	echo "Old compilation method"
#	f2py3 -c NumBAT.pyf $(OBJS) --verbose --fcompiler=$(F_COMP_EXE) --f77exec=$(F_COMP_EXE) \
#	$(DEFS)  $(UMFPACK_WRAPPER) \
##	$(MATHINCL) $(MATHLIB)  $(UMFPACK_LIB) --build-dir tbuild


#F77FLAGS=-Wall  -fPIC -Ofast -funroll-loops
#CFLAGS=-Wall -Wno-format -g0 -Ofast -fwrapv


#DEFS = -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
