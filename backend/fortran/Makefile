# Makefile for NumBAT. May require some adjustment to reflect your installation.

NumBAT_MOD=NumBAT

all: $(NumBAT_MOD).so

## Choose your platform

MYPLAT="linux_inteloneapi"
#MYPLAT="linux_vanilla"
#MYPLAT="MacOS"
#MYPLAT="myplatform"

PLAT = $(MYPLAT)

ifeq ($(origin NUMBAT_PLATFORM), environment)
PLAT=$(NUMBAT_PLATFORM)
endif

ifndef PLAT
$(error 'Error: Must define MYPLAT variable in Makefile or environment variable NUMBAT_PLATFORM')
endif

UNAME_S := $(shell uname -s)

ifeq ($(PLAT), "linux_vanilla")

ifneq ($(UNAME_S), Linux)
$(error 'We do not seem to be building on Linux. Check the value of MYPLAT in backend/fortran/Makefile')
endif

F_COMPILER = g95
F_COMP_EXE = f95
F2PY = f2py3
C_COMP_EXE = gcc

# BLAS and LAPACK from generic GNU-Linux
#
ARPACK_INCL = -I/usr/local/include
ARPACK_LIB = -L/usr/local/lib -larpack

UMFPACK_INC = -I/usr/include/suitesparse
UMFPACK_LIB = -L/usr/local/lib -lumfpack -lcholmod -lamd

MATHINC = $(ARPACK_INCL)
MATHLIB = $(ARPACK_LIB) $(UMFPACK)_LIB -llapack -lblas
GENLIB = -lrt

PYVERMAJ=3                                         ## Check me
PYVERMAJMIN=3.11                                   ## Check me
PYROOT = /usr/lib/python$(PYVERMAJ)                ## Check me
PYINCCORE = -I /usr/include/python$(PYVERMAJMIN)   ## Check me

NUMPYDIR = $(PYROOT)/dist-packages/numpy
PYINC = -I $(PYINCCORE) -I $(NUMPYDIR)/f2py/src/


else ifeq ($(PLAT), "linux_inteloneapi")
# Uses BLAS and LAPACK from free Intel OneApi MKL library

ifneq ($(UNAME_S), Linux)
$(error 'We do not seem to be building on Linux. Check the value of MYPLAT in backend/fortran/Makefile')
endif

F_COMPILER = g95
F_COMP_EXE = f95
F2PY = f2py3
C_COMP_EXE = gcc

## Check from here until EndDefs
ARPACK_INCL = -I /usr/local/include
ARPACK_LIB  = -L /usr/local/lib -larpack

MKLROOT = /opt/intel/oneapi/mkl/latest

UMFPACK_INC = -I/usr/include/suitesparse
UMFPACK_LIB = -L/usr/local/lib -lumfpack -lcholmod -lamd

PYVERMAJ=3
PYVERMAJMIN=3.11
PYROOT = /usr/lib/python$(PYVERMAJ)
PYINCCORE = /usr/include/python$(PYVERMAJMIN)
## EndDefs



MKL_INCL = -I $(MKLROOT)/include
MKL_LIB  = -L $(MKLROOT)/lib/intel64 -lmkl_rt

MATHINC = $(ARPACK_INCL) $(MKL_INCL)
MATHLIB = $(ARPACK_LIB) $(UMFPACK_LIB) $(MKL_LIB)
GENLIB = -lrt

NUMPYDIR = $(PYROOT)/dist-packages/numpy
PYINC = -I $(PYINCCORE) -I $(NUMPYDIR)/f2py/src/


else ifeq ($(PLAT), "MacOS")

ifneq ($(UNAME_S), Darwin)
$(error 'We do not seem to be building on Mac OS. Check the value of MYPLAT in backend/fortran/Makefile')
endif

GCCROOT = /opt/local
GCCBIN = $(GCCROOT)/bin
GCCVER = mp-13
F_COMPILER = $(GCCBIN)/gfortran-$(GCCVER)
F_COMP_EXE = $(GCCBIN)/gfortran-$(GCCVER)
C_COMP_EXE = $(GCCBIN)/gcc-$(GCCVER)

# standard lapack, blas and arpack
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

#DEFS = -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

PYVER = 3.12
PYNAME = python$(PYVER)
PYROOT = /Library/Frameworks/Python.framework/Versions/$(PYVER)
NUMPYDIR = $(PYROOT)/lib/$(PYNAME)/site-packages/numpy
PYINC = -I $(PYROOT)/include/$(PYNAME) -I $(NUMPYDIR)/core/include -I $(NUMPYDIR)/f2py/src/
PYLIB = -L $(PYROOT)/lib -l$(PYNAME)
OSF77FLAGS = -D __APPLE__
F2PY = $(PYROOT)/bin/f2py

else ifeq ($(PLAT), "myplatform")
# You complete me!
MATHINC =
MATHLIB =
GENLIB = -lrt

endif


SRCS = $(wildcard src/*.f)
OBJSA = $(SRCS:.f=.o)
OBJS = $(subst src,build,$(OBJSA))

PYINTERFACE_SRCSA = conv_gmsh.f py_calc_modes_AC.f py_calc_modes.f \
	EM_mode_energy_int_v2_Ez.f \
  EM_mode_energy_int_Ez.f photoelastic_int.f photoelastic_int_v2.f \
  AC_mode_power_int.f AC_mode_power_int_v2.f AC_mode_power_int_v4.f \
  AC_alpha_int.f AC_alpha_int_v2.f moving_boundary.f array_material_AC.f \
  AC_mode_elastic_energy_int.f AC_mode_elastic_energy_int_v4.f EM_mode_E_energy_int.f \
  H_mode_field_Ez.f array_size.f

PYINTERFACE_SRCS = $(addprefix src/, $(PYINTERFACE_SRCSA))

#SuiteSparse details
SS_WRAP = lib
UMFPACK_OBJ = umf4_f77zwrapper64.o
UMFPACK_WRAPPER = $(SS_WRAP)/$(UMFPACK_OBJ)
#SSP_VER=3.7.1
#SSP_VER=4.2.0
SS_VER=7.1.0_linux
UMFPACK_SRC = $(SS_WRAP)/umf4_f77zwrapper_SS_$(SS_VER).c



#DEFS = -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION


#Flags come in a Distribution object passed through as attributes
#F77FLAGS=-Wall -Wno-surprising -g0 -fPIC -Ofast -funroll-loops -ffixed-form -fno-second-underscore -O3 -cpp
F77FLAGS=-Wall -Wno-surprising -g0 -fPIC -Ofast -funroll-loops -fimplicit-none -fno-second-underscore -O3 -cpp $(OSF77FLAGS)
CFLAGS=-Wall -g0 -fPIC -Ofast -funroll-loops

#F77FLAGS=-Wall  -fPIC -Ofast -funroll-loops
#CFLAGS=-Wall -Wno-format -g0 -Ofast -fwrapv

build/%.o : src/%.f
	$(F_COMP_EXE) -c $(F77FLAGS) $(MATHINCL) $< -o $@

build/%.o : src/%.c
	$(C_COMP_EXE) -c $(CFLAGS) $(PYINC) $< -o $@

build/%.o : build/%.f
	$(F_COMP_EXE) -c $(F77FLAGS) $(MATHINCL) $< -o $@

build/%.o : build/%.c
	$(C_COMP_EXE) -c $(CFLAGS) $(PYINC) $< -o $@

#NumBAT.so: NumBAT.pyf $(SRCS)  $(SS_LIB_LOCATION)/$(UMFPACK_NAME)
#	f2py3 -c NumBAT.pyf $(SRCS) --verbose --fcompiler=$(F_COMPILER) --f77exec=$(F_COMP_EXE) \
#	$(DEFS)  $(SS_LIB_LOCATION)/$(UMFPACK_NAME)   \
#	$(MATHINCL) $(MATHLIB)  $(UMFPACK_LIB)

F2PYLINK=Ccomp
#F2PYLINK=Fcomp

ifeq ($(F2PYLINK), Ccomp)
FMETHOD=null
CMETHOD=$(NumBAT_MOD)
else
FMETHOD=$(NumBAT_MOD)
CMETHOD=null
endif



check:
	echo $(F2PYLINK) $(FMETHOD) $(CMETHOD)

$(FMETHOD).so: $(FMETHOD).pyf $(SRCS)  $(UMFPACK_WRAPPER)
	echo "Old compilation method Step B"
	$(F2PY) -c NumBAT.pyf $(SRCS) --verbose --fcompiler=$(F_COMPILER) --f77exec=$(F_COMP_EXE) \
					--f77flags=$(F77FLAGS) --f90flags=$(F77FLAGS)  $(DEFS) \
			 	$(UMFPACK_WRAPPER) $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) $(GENLIB) --build-dir tbuild

$(FMETHOD).pyf: $(SRCS)
	echo "Old compilation method Step A"
	$(F2PY)  -h $(FMETHOD).pyf $(PYINTERFACE_SRCS) --overwrite-signature  --lower
	$(F2PY) -m $(FMETHOD)  $(PYINTERFACE_SRCS) --overwrite-signature  --lower  --build-dir build
	#mv $(NumBAT_MOD)-f2pywrappers.f $(NumBAT_MOD)module.c build

#$(FMETHOD).so: NumBAT.pyf $(OBJS)  $(UMFPACK_WRAPPER)
#	echo "Old compilation method"
#	f2py3 -c NumBAT.pyf $(OBJS) --verbose --fcompiler=$(F_COMPILER) --f77exec=$(F_COMP_EXE) \
#	$(DEFS)  $(UMFPACK_WRAPPER) \
##	$(MATHINCL) $(MATHLIB)  $(UMFPACK_LIB) --build-dir tbuild



F2PY_SRCS = build/$(NumBAT_MOD)module.c build/$(NumBAT_MOD)-f2pywrappers.f build/fortranobject.f
F2PY_OBJSA = $(F2PY_SRCS:.c=.o)
F2PY_OBJS = $(F2PY_OBJSA:.f=.o)

$(CMETHOD).so: $(CMETHOD).pyf $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)
	echo "New compilation method Step B"
	#gcc -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER) $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) -lgfortran -lm -o $(NumBAT_MOD).so
	#gfortran -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER) $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) $(GENLIB) -lgfortran -lm -o $(NumBAT_MOD).so
	$(C_COMP_EXE) -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER)  $(MATHLIB) $(GENLIB) $(PYLIB) -lgfortran -lm -o $(NumBAT_MOD).so

$(UMFPACK_WRAPPER): $(UMFPACK_SRC)
	$(C_COMP_EXE) -O5 -DZLONG $(UMFPACK_INC) -c $(UMFPACK_SRC) \
		-fPIC -o $(UMFPACK_WRAPPER)
test:
	echo $(NUMPYDIR)
$(CMETHOD).pyf: $(SRCS)
	echo "New compilation method Step A"
	-rm $(F2PY_SRCS)
	-mkdir build
	$(F2PY) -m $(NumBAT_MOD) -h $(NumBAT_MOD).pyf $(PYINTERFACE_SRCS) --overwrite-signature --lower
	$(F2PY) $(PYINTERFACE_SRCS) -m $(NumBAT_MOD) --build-dir build --lower
	cp $(NUMPYDIR)/f2py/src/fortranobject.c build
	#mv $(NumBAT_MOD)-f2pywrappers.f $(NumBAT_MOD)module.c build



clean:
	-rm -f $(SS_WRAP)/*.o
	-rm -f build/*.o build/*.f build/*.c
	-rm -f $(NumBAT_MOD).pyf $(NumBAT_MOD).so $(NumBAT_MOD).*.so
