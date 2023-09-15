# Makefile for NumBAT. May require some adjustment to reflect your installation.

#F_COMPILER = g95
F_COMPILER = gnu95
#F_COMP_EXE = f77
F_COMP_EXE = f95

#FORT_EXE = intelem
C_COMP_EXE = gcc

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
SS_VER=7.1.0
UMFPACK_SRC = $(SS_WRAP)/umf4_f77zwrapper_SS_$(SS_VER).c

MATH_PLAT="inteloneapi"
#MATH_PLAT="vanilla"
#MATH_PLAT="myplatform"

ifeq ($(MATH_PLAT), "vanilla")
# BLAS and LAPACK from generic GNU-Linux 
#
ARPACK_INCL = -I/usr/local/include
ARPACK_LIB = -L/usr/local/lib -larpack
UMFPACK_LIB = -L/usr/local/lib -lumfpack -lcholmod -lamd
MATHINC = $(ARPACK_INCL) 
MATHLIB = $(ARPACK_LIB) $(UMFPACK)_LIB) -llapack -lblas 
GENLIB = -lrt

else ifeq ($(MATH_PLAT), "inteloneapi")
# BLAS and LAPACK from free Intel OneApi MKL library
#
ARPACK_INCL = -I/usr/local/include
ARPACK_LIB = -L/usr/local/lib -larpack
MKLROOT = /opt/intel/oneapi/mkl/latest
MKL_INCL = -I$(MKLROOT)/include
MKL_LIB = -L$(MKLROOT)/lib/intel64 -lmkl_rt
UMFPACK_LIB = -L/usr/local/lib -lumfpack -lcholmod -lamd
MATHINC = $(MKL_INCL) $(ARPACK_INCL) 
MATHLIB = $(ARPACK_LIB) $(UMFPACK_LIB) $(MKL_LIB)
GENLIB = -lrt

else ifeq ($(MATH_PLAT), "myplatform")
# You complete me!
MATHINC = 
MATHLIB = 
GENLIB = -lrt

endif

#DEFS = -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
PYINC = -I /usr/include/python3.11 -I /usr/lib/python3/dist-packages/numpy/f2py/src/


#Flags come in a Distribution object passed through as attributes
F77FLAGS=-Wall -Wno-surprising -g0 -fPIC -Ofast -funroll-loops -ffixed-form -fno-second-underscore -O3
CFLAGS=-Wall -g0 -fPIC -Ofast -funroll-loops -ffixed-form 

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

NumBAT_MOD=NumBAT

F2PYLINK=Ccomp
#F2PYLINK=Fcomp

ifeq ($(F2PYLINK), Ccomp)
FMETHOD=null
CMETHOD=$(NumBAT_MOD)
else
FMETHOD=$(NumBAT_MOD)
CMETHOD=null
endif

all: $(NumBAT_MOD).so

check: 
	echo $(F2PYLINK) $(FMETHOD) $(CMETHOD)

$(FMETHOD).so: $(FMETHOD).pyf $(SRCS)  $(UMFPACK_WRAPPER)
	echo "Old compilation method Step B"
	f2py3 -c NumBAT.pyf $(SRCS) --verbose --fcompiler=$(F_COMPILER) --f77exec=$(F_COMP_EXE) \
					--f77flags=$(F77FLAGS) --f90flags=$(F77FLAGS)  $(DEFS) \
			 	$(UMFPACK_WRAPPER) $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) $(GENLIB) --build-dir tbuild 

$(FMETHOD).pyf: $(SRCS)
	echo "Old compilation method Step A"
	f2py3  -h $(FMETHOD).pyf $(PYINTERFACE_SRCS) --overwrite-signature  --lower 
	f2py3 -m $(FMETHOD)  $(PYINTERFACE_SRCS) --overwrite-signature  --lower  --build-dir build
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
	gcc -shared -Wall --verbose $(DEFS)  $(OBJS) $(F2PY_OBJS) $(UMFPACK_WRAPPER) $(UMFPACK_LIB) $(MATHINCL) $(MATHLIB) $(GENLIB) -lgfortran -lm -o $(NumBAT_MOD).so

$(UMFPACK_WRAPPER): $(UMFPACK_SRC)
	$(C_COMP_EXE) -O5 -DZLONG -I/usr/include/suitesparse -c $(UMFPACK_SRC) \
		-fPIC -o $(UMFPACK_WRAPPER)

$(CMETHOD).pyf: $(SRCS)
	echo "New compilation method Step A"
	-rm $(F2PY_SRCS)
	-mkdir build
	f2py3 -m $(NumBAT_MOD) -h $(NumBAT_MOD).pyf $(PYINTERFACE_SRCS) --overwrite-signature --lower 
	f2py3 $(PYINTERFACE_SRCS) -m $(NumBAT_MOD) --build-dir build --lower
	cp /usr/lib/python3/dist-packages/numpy/f2py/src/fortranobject.c build
	#mv $(NumBAT_MOD)-f2pywrappers.f $(NumBAT_MOD)module.c build


#purge:
#	-rm $(NumBAT_MOD).so $(NumBAT_MOD).*.so

clean:
	-rm -f $(SS_WRAP)/*.o
	-rm -rf tbuild/*
	-rm -f build/*.o build/*.f build/*.c
	-rm -f $(NumBAT_MOD).pyf $(NumBAT_MOD).so $(NumBAT_MOD).*.so
