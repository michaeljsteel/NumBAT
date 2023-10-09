# Makefile for NumBAT 
#


build:
	export FFLAGS=-fallow-argument-mismatch  # handle stricter gfortran rules in GCC10, comment out for earlier
	cd backend/fortran && make 

installdeps:
	apt-get install -y python3-numpy python3-dev python3-scipy python3-nose python3-pip gfortran make gmsh libatlas-base-dev libblas-dev liblapack-dev libsuitesparse-dev
	pip3 install matplotlib 

examples:
	cd tutorials && make 
	cd lit_examples && make 
	cd josab_tutorial && make 

fastexamples:
	cd tutorials && make allfast
	cd lit_examples && make allfast
	cd josab_tutorial && make allfast

testinstall:
	cd tests && nosetests3

manual:
	-make -C docs latexpdf
	mv docs/build/latex/NumBAT.pdf .

clean:
		make -C backend/fortran clean
		make -C tutorials clean
		make -C lit_examples clean
		make -C JOSAB_tutorial clean
		make -C docs clean

cleanmesh:
	-rm backend/fortran/msh/build/*.msh
	-rm backend/fortran/msh/build/*.geo
	-rm backend/fortran/msh/build/*.mail
