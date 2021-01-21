
PYTHON=python3

all: tut01 tut02 tut03a tut03b tut04 tut05 tut06 tut06 tut07 tut08 tut09 tut10 tut11a 
allfast: tut01 tut02 tut03a tut03b tut04 tut05 ftut06 ftut06 ftut07 ftut08 ftut09 ftut10 ftut11a 


tut01: 
	$(PYTHON) simo-tut_01-first_calc.py 
  
tut02: 
	$(PYTHON) simo-tut_02-gain_spectra-npsave.py 

tut03a: 
	$(PYTHON) simo-tut_03_1-dispersion-npload.py 

tut03b: 
	$(PYTHON) simo-tut_03_2-dispersion-multicore.py 

tut04: 
	$(PYTHON) simo-tut_04-scan_widths.py 

tut05: 
	$(PYTHON) simo-tut_05-convergence_study.py 

tut06: 
	$(PYTHON) simo-tut_06-silica_nanowire.py 
ftut06: 
	$(PYTHON) simo-tut_06-silica_nanowire.py fast=1

tut07: 
	$(PYTHON) simo-tut_07-slot.py 
ftut07: 
	$(PYTHON) simo-tut_07-slot.py fast=1

tut08: 
	$(PYTHON) simo-tut_08-slot_coated-scan.py 
ftut08: 
	$(PYTHON) simo-tut_08-slot_coated-scan.py fast=1

tut09: 
	$(PYTHON) simo-tut_09-anisotropy.py
ftut09: 
	$(PYTHON) simo-tut_09-anisotropy.py fast=1

tut10: 
	$(PYTHON) simo-tut_10-onion.py
ftut10: 
	$(PYTHON) simo-tut_10-onion.py fast=1

tut11a: 
	$(PYTHON) simo-tut_11a-onion2.py
ftut11a: 
	$(PYTHON) simo-tut_11a-onion2.py fast=1

tut11b: 
	$(PYTHON) simo-tut_11b-onion3.py
ftut11b: 
	$(PYTHON) simo-tut_11b-onion3.py fast=1

tut12: 
	$(PYTHON) simo-tut_12-smf28.py
ftut12: 
	$(PYTHON) simo-tut_12-smf28.py fast=1

clean:
	rm *.png *.npz
	rm -rf tut*fields ftut*fields
