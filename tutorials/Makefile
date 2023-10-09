
PYTHON=python3

all: tut01 tut02 tut03a tut03b tut04 tut05 tut06 tut07 \
				tut08 tut09 tut11a  tut11b tut12 tut13 tut14

allfast: tut01 tut02 tut03a ftut03b ftut04 ftut05 ftut06 ftut07 \
				ftut08 ftut09 ftut11a  ftut12 tut13 tut14

broken: tut10

slow: tut04 tut05 tut06 tut07 \
				tut08 tut09 tut11a  tut11b tut12 tut13 tut14
	

#fast1: tut01 ftut02 tut03a ftut03b ftut04 ftut05 ftut06 ftut07 
fast1: ftut02 ftut03a ftut03b ftut04 ftut05 ftut06 ftut07 
fast2: ftut08 ftut09 ftut10 ftut11a  ftut11b ftut12 ftut14
fast3: ftut10 ftut11a  ftut11b 
	#ftut12 ftut14


tut01: 
	$(PYTHON) simo-tut_01-first_calc.py 
  
tut02: 
	$(PYTHON) simo-tut_02-gain_spectra-npsave.py 
ftut02: 
	$(PYTHON) simo-tut_02-gain_spectra-npsave.py  fast=1

tut03a: 
	$(PYTHON) simo-tut_03_1-dispersion-npload.py 
ftut03a: 
	$(PYTHON) simo-tut_03_1-dispersion-npload.py  fast=1

tut03b: 
	$(PYTHON) simo-tut_03_2-dispersion-multicore.py 
ftut03b: 
	$(PYTHON) simo-tut_03_2-dispersion-multicore.py fast=1

tut04: 
	$(PYTHON) simo-tut_04-scan_widths.py 
ftut04: 
	$(PYTHON) simo-tut_04-scan_widths.py fast=1

tut05: 
	$(PYTHON) simo-tut_05-convergence_study.py 
ftut05: 
	$(PYTHON) simo-tut_05-convergence_study.py fast=1

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
	$(PYTHON) simo-tut_12.py
ftut12: 
	$(PYTHON) simo-tut_12.py fast=1

tut13: 
	$(PYTHON) simo-tut_13.py
ftut13: 
	$(PYTHON) simo-tut_13.py fast=1

tut14: 
	$(PYTHON) simo-tut_14-multilayer-fibre.py
ftut14: 
	$(PYTHON) simo-tut_14-multilayer-fibre.py fast=1

clean:
	-rm *.png *.npz
	-rm -rf tut*fields ftut*fields
