
PYTHON=python

all: tut01 tut02 tut03a tut03b tut04 tut05 tut06 tut07 \
				tut08 tut09 tut11a  tut11b tut12 tut13 tut14

allfast: tut01 tut02 tut03a ftut03b ftut04 ftut05 ftut06 ftut07 \
				ftut08 ftut09 ftut11a  ftut12 ftut13 
allfastb: tut01 tut02 tut03a ftut03b ftut04 ftut05 ftut06 ftut07 \
				ftut08 ftut09 ftut11a  ftut12 ftut13 

broken: tut10

slow: tut04 tut05 tut06 tut07 tut08 tut09 tut11a  tut11b tut12 tut13 tut14

goodtest: tut01 ftut02 ftut06 ftut08


tut01:
	$(PYTHON) sim-tut_01-first_calc.py

tut02:
	$(PYTHON) sim-tut_02-gain_spectra-npsave.py
ftut02:
	$(PYTHON) sim-tut_02-gain_spectra-npsave.py  fast=1

tut03a:
	$(PYTHON) sim-tut_03_1-dispersion-npload.py
ftut03a:
	$(PYTHON) sim-tut_03_1-dispersion-npload.py  fast=1

tut03b:
	$(PYTHON) sim-tut_03_2-dispersion-multicore.py
ftut03b:
	$(PYTHON) sim-tut_03_2-dispersion-multicore.py fast=1

tut04:
	$(PYTHON) sim-tut_04-scan_widths.py
ftut04:
	$(PYTHON) sim-tut_04-scan_widths.py fast=1

tut05:
	$(PYTHON) sim-tut_05-convergence_study.py
ftut05:
	$(PYTHON) sim-tut_05-convergence_study.py fast=1

tut06:
	$(PYTHON) sim-tut_06-silica_nanowire.py
ftut06:
	$(PYTHON) sim-tut_06-silica_nanowire.py fast=1

tut06b:
	$(PYTHON) sim-tut_06b-silicon_triangular_nanowire.py
ftut06b:
	$(PYTHON) sim-tut_06b-silicon_triangular_nanowire.py fast=1

tut07:
	$(PYTHON) sim-tut_07-slot.py
ftut07:
	$(PYTHON) sim-tut_07-slot.py fast=1

tut08:
	$(PYTHON) sim-tut_08-slot_coated-scan.py
ftut08:
	$(PYTHON) sim-tut_08-slot_coated-scan.py fast=1

tut09:
	$(PYTHON) sim-tut_09-anisotropy.py
ftut09:
	$(PYTHON) sim-tut_09-anisotropy.py fast=1

tut10:
	$(PYTHON) sim-tut_10-onion.py
ftut10:
	$(PYTHON) sim-tut_10-onion.py fast=1

tut11a:
	$(PYTHON) sim-tut_11a-onion2.py
ftut11a:
	$(PYTHON) sim-tut_11a-onion2.py fast=1

tut11b:
	$(PYTHON) sim-tut_11b-onion3.py
ftut11b:
	$(PYTHON) sim-tut_11b-onion3.py fast=1

tut12:
	$(PYTHON) sim-tut_12.py
	$(PYTHON) sim-tut_12.py smf28
ftut12:
	$(PYTHON) sim-tut_12.py fast=1
	$(PYTHON) sim-tut_12.py fast=1 smf28

tut13:
	$(PYTHON) sim-tut_13.py
ftut13:
	$(PYTHON) sim-tut_13.py fast=1

tut14:
	$(PYTHON) sim-tut_14-multilayer-fibre.py
ftut14:
	$(PYTHON) sim-tut_14-multilayer-fibre.py fast=1

tut17:
	$(PYTHON) sim-tut_17-rotations.py
ftut17:
	$(PYTHON) sim-tut_17-rotations.py fast=1

tut18:
	$(PYTHON) sim-tut_18-sbs-in-smf28.py
ftut18:
	$(PYTHON) sim-tut_18-sbs-in-smf28.py fast=1


tut19:
	$(PYTHON) sim-tut_19_waveguide_profiles.py
ftut19:
	$(PYTHON) sim-tut_19_waveguide_profiles.py fast=1


clean:
	-rm *.png *.npz
	-rm -rf tut*fields ftut*fields __pycache__
	-rm -rf tut*out tut*.txt ftut*.txt

