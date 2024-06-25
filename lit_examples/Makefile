
PYTHON=python3

all: lit01 lit02 lit03 lit04a lit04b lit05 lit06a lit06b lit07 lit08 lit09 lit10
slow: lit04a lit04b lit05 lit06a lit06b lit07 lit08 lit09 lit10

allfast: flit01 flit02 flit03 flit04a flit04b flit05 flit06a flit06b flit07 flit08 flit09 flit10
#allfast: flit03 flit04a flit04b flit05 flit06a flit06b flit07 flit08 flit09 flit10
fast1: flit08 flit09 flit10

lit01: 
	$(PYTHON) sim-lit_01-Laude-AIPAdv_2013-silica.py
flit01: 
	$(PYTHON) sim-lit_01-Laude-AIPAdv_2013-silica.py fast=1

lit02: 
	$(PYTHON) sim-lit_02-Laude-AIPAdv_2013-silicon.py
flit02: 
	$(PYTHON) sim-lit_02-Laude-AIPAdv_2013-silicon.py fast=1

lit03: 
	$(PYTHON) sim-lit_03-Beugnot-NatComm_2014.py
	$(PYTHON) sim-lit_03-Beugnot-NatComm_2014.py widerange
flit03: 
	$(PYTHON) sim-lit_03-Beugnot-NatComm_2014.py fast=1
	$(PYTHON) sim-lit_03-Beugnot-NatComm_2014.py fast=1 widerange

lit04a: 
	$(PYTHON) sim-lit_04-no_pillar-Van_Laer-NP_2015.py
flit04a: 
	$(PYTHON) sim-lit_04-no_pillar-Van_Laer-NP_2015.py fast=1

lit04b: 
	$(PYTHON) sim-lit_04-pillar-Van_Laer-NP_2015.py
flit04b: 
	$(PYTHON) sim-lit_04-pillar-Van_Laer-NP_2015.py fast=1

lit05: 
	$(PYTHON) sim-lit_05-Van_Laer-NJP_2015.py
flit05: 
	$(PYTHON) sim-lit_05-Van_Laer-NJP_2015.py fast=1

lit06a: 
	$(PYTHON) sim-lit_06_1-Florez-NatComm_2016-d550nm.py
flit06a: 
	$(PYTHON) sim-lit_06_1-Florez-NatComm_2016-d550nm.py fast=1

lit06b: 
	$(PYTHON) sim-lit_06_2-Florez-NatComm_2016-d1160nm.py
flit06b: 
	$(PYTHON) sim-lit_06_2-Florez-NatComm_2016-d1160nm.py fast=1

lit07: 
	$(PYTHON) sim-lit_07-Kittlaus-NatPhot_2016.py
flit07: 
	$(PYTHON) sim-lit_07-Kittlaus-NatPhot_2016.py fast=1

lit08: 
	$(PYTHON) sim-lit_08-Kittlaus-NatComm_2017.py
flit08: 
	$(PYTHON) sim-lit_08-Kittlaus-NatComm_2017.py fast=1

lit09: 
	$(PYTHON) sim-lit_09-Morrison-Optica_2017.py
flit09: 
	$(PYTHON) sim-lit_09-Morrison-Optica_2017.py fast=1

lit10: 
	$(PYTHON) sim-lit_10a-Wolff-OptExpress-2014.py
	$(PYTHON) sim-lit_10b-Wolff-OptExpress-2014.py
flit10: 
	$(PYTHON) sim-lit_10a-Wolff-OptExpress-2014.py fast=1
	$(PYTHON) sim-lit_10b-Wolff-OptExpress-2014.py fast=1

clean:
	-rm -rf lit_*-fields *-out __pycache__
	-rm -rf flit_*-fields
	-rm lit* flit*
