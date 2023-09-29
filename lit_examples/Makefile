
PYTHON=python3

all: lit01 lit02 lit03 lit04a lit04b lit05 lit06a lit06b lit07 lit08 lit09 lit10
slow: lit04a lit04b lit05 lit06a lit06b lit07 lit08 lit09 lit10

#allfast: flit01 flit02 flit03 flit04-1 flit04-2 flit05 flit06a flit06b flit07 flit08 flit09 flit10
allfast: flit03 flit04a flit04b flit05 flit06a flit06b flit07 flit08 flit09 flit10
fast1: flit06a flit06b flit07 flit08 flit09 flit10

lit01: 
	$(PYTHON) simo-lit_01-Laude-AIPAdv_2013-silica.py
flit01: 
	$(PYTHON) simo-lit_01-Laude-AIPAdv_2013-silica.py fast=1

lit02: 
	$(PYTHON) simo-lit_02-Laude-AIPAdv_2013-silicon.py
flit02: 
	$(PYTHON) simo-lit_02-Laude-AIPAdv_2013-silicon.py fast=1

lit03: 
	$(PYTHON) simo-lit_03-Beugnot-NatComm_2014.py
flit03: 
	$(PYTHON) simo-lit_03-Beugnot-NatComm_2014.py fast=1

lit04a: 
	$(PYTHON) simo-lit_04-no_pillar-Van_Laer-NP_2015.py
flit04a: 
	$(PYTHON) simo-lit_04-no_pillar-Van_Laer-NP_2015.py fast=1

lit04b: 
	$(PYTHON) simo-lit_04-pillar-Van_Laer-NP_2015.py
flit04b: 
	$(PYTHON) simo-lit_04-pillar-Van_Laer-NP_2015.py fast=1

lit05: 
	$(PYTHON) simo-lit_05-Van_Laer-NJP_2015.py
flit05: 
	$(PYTHON) simo-lit_05-Van_Laer-NJP_2015.py fast=1

lit06a: 
	$(PYTHON) simo-lit_06_1-Florez-NatComm_2016-d550nm.py
flit06a: 
	$(PYTHON) simo-lit_06_1-Florez-NatComm_2016-d550nm.py fast=1

lit06b: 
	$(PYTHON) simo-lit_06_2-Florez-NatComm_2016-d1160nm.py
flit06b: 
	$(PYTHON) simo-lit_06_2-Florez-NatComm_2016-d1160nm.py fast=1

lit07: 
	$(PYTHON) simo-lit_07-Kittlaus-NatPhot_2016.py
flit07: 
	$(PYTHON) simo-lit_07-Kittlaus-NatPhot_2016.py fast=1

lit08: 
	$(PYTHON) simo-lit_08-Kittlaus-NatComm_2017.py
flit08: 
	$(PYTHON) simo-lit_08-Kittlaus-NatComm_2017.py fast=1

lit09: 
	$(PYTHON) simo-lit_09-Morrison-Optica_2017.py
flit09: 
	$(PYTHON) simo-lit_09-Morrison-Optica_2017.py fast=1

lit10: 
#$(PYTHON) simo-lit_10-Wolff-OptExpress-2014.py
	$(PYTHON) simo-lit_10a-Wolff-OptExpress-2014.py
flit10: 
	$(PYTHON) simo-lit_10-Wolff-OptExpress-2014.py fast=1
	#$(PYTHON) simo-lit_10a-Wolff-OptExpress-2014.py fast=1

clean:
	-rm -rf lit_*-fields
	-rm -rf flit_*-fields
	-rm lit* flit*
