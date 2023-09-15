
PYTHON=python3

.PHONY: all
all: josab01 josab02 josab03 josab04 josab05 josab06 
all2: josab02 josab03 josab04 josab05 josab06 

.PHONY: allfast
allfast: fjosab01 fjosab02 fjosab03 fjosab04 fjosab05 fjosab06

josab01: 
	$(PYTHON) simo-josab-01-BSBS-1umcylwg-SiO2.py
fjosab01: 
	$(PYTHON) simo-josab-01-BSBS-1umcylwg-SiO2.py fast=1

josab02: 
	#$(PYTHON) simo-josab-02-BSBS-450x200nmrectwg-Si.py
	$(PYTHON) simo-josab-02b-BSBS-acbands-450x200nmrectwg-Si.py
fjosab02: 
	$(PYTHON)  simo-josab-02-BSBS-450x200nmrectwg-Si.py fast=1
	$(PYTHON) simo-josab-02b-BSBS-acbands-450x200nmrectwg-Si.py fast=1

josab03: 
	$(PYTHON) simo-josab-03-FSBS-1umcylwg-SiO2.py
fjosab03: 
	$(PYTHON) simo-josab-03-FSBS-1umcylwg-SiO2.py fast=1

josab04: 
	$(PYTHON) simo-josab-04-FSBS-450x200nmrectwg-Si.py
fjosab04: 
	$(PYTHON) simo-josab-04-FSBS-450x200nmrectwg-Si.py fast=1

josab05: 
	$(PYTHON) simo-josab-05-IFSBS-1umcylwg-SiO2.py
fjosab05: 
	$(PYTHON) simo-josab-05-IFSBS-1umcylwg-SiO2.py fast=1

josab06: 
	$(PYTHON) simo-josab-06-IFSBS-450x200nmrectwg-Si.py
fjosab06: 
	$(PYTHON) simo-josab-06-IFSBS-450x200nmrectwg-Si.py fast=1

clean:
	-rm -rf josab* fjosab*

