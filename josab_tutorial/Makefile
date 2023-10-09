
PYTHON=python3

.PHONY: all
all: josab01 josab02 josab03 josab04 josab05 josab06 
all2: josab02 josab03 josab04 josab05 josab06 

.PHONY: allfast
allfast: fjosab01 fjosab02 fjosab03 fjosab04 fjosab05 fjosab06

josab01: 
	$(PYTHON) simo-josab-01.py
fjosab01: 
	$(PYTHON) simo-josab-01.py fast=1

josab02: 
	$(PYTHON)  simo-josab-02.py 
	$(PYTHON) simo-josab-02b-acdisp.py
fjosab02: 
	$(PYTHON)  simo-josab-02.py fast=1
	$(PYTHON) simo-josab-02b-acdisp.py fast=1

josab03: 
	$(PYTHON) simo-josab-03.py
fjosab03: 
	$(PYTHON) simo-josab-03.py fast=1

josab04: 
	$(PYTHON) simo-josab-04i.py
fjosab04: 
	$(PYTHON) simo-josab-04.py fast=1

josab05: 
	$(PYTHON) simo-josab-05.py
fjosab05: 
	$(PYTHON) simo-josab-05.py fast=1

josab06: 
	$(PYTHON) simo-josab-06.py
fjosab06: 
	$(PYTHON) simo-josab-06.py fast=1

clean:
	-rm -rf josab* fjosab*

fjosab11:
	$(PYTHON) simo-josab-fig13-disp.py fast=1

