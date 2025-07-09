
PYTHON=python3

.PHONY: all
all: josab01 josab02 josab03 josab04 josab05 josab06 
all2: josab04 josab05 josab06 

.PHONY: allfast
allfast: fjosab01 fjosab02 fjosab03 fjosab04 fjosab05 fjosab06
allfastb: fjosab04 fjosab05 fjosab06

josab01: 
	$(PYTHON) sim-josab-01.py
fjosab01: 
	$(PYTHON) sim-josab-01.py fast=1

josab02: 
	$(PYTHON) sim-josab-02.py 
	$(PYTHON) sim-josab-02b-acdisp.py
fjosab02: 
	$(PYTHON) sim-josab-02.py fast=1
	$(PYTHON) sim-josab-02b-acdisp.py fast=1

josab03: 
	$(PYTHON) sim-josab-03.py
fjosab03: 
	$(PYTHON) sim-josab-03.py fast=1

josab04: 
	$(PYTHON) sim-josab-04.py
fjosab04: 
	$(PYTHON) sim-josab-04.py fast=1

josab05: 
	$(PYTHON) sim-josab-05.py
fjosab05: 
	$(PYTHON) sim-josab-05.py fast=1

josab06: 
	$(PYTHON) sim-josab-06.py
fjosab06: 
	$(PYTHON) sim-josab-06.py fast=1

clean:
	-rm -rf josab* fjosab*

fjosab11:
	$(PYTHON) sim-josab-fig13-disp.py fast=1

