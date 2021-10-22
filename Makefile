Rfiles= plots.R
MAIN = sim.py create.py est.py results.py
PY = python 
RENV=R
RFLAGS=CMD BATCH --vanilla
R_FILES=plots.R 
all: simulation create estimation results  

simulation: main/sim.py
	$(PY) main/sim.py

create: main/create.py
	$(PY) main/create.py

estimation: main/est.py
	$(PY) main/est.py

results: main/results.py
	$(PY) main/results.py

graphs: $(R_FILES)
	$(RENV) $(RFLAGS) $(R_FILES) 

clean: main/clean.sh
	bash main/clean.sh 

RUNALL: main/$(MAIN)
	$(PY)  main/$(MAIN) 
