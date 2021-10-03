RENG=(R CMD BATCH --vanilla) 
Rfiles=
all: simulation create estimation results clean 

simulation: main/sim.py
	python main/sim.py

create: main/create.py
	python main/create.py

estimation: main/est.py
	python main/est.py

results: main/results.py
	python main/results.py

graphs:

clean: main/clean.sh
	bash main/clean.sh 
