all: simulation create

simulation: main/sim.py
	python main/sim.py

create: main/create.py
	python main/create.py


