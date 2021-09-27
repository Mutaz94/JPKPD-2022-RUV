
run: main.py
	python main.py
graphs: src/help/plots.R
	R CMD BATCH --vanilla plots.R
clean: clean.sh
	bash clean.sh 


