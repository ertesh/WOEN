all: compile run graph doc

compile: prog.cpp
	g++ prog.cpp -o prog -O2 -Wall

run: prog
	./prog

graph:
	gnuplot prog.gplt

doc:
	pdflatex raport1.tex

clean:
	rm -f prog *.aux *.log *.toc

mrproper: clean
	rm -f *.dat prog.*.tex
