N=1200
M=50000

all: compile run graph doc

compile: prog.cpp
	g++ prog.cpp -o prog -O2

run: prog
	echo "$N $M" | ./prog > data

graph:
	gnuplot prog.gplt

doc:
	pdflatex raport1.tex

clean:
	rm -f prog *.aux *.log *.toc

mrproper: clean
	rm -f *.dat prog.*.tex
