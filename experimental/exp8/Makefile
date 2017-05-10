CC=gcc

# The /apps directory is for compatibility with the cluster
MINOTAURO_CFLAGS=-I/apps/MPFR/3.1.5/include -L/apps/MPFR/3.1.5/lib

CFLAGS=-O2 -Wall $(MINOTAURO_CFLAGS)

# More minotauro problems
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/MPFR/3.1.5/lib/

LDFLAGS=-lm -lmpfr -lgmp

OBJECTS=8a 8b 8c slides.pdf 8d


all: $(OBJECTS)

clean: clean-data
	rm -f *.o
	rm -f $(OBJECTS)

clean-data:
	rm -f data/*

%.pdf: %.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode -shell-escape" -use-make $<


