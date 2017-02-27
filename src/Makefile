CC=gcc -g

LDFLAGS=-lm -lmpfr -lgmp

all: error big test precision

test: test_mpvector test_householder

test_mpvector: err.o mpvector.o

test_householder: err.o mpvector.o reduction.o

mpvector.o: err.o

precision: precision.c

errhh: mpvector.o reduction.o err.o

%.png:
	gnuplot $(basename $@).gp
	sxiv $@

error: pca.o error.c

big: pca.o big.c

test2: pca.o test2.c

clean:
	rm -f pca error big precision *.o
