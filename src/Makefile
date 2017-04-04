CC=gcc

CFLAGS=-g -Wall

LDFLAGS=-lm -lmpfr -lgmp

OBJECTS=error big test precision test_ql test_mpvector test_householder errhh \
errcond errn exp6 exp7 test2

all: error big test precision

test: test_mpvector test_householder test_ql

test_mpvector: err.o mpvector.o

test_householder: err.o mpvector.o reduction.o

test_ql: err.o mpvector.o reduction.o

mpvector.o: err.o

precision: precision.c

errhh: mpvector.o reduction.o err.o

errcond: mpvector.o reduction.o err.o

errn: mpvector.o reduction.o err.o

exp6: mpvector.o reduction.o err.o

exp7: mpvector.o reduction.o err.o

%.png:
	gnuplot $(basename $@).gp
	sxiv $@

error: pca.o error.c

big: pca.o big.c

test2: pca.o test2.c

clean:
	rm -f *.o
	rm -f $(OBJECTS)
