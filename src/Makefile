CC=gcc

CFLAGS=-O2 -Wall -I/apps/MPFR/3.1.5/include -L/apps/MPFR/3.1.5/lib

LDFLAGS=-lm -lmpfr -lgmp

OBJECTS=error big test precision test_ql test_mpvector test_householder errhh \
errcond errn exp6 exp7 test2 exp7a exp7b exp8a exp8b exp8c

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/MPFR/3.1.5/lib/

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

exp7a: mpvector.o reduction.o err.o

exp7b: mpvector.o reduction.o err.o

exp8a: mpvector.o reduction.o err.o

exp8b: mpvector.o reduction.o err.o

exp8c: mpvector.o reduction.o err.o

%.png:
	gnuplot $(basename $@).gp
	sxiv $@

error: pca.o error.c

big: pca.o big.c

test2: pca.o test2.c

clean:
	rm -f *.o
	rm -f $(OBJECTS)
