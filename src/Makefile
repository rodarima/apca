CC=gcc -g

all: error big test precision


precision: precision.c
	$(CC) precision.c -lmpfr -o precision

pca.o:
	$(CC) -c -lm pca.c

error: pca.o error.c
	$(CC) -o error -lm pca.o error.c

big: pca.o big.c
	$(CC) -o big -lm pca.o big.c

test: pca.o test.c
	$(CC) -o test -lm pca.o test.c

clean:
	rm -f pca error big precision *.o
