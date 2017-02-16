pca:
	$(CC) -o pca -lm pca.c

error: pca.o error.c
	$(CC) -o error -lm pca.o error.c

big: pca.o big.c
	$(CC) -o big -lm pca.o big.c

clean:
	rm -f pca
