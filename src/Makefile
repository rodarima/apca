pca:
	$(CC) -o pca -lm pca.c

error: pca.o error.c
	$(CC) -o error -lm pca.o error.c

clean:
	rm -f pca
