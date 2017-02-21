#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "pca.h"
#include <math.h>

//#define N 10

float norm2(float **a, int n, int m)
{
	int i, j;
	float accum = 0.0;

	for(i=1; i<=n; i++)
	{
		for(j=1; j<=m; j++)
		{
			accum += a[i][j] * a[i][j];
		}
	}
	return sqrt(accum);
}

float frand(float lmin, float lmax)
{
	return (lmax - lmin) * ((float)rand() / RAND_MAX) + lmin;
}

void rand_symmat(float **a, int n, float lmin, float lmax)
{
	int i, j;

	for(i = 1; i <= n; i++)
	{
		for(j = 1; j <= n; j++)
		{
			if(j >= i)
				a[i][j] = frand(lmin, lmax);
			else
				a[i][j] = a[j][i];
		}
	}
	//print_matrix(a, n, n);
}

float **symmat(int n)
{
	return matrix(n, n);
}

void free_symmat(float **a, int n)
{
	free_matrix(a, n, n);
}

void symmat_sum(float **a, float **b, float **c, int n)
{
	int i, j;

	for(i = 1; i <= n; i++)
	{
		for(j = 1; j <= n; j++)
		{
			if(j >= i)
				c[i][j] = a[i][j] + b[i][j];
			else
				c[i][j] = c[j][i];
		}
	}
}

#define ITER 20

struct eig_t
{
	float **A;
	float *diag;
	float *offdiag;
	int n;
};

int test(int n)
{
	struct eig_t eig, eigerr;
	float **noise;
	int i, j, k;
	float enorm, error, error2;

	eig.A = symmat(n);
	eigerr.A = symmat(n);
	noise = symmat(n);

	if(!eig.A || !eigerr.A || !noise) exit(1);

	eig.diag = vector(n);
	eig.offdiag = vector(n);
	eigerr.diag = vector(n);
	eigerr.offdiag = vector(n);

	for(k = 0; k < ITER; k++)
	{
		/* Genererate random matrix A */
		rand_symmat(eig.A, n, -1.0, +1.0);

		/* Some noise */
		rand_symmat(noise, n, -0.1, +0.1);

		/* Compute the sum */
		symmat_sum(eig.A, noise, eigerr.A, n);

		/* Compute bound as norm 2 */
		enorm = norm2(noise, n, n);
		
		/* Tridiagonalization of a */
		tred2(eig.A, n, eig.diag, eig.offdiag);

		/* now with noise a+e */
		tred2(eigerr.A, n, eigerr.diag, eigerr.offdiag);

		error = 0.0;
		/* Compute the difference */
		for(i = 1; i <= n; i++)
		{
			error += (eig.diag[i] - eigerr.diag[i]) * 
				(eig.diag[i] - eigerr.diag[i]);
		}
		error = sqrt(error);


		tqli(eig.diag, eig.offdiag, n, eig.A);
		tqli(eigerr.diag, eigerr.offdiag, n, eigerr.A);

		error2 = 0.0;
		/* Compute the difference */
		for(i = 1; i <= n; i++)
		{
			error2 += (eig.diag[i] - eigerr.diag[i]) * 
				(eig.diag[i] - eigerr.diag[i]);
		}

		error2 = sqrt(error2);
		printf("%03d: error HH:%10.5g  QL:%10.5g  Enorm:%g\n", k, error, error2, enorm);
	}

	free_vector(eigerr.diag, n);
	free_vector(eigerr.offdiag, n);
	free_vector(eig.diag, n);
	free_vector(eig.offdiag, n);

	free_symmat(noise, n);
	free_symmat(eigerr.A, n);
	free_symmat(eig.A, n);

	return 0;
}

void usage(int argc, char *argv[])
{
	printf("Usage: %s n\n", argv[0]);
	printf("Where the matrix is n x n\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc != 2)
		usage(argc, argv);
	
	srand(time(NULL));

	int n = atoi(argv[1]);

	test(n);

	return 0;
}
