#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "pca.h"
#include <math.h>

//#define N 10
int N = 10;

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

#define ITER 100

struct eig_t
{
	float **A;
	float *diag;
	float *offdiag;
	int n;
};


int test_hh()
{
	float **a, **e, **s;
	float *a_diag, *a_offdiag, *s_diag, *s_offdiag;
	int i, j, k;
	float enorm, error;

	a = symmat(N);
	e = symmat(N);
	s = symmat(N);

	if(!a || !e || !s) exit(1);

	a_diag = vector(N);
	a_offdiag = vector(N);
	s_diag = vector(N);
	s_offdiag = vector(N);

	for(k = 0; k < ITER; k++)
	{
		/* Genererate random matrix A */
		rand_symmat(a, N, -1.0, +1.0);

		/* Some noise */
		rand_symmat(e, N, -0.1, +0.1);

		/* Compute the sum */
		symmat_sum(a, e, s, N);

		/* Compute bound as norm 2 */
		//enorm = norm2(e, n, n);
		
		/* Tridiagonalization of a */
		tred2(a, N, a_diag, a_offdiag);

		/* Now with noise a+e */
		tred2(s, N, s_diag, s_offdiag);

		error = 0.0;
		/* Compute the difference */
		for(i = 1; i <= N; i++)
		{
			error += (a_diag[i] - s_diag[i]) * 
				(a_diag[i] - s_diag[i]);
		}
		error = sqrt(error);

		printf("%03d iteration, error %g\n", k, error);
	}

	free_vector(s_diag, N);
	free_vector(s_offdiag, N);
	free_vector(a_diag, N);
	free_vector(a_offdiag, N);

	free_symmat(s, N);
	free_symmat(e, N);
	free_symmat(a, N);
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

	N = atoi(argv[1]);

	test_hh();

	return 0;
}
