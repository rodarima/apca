#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "pca.h"
#include <math.h>

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


float frand()
{
	return ((float) rand()) / ((float) RAND_MAX) * 2.0 - 1.0;
}

int main(int argc, char *argv[])
{
	float **A, **E, **SUM;
	float *evals, *interm;
	float Enorm;
	int i, j;

	srand(time(NULL));

	A = matrix(2, 2);
	E = matrix(2, 2);
	SUM = matrix(2, 2);

    	evals = vector(2);     /* Storage alloc. for vector of eigenvalues */
    	interm = vector(2);    /* Storage alloc. for 'intermediate' vector */

	float sA[2][2] = {{6.8, 2.4}, {2.4, 8.2}};
	float sE[2][2] = {{.002, .003}, {.003, .001}};

	for (i = 0; i < 2; i++)
	{
		for(j = 0; j < 2; j++)
		{
			A[i+1][j+1] = sA[i][j];
			if(i >= j)
				E[i+1][j+1] = frand() * 0.005;
		}
	}
	for (i = 0; i < 2; i++)
	{
		for(j = 0; j < 2; j++)
		{
			if(i < j)
				E[i+1][j+1] = E[j+1][i+1];
		}
	}
	for (i = 0; i < 2; i++)
	{
		for(j = 0; j < 2; j++)
		{
			SUM[i+1][j+1] = A[i+1][j+1] + E[i+1][j+1];
		}
	}

	Enorm = norm2(E, 2, 2);

	printf("Matrix A:\n");
	print_matrix(A, 2, 2);
	printf("Matrix E:\n");
	print_matrix(E, 2, 2);
	printf("Matrix E norm2 = %f:\n", Enorm);
	printf("Matrix SUM:\n");
	print_matrix(SUM, 2, 2);

	printf("Householder tridiagonalization\n");
    	tred2(SUM, 2, evals, interm);  /* Triangular decomposition */
	printf("diagonal:\n");
	print_vector(evals, 2);
	printf("off-diagonal:\n");
	print_vector(interm, 2);
	printf("Q:\n");
	print_matrix(SUM, 2, 2);
	printf("Implicit QR diagonalization\n");
	tqli(evals, interm, 2, SUM);   /* Reduction of sym. trid. matrix */
	printf("Eigenvalues:\n");
	for (j = 2; j >= 1; j--) {
		printf("%60.50f\n", evals[j]);
	}

	/*
	double eigs[] = {
		4.998759651901539591278833540854975581169128417968750,
		10.0042403480984596342295844806358218193054199218750
	};*/
	double eigs[] = { 5.0, 10.0 };
	double eigerr;
	
	printf("Error in eigenvalues: (should be <= %f)\n", Enorm);
	for (j = 2; j >= 1; j--) {
		eigerr = eigs[j-1] - ((double) evals[j]);
		printf("%20.10G", eigerr);
		if(eigerr > Enorm) printf(" <- ERROR, eigenvalue exeeds maximum");
		printf("\n");
	}

	free_vector(interm, 2);
	free_vector(evals, 2);
	free_matrix(SUM, 2, 2);
	free_matrix(E, 2, 2);
	free_matrix(A, 2, 2);

	return 0;
}
