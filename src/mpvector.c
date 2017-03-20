#include "mpvector.h"

#include <stdlib.h>
#include <mpfr.h>
#include <stdio.h>
#include "err.h"

struct mptridiag_t *mptridiag_init(int n, mpfr_prec_t *prec, mpfr_rnd_t rnd)
{
	struct mptridiag_t *td = malloc(sizeof(*td));
	if(!td) erhand("Allocation failure in mptridiag_init().");
	td->A = mpmatrix_init(n, n, prec[0]);
	td->diag = mpvector_init(n, prec[1]);
	td->offdiag = mpvector_init(n, prec[2]);
	td->n = n;
	td->prec = prec[0];
	td->rnd = rnd;

	return td;
}

void mptridiag_free(struct mptridiag_t *td)
{
	mpmatrix_free(td->A, td->n, td->n);
	mpvector_free(td->diag, td->n);
	mpvector_free(td->offdiag, td->n);
	free(td);
}

/**  Allocation of vector storage  ***********************************/

/* Allocates a float vector with range [1..n]. */
mpfr_t *mpvector_init(int n, mpfr_prec_t prec)
{

	mpfr_t *v;
	int i;

	v = (mpfr_t *) malloc ((unsigned) n * sizeof(mpfr_t));
	if (!v) erhand("Allocation failure in vector().");
	v--;
	for(i = 1; i <= n; i++)
	{
		mpfr_init2(v[i], prec);
	}
	return v;

}

void mpvector_diff(mpfr_t err, mpfr_t *a, mpfr_t *b, int n, mpfr_prec_t prec, mpfr_rnd_t rnd)
{
	int i;
	mpfr_t tmp1;

	mpfr_init2(tmp1, prec);
	
	mpfr_set_d(err, 0.0, rnd);

	for(i = 1; i<=n; i++)
	{
		mpfr_sub(tmp1, a[i], b[i], rnd);
		mpfr_mul(tmp1, tmp1, tmp1, rnd);
		mpfr_add(err, err, tmp1, rnd);
	}
	mpfr_sqrt(err, err, rnd);

	mpfr_clear(tmp1);
}
/**  Allocation of float matrix storage  *****************************/

/* Allocate a float matrix with range [1..n][1..m]. */
mpfr_t **mpmatrix_init(int n, int m, mpfr_prec_t prec)
{
	int i, j;
	mpfr_t **mat;
	mpfr_t *ptr;

	/* Allocate pointers to rows. */
	mat = (mpfr_t **) malloc((unsigned) (n)*sizeof(mpfr_t*));
	if (!mat) erhand("Allocation failure 1 in matrix().");
	mat--;

	/* Allocate rows and set pointers to them. */
	for (i = 1; i <= n; i++)
	{
		mat[i] = (mpfr_t *) malloc((unsigned) (m)*sizeof(mpfr_t));
		if (!mat[i]) erhand("Allocation failure 2 in matrix().");
		mat[i]--;
		for (j = 1; j <= m; j++)
		{
			mpfr_init2(mat[i][j], prec);
		}
	}

	/* Return pointer to array of pointers to rows. */
	return mat;

}

/**  Deallocate vector storage  *********************************/

/* Free a float vector allocated by vector(). */
void mpvector_free(mpfr_t *v, int n)
{
	int i;

	for(i = 1; i <= n; i++)
	{
		mpfr_clear(v[i]);
	}

	free((char*) (v+1));
}

/**  Deallocate float matrix storage  ***************************/

/* Free a float matrix allocated by matrix(). */
void mpmatrix_free(mpfr_t **mat, int n, int m)
{
	int i, j;

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= m; j++)
		{
			mpfr_clear(mat[i][j]);
		}
		free ((char*) (mat[i]+1));
	}
	free ((char*) (mat+1));
}


void mpmatrix_print(mpfr_t **mat, int n, int m)
{
	int i, j;

	for(i=1; i<=n; i++)
	{
		for(j=1; j<=m; j++)
		{
			mpfr_printf("%15.12Rf ", mat[i][j]);
		}
		printf("\n");
	}
}

void mpvector_print(mpfr_t *vec, int n)
{
	int i;

	for(i=1; i<=n; i++)
	{
		mpfr_printf("%.15Rf ", vec[i]);
	}
	printf("\n");
}

/**  Allocation of vector storage  ***********************************/

float *vector_init(n)
int n;
/* Allocates a float vector with range [1..n]. */
{

    float *v;

    v = (float *) malloc ((unsigned) n*sizeof(float));
    if (!v) erhand("Allocation failure in vector().");
    return v-1;

}

/**  Allocation of float matrix storage  *****************************/

float **matrix_init(n,m)
int n, m;
/* Allocate a float matrix with range [1..n][1..m]. */
{
    int i;
    float **mat;

    /* Allocate pointers to rows. */
    mat = (float **) malloc((unsigned) (n)*sizeof(float*));
    if (!mat) erhand("Allocation failure 1 in matrix().");
    mat -= 1;

    /* Allocate rows and set pointers to them. */
    for (i = 1; i <= n; i++)
        {
        mat[i] = (float *) malloc((unsigned) (m)*sizeof(float));
        if (!mat[i]) erhand("Allocation failure 2 in matrix().");
        mat[i] -= 1;
        }

     /* Return pointer to array of pointers to rows. */
     return mat;

}

/**  Deallocate vector storage  *********************************/

void vector_free(v,n)
float *v;
int n;
/* Free a float vector allocated by vector(). */
{
   free((char*) (v+1));
}

/**  Deallocate float matrix storage  ***************************/

void matrix_free(mat,n,m)
float **mat;
int n, m;
/* Free a float matrix allocated by matrix(). */
{
   int i;

   for (i = n; i >= 1; i--)
       {
       free ((char*) (mat[i]+1));
       }
   free ((char*) (mat+1));
}

void matrix_print(float **mat, int n, int m)
{
	int i, j;

	for(i=1; i<=n; i++)
	{
		for(j=1; j<=m; j++)
		{
			printf("%15.12f ", mat[i][j]);
		}
		printf("\n");
	}
}

void vector_print(float *vec, int n)
{
	int i;

	for(i=1; i<=n; i++)
	{
		printf("%.15f ", vec[i]);
	}
	printf("\n");
}
