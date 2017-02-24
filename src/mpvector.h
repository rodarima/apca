#ifndef _MPVECTOR_H_
#define _MPVECTOR_H_

#include <mpfr.h>

struct mptridiag_t
{
	mpfr_t **A;
	mpfr_t *diag;
	mpfr_t *offdiag;
	int n;
	mpfr_prec_t prec;
	mpfr_rnd_t rnd;
};

struct mptridiag_t *mptridiag_init(int n, mpfr_prec_t prec, mpfr_rnd_t rnd);
void mptridiag_free(struct mptridiag_t *td);

/**  Allocation of vector storage  ***********************************/

/* Allocates a float vector with range [1..n]. */
mpfr_t *mpvector_init(int n, mpfr_prec_t prec);

void mpvector_diff(mpfr_t err, mpfr_t *a, mpfr_t *b, int n, mpfr_prec_t prec, mpfr_rnd_t rnd);

/* Allocate a float matrix with range [1..n][1..m]. */
mpfr_t **mpmatrix_init(int n, int m, mpfr_prec_t prec);

/**  Deallocate vector storage  *********************************/

/* Free a float vector allocated by vector(). */
void mpvector_free(mpfr_t *v, int n);

/**  Deallocate float matrix storage  ***************************/

/* Free a float matrix allocated by matrix(). */
void mpmatrix_free(mpfr_t **mat, int n, int m);


void mpmatrix_print(mpfr_t **mat, int n, int m);
void mpvector_print(mpfr_t *vec, int n);






/* Allocates a float vector with range [1..n]. */
float *vector_init(int n);

/* Allocate a float matrix with range [1..n][1..m]. */
float **matrix_init(int n, int m);

/* Free a float vector allocated by vector(). */
void vector_free(float *v, int n);

/**  Deallocate float matrix storage  ***************************/

/* Free a float matrix allocated by matrix(). */
void matrix_free(float **mat, int n, int m);


void matrix_print(float **mat, int n, int m);
void vector_print(float *vec, int n);
#endif /* _MPVECTOR_H_ */
