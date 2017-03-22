#include "reduction.h"
#include "mpvector.h"

#include <mpfr.h>
#include <stdio.h>

#define N 5

void test(mpfr_prec_t prec, mpfr_rnd_t rnd)
{
	mpfr_t **A, *diag, *offdiag, err, tmp1;
	mpfr_prec_t prec_vec[20];
	float **_A, *_diag, *_offdiag;
	int i, j;

	float AA[N][N] = {
		{5,4,3,2,1},
		{4,6,0,4,3},
		{3,0,7,6,5},
		{2,4,6,8,7},
		{1,3,5,7,9},
	};

	for(i = 0; i<20; i++)
		prec_vec[i] = prec;

	mpfr_init2(tmp1, prec);
	mpfr_init2(err, prec);

	A = mpmatrix_init(N, N, prec);
	diag = mpvector_init(N, prec);
	offdiag = mpvector_init(N, prec);

	_A = matrix_init(N, N);
	_diag = vector_init(N);
	_offdiag = vector_init(N);

	for(i=1; i<=N; i++)
	{
		for(j=1; j<=N; j++)
		{
			mpfr_set_flt(A[i][j], AA[i-1][j-1], rnd);
			_A[i][j] = AA[i-1][j-1];
		}
	}

	printf("\n-------- Float precision ---------\n");
//	
//	printf("\nMatrix _A:\n");
//	matrix_print(_A, N, N);

	tred2(_A, N, _diag, _offdiag);

	printf("\nDiagonal:\n");
	vector_print(_diag, N);
	printf("\nOffdiagonal:\n");
	vector_print(_offdiag, N);

	tqli(_diag, _offdiag, N, _A);

	printf("\nEigenvalues:\n");
	vector_print(_diag, N);

	printf("-------- MPFR precision (%d bits) ----------\n", prec);
//	printf("Matrix A:\n");
//	mpmatrix_print(A, N, N);

	mp_tred2(A, N, diag, offdiag, prec_vec, rnd);

	printf("\nDiagonal:\n");
	mpvector_print(diag, N);
	printf("\nOffdiagonal:\n");
	mpvector_print(offdiag, N);
//
	mp_tqli(diag, offdiag, N, A, prec_vec, rnd);

	printf("\nEigenvalues:\n");
	mpvector_print(diag, N);

	//printf("\n-------- Differences ---------\n");
	mpfr_set_flt(err, 0.0, rnd);
	for(i=1; i<=N; i++)
	{
		mpfr_set_flt(tmp1, _diag[i], rnd);
		mpfr_sub(tmp1, tmp1, diag[i], rnd);
		mpfr_mul(tmp1, tmp1, tmp1, rnd);
		mpfr_add(err, err, tmp1, rnd);
		//mpfr_printf("%.15Rf ", tmp1);
	}
	mpfr_sqrt(err, err, rnd);
	//mpfr_printf("\n");
	mpfr_printf("%3d bits, error %.15Rg\n", prec, err);

	mpfr_clear(tmp1);
	mpfr_clear(err);

	mpvector_free(diag, N);
	mpvector_free(offdiag, N);
	mpmatrix_free(A, N, N);

	vector_free(_diag, N);
	vector_free(_offdiag, N);
	matrix_free(_A, N, N);

}

#define PREC_MIN 24
#define PREC_MAX 100
#define PREC_GOLD 500
#define RND_MODE MPFR_RNDN

int main()
{
	mpfr_prec_t b;

	for(b = PREC_MIN; b <= PREC_MAX; b+=10)
	{
		test(b, RND_MODE);
	}

	return 0;
}
