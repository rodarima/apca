#include "reduction.h"
#include "mpvector.h"

#include <mpfr.h>
#include <stdio.h>

#define N 5
#define PREC 24 /* this machine is NOT IEEE 754 compliant */

struct rnd_t
{
	mpfr_rnd_t rnd;
	char *rnd_str;
};

#define N_RND_MODES 5
struct rnd_t rnd_modes[N_RND_MODES] = {
	{MPFR_RNDN, "N"},
	{MPFR_RNDZ, "Z"},
	{MPFR_RNDU, "U"},
	{MPFR_RNDD, "D"},
	{MPFR_RNDA, "A"},
};

void test(mpfr_prec_t prec, mpfr_rnd_t rnd, char *rnd_str)
{
	mpfr_t **A, *diag, *offdiag, err, tmp1;
	float **_A, *_diag, *_offdiag;
	int i, j;

	float AA[N][N] = {
		{5,4,3,2,1},
		{4,6,0,4,3},
		{3,0,7,6,5},
		{2,4,6,8,7},
		{1,3,5,7,9},
	};

	A = mpmatrix_init(N, N, prec);
	diag = mpvector_init(N, prec);
	offdiag = mpvector_init(N, prec);
	mpfr_init2(tmp1, prec);
	mpfr_init2(err, prec);
	
	_A = matrix_init(N, N);
	_diag = vector_init(N);
	_offdiag = vector_init(N);

	for(i=1; i<=N; i++)
		for(j=1; j<=N; j++)
		{
			mpfr_set_flt(A[i][j], AA[i-1][j-1], rnd);
			_A[i][j] = AA[i-1][j-1];
		}
	
//	printf("-------- MPFR precision (%d bits) ----------\n", prec);
//	printf("Matrix A:\n");
//	mpmatrix_print(A, N, N);

	mp_tred2(A, N, diag, offdiag, prec, rnd);

//	printf("\nDiagonal:\n");
//	mpvector_print(diag, N);
//	printf("\nOffdiagonal:\n");
//	mpvector_print(offdiag, N);
//
//	printf("\n-------- Float precision ---------\n");
//	
//	printf("\nMatrix _A:\n");
//	matrix_print(_A, N, N);

	tred2(_A, N, _diag, _offdiag);

//	printf("\nDiagonal:\n");
//	vector_print(_diag, N);
//	printf("\nOffdiagonal:\n");
//	vector_print(_offdiag, N);

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
	mpfr_printf("%3d bits, error %.15Rg %s\n", prec, err, rnd_str);

	mpfr_clear(tmp1);
	mpfr_clear(err);

	mpvector_free(diag, N);
	mpvector_free(offdiag, N);
	mpmatrix_free(A, N, N);
	
	vector_free(_diag, N);
	vector_free(_offdiag, N);
	matrix_free(_A, N, N);

}

int main()
{
	int i, j;

	for(j = 15; j <= 30; j++)
		for(i = 0; i < N_RND_MODES; i++)
			test(j, rnd_modes[i].rnd, rnd_modes[i].rnd_str);

	return 0;
}
