#include "reduction.h"
#include "mpvector.h"

#include <mpfr.h>
#include <stdio.h>

#define N 5
#define PREC 24 /* this machine is NOT IEEE 754 compliant */
#define RND_DIR MPFR_RNDN

#define GOLD_PREC 500
#define PREC_MIN 2
#define PREC_MAX 100
#define RUNS 10000

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

void hh(struct mptridiag_t *td)
{
	mp_tred2(td->A, td->n, td->diag, td->offdiag, td->prec, td->rnd);
}

void errhh(struct mptridiag_t *gold, struct mptridiag_t *target)
{
	mpfr_t err, err1;
	hh(target);
	
	mpfr_init2(err1, gold->prec);
	mpfr_init2(err, gold->prec);

	mpfr_sub(err, gold->diag[1], target->diag[1], gold->rnd);
	mpfr_abs(err1, err, gold->rnd);
	
	//mpvector_diff(err, gold->diag, target->diag, gold->n, gold->prec, gold->rnd);

	//mpfr_printf("errhh:%d bits, error %.20Rg\n", target->prec, err);
	mpfr_printf("%d\t%.20Rg\t%.20Rg\n", target->prec, err1, err);

	mpfr_clear(err1);
	mpfr_clear(err);
}

void rand_mpmat(struct mptridiag_t *rop, gmp_randstate_t state)
{
	int i, j;
	mpfr_t **a = rop->A;

	for(i = 1; i <= rop->n; i++)
	{
		for(j = 1; j <= rop->n; j++)
		{
			if(j >= i)
			{
				mpfr_urandom(a[i][j], state, rop->rnd);
			}
			else
			{
				mpfr_set(a[i][j], a[j][i], rop->rnd);
			}
		}
	}
}

void copy_mpmat(struct mptridiag_t *rop, struct mptridiag_t *src)
{
	int i,j;

	for(i = 1; i <= rop->n; i++)
		for(j = 1; j <= rop->n; j++)
			mpfr_set(rop->A[i][j], src->A[i][j], rop->rnd);
}

void load_matrix(struct mptridiag_t *rop, double A[N][N])
{
	int i,j;

	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			mpfr_set_d(rop->A[i+1][j+1], A[i][j], rop->rnd);
}

void test_rnd(gmp_randstate_t state)
{
	int i, n = N;
	struct mptridiag_t *gold, *target, *ref;

	gold = mptridiag_init(n, GOLD_PREC, RND_DIR);
	ref = mptridiag_init(n, GOLD_PREC, RND_DIR);
	//load_matrix(gold, AA);
	rand_mpmat(ref, state);
	copy_mpmat(gold, ref);

	hh(gold);
	//printf("Gold diagonal:\n");
	//mpvector_print(gold->diag, gold->n);

	for(i = PREC_MIN; i<=PREC_MAX; i++)
	{
		target = mptridiag_init(n, i, RND_DIR);
		//load_matrix(target, AA);
		copy_mpmat(target, ref);
		errhh(gold, target);
		//printf("Target diagonal:\n");
		//mpvector_print(target->diag, target->n);
		mptridiag_free(target);
	}
	
	mptridiag_free(gold);

}

int main()
{
	int i;
	gmp_randstate_t state;

	gmp_randinit_default(state);

	for(i = 0; i < RUNS; i++)
	{
		test_rnd(state);
		if((i % (RUNS/10)) == 0)
		{
			fprintf(stderr, "Step %d of %d\n", i, RUNS);
		}
	}
	return 0;
}
