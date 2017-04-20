#include "reduction.h"
#include "mpvector.h"

#include <math.h>
#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N 10
#define PREC 24 /* this machine is NOT IEEE 754 compliant */
#define RND_DIR MPFR_RNDN

#define GOLD_PREC 500
#define PREC_MIN 10
#define PREC_MAX 100
#define PREC_SUM 10
#define RUNS 10

#define CSV_HEADER "n bits err time time_gold"


#include <linux/version.h>

#if LINUX_VERSION_CODE <= KERNEL_VERSION(3,7,0)

#define RANDOM_DEV "/dev/urandom"

/* getrandom is not yet available on 2.6 Linux.
 * This is only a mock. Flags are ignored. */
int getrandom(void *buf, size_t buflen, unsigned int flags)
{
	FILE *f;
	size_t r;

	f = fopen(RANDOM_DEV, "r");

	if(!f) return -1;

	r = fread(buf, 1, buflen, f);

	fclose(f);

	return r;
}

#else
#include <linux/random.h>
#endif


FILE *outfile;

double sec()
{
	struct timeval t;
	if (gettimeofday(&t, 0) < 0 ) return 0.0;
	return (t.tv_usec / 1000000.0 + t.tv_sec);
}

double hh(struct mptridiag_t *td, mpfr_prec_t *prec)
{
	double t = sec();
	mp_tred2(td->A, td->n, td->diag, td->offdiag, prec, td->rnd);
	t = sec() - t;
	return t;
}

void ql(struct mptridiag_t *td, mpfr_prec_t *prec)
{
	//int i;
	//for(i=9; i < 19; i++)
	//	printf("prec[%d]=%d\n", i, prec[i]);
	mp_tqli(td->diag, td->offdiag, td->n, td->A, prec, td->rnd);
}

void set_prec(mpfr_prec_t *prec_vec, int n, mpfr_prec_t prec)
{
	int i;

	for(i=0; i<n; i++)
	{
		prec_vec[i] = prec;
	}
}

void errcond(
	struct mptridiag_t *gold, struct mptridiag_t *target,
	mpfr_prec_t *prec, int var, double time_gold)
{
	double t, e;
	int bits, n;
	mpfr_t err, err1, norm2_err, norm2_rel, cond;
	mpfr_prec_t *ptri, *phh;

	ptri = prec;
	phh = ptri + VARS_MPTRIDIAG;

	t = hh(target, phh);
	n = gold->n;
	bits = phh[var];

	mpfr_init2(err1, gold->prec);
	mpfr_init2(err, gold->prec);
	mpfr_init2(norm2_err, gold->prec);
	mpfr_init2(norm2_rel, gold->prec);
	mpfr_init2(cond, gold->prec);

	mpfr_sub(err, gold->diag[1], target->diag[1], gold->rnd);
	mpfr_abs(err1, err, gold->rnd);

	mpvector_diff(norm2_err, gold->diag, target->diag, gold->n, gold->prec, gold->rnd);
	mpvector_norm2(norm2_rel, target->diag, gold->n, gold->prec, gold->rnd);

	mpfr_div(norm2_rel, norm2_err, norm2_rel, gold->rnd);

	//ql(target, pql); //XXX gold prec?
	//mp_eigcond(cond, target->diag, target->n, gold->prec, gold->rnd);

	//mpfr_printf("errcond:%d bits, error %.20Rg\n", target->prec, err);


	e = mpfr_get_d(norm2_err, gold->rnd);

	fprintf(outfile, "%d %d %f %e %e\n", n, bits, bits+log2(e), t, time_gold);

//	mpfr_printf("%d\t%d\t%d\t%.20Rg\t%.20Rg\t%.20Rg\t%.20Rg\n",
//		phh[var], var, gold->n, err1, err, norm2_err, norm2_rel);

	mpfr_clear(cond);
	mpfr_clear(norm2_err);
	mpfr_clear(norm2_rel);
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

#define SIZE_PREC_SET 10
#define SIZE_CONF_VAR 9
#define VARS_TOTAL (VARS_MPTRIDIAG + VARS_TRED2 + VARS_TQLI)

void test_rnd(gmp_randstate_t state, int n, int bits)
{
	int i;
	struct mptridiag_t *gold, *target, *ref;
	double t;

	mpfr_prec_t prec_vec[VARS_TOTAL];

	mpfr_prec_t *ptri, *phh;
	ptri = prec_vec;
	phh = ptri + VARS_MPTRIDIAG;

	set_prec(ptri, VARS_MPTRIDIAG, GOLD_PREC);
	set_prec(phh, VARS_TRED2, GOLD_PREC);

	gold = mptridiag_init(n, ptri, RND_DIR);
	ref  = mptridiag_init(n, ptri, RND_DIR);
	//load_matrix(gold, AA);
	rand_mpmat(ref, state);
	copy_mpmat(gold, ref);

	t = hh(gold, phh);
	//fprintf(stderr, "Gold HH took %f seconds.\n", t);
	//fprintf(stderr, "%s\n", CSV_HEADER);
	//printf("Gold diagonal:\n");
	//mpvector_print(gold->diag, gold->n);

	if (bits != 0)
	{
		set_prec(phh, VARS_TRED2, bits);
		target = mptridiag_init(n, ptri, RND_DIR);
		copy_mpmat(target, ref);
		errcond(gold, target, prec_vec, 0, t);
		mptridiag_free(target);
	}
	else
	{

		for(i = PREC_MIN; i <= PREC_MAX; i+=PREC_SUM)
		{
			set_prec(phh, VARS_TRED2, i);
			target = mptridiag_init(n, ptri, RND_DIR);
			//load_matrix(target, AA);
			copy_mpmat(target, ref);
			errcond(gold, target, prec_vec, 0, t);
			//printf("Target diagonal:\n");
			//mpvector_print(target->diag, target->n);
			mptridiag_free(target);
			//phh[j] = GOLD_PREC;
		}
	}

	mptridiag_free(gold);

}

#define N_MIN 500
#define N_MAX 5000
#define N_MUL 1.5
#define N_SUM 5

void test_size(gmp_randstate_t st, int n, int bits)
{
	test_rnd(st, n, bits);
	fflush(outfile);
}

void usage(int argc, char *argv[])
{
	printf("Usage: %s n [bits [runs]]\n", argv[0]);
	printf("\n");
	printf("\n");
	printf("This program computes the Householder tridiagonalization for a n-by-n matrix.\n");
	printf("First a gold value is obtained using a precision of %d bits.\n", GOLD_PREC);
	printf("Then new computations are performed from %d to %d bits, each %d.\n",
		PREC_MIN, PREC_MAX, PREC_SUM);
	printf("If 'bits' is specified and is not 0, only one computation is performed\n");
	printf("with the specified bit-width.\n");
	printf("\n");
	printf("If 'runs' is specified, the experiment is repeated this amount of times.\n");
	printf("\n");
	printf("The results are measured and printed out in CSV with the following format:\n");
	printf("%s\n", CSV_HEADER);
	printf("n         : The size of the n-by-n matrix.\n");
	printf("bits      : The number of bits of precision in the computation.\n");
	printf("err       : The norm-2 error measured ONLY on the diagonal.\n");
	printf("time      : The time in seconds of the computation.\n");
	printf("time_gold : The time in seconds computing the gold value.\n");
	printf("\n");
	printf("Note: The header is also included in the first line.\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	int n, bits = 0, runs = 1, i;
	char *filename = NULL;
	gmp_randstate_t state;
	unsigned long int seed;


	if(argc < 2 || argc > 5) usage(argc, argv);

	n = (int) atoi(argv[1]);

	if(argc > 2) bits = (int) atoi(argv[2]);

	if(argc > 3) runs = (int) atoi(argv[3]);

	if(argc > 4) filename = argv[4];

	if(n <= 0) usage(argc, argv);
	if(bits < 0) usage(argc, argv);
	if(runs <= 0) usage(argc, argv);

	if(filename)
	{
		outfile = fopen(filename, "w");
		if(!outfile)
		{
			perror("fopen");
			exit(1);
		}
	}
	else
	{
		outfile = stdout;
	}

	gmp_randinit_default(state);

	while(getrandom(&seed, sizeof(seed), 0) != sizeof(seed));

	gmp_randseed_ui(state, seed);

	fprintf(outfile, "%s\n", CSV_HEADER);
	for(i = 0; i < runs; i++)
	{
		//fprintf(stderr, "Step %d of %d\n", i, runs);
		test_size(state, n, bits);
	}

	fclose(outfile);
	return 0;
}
