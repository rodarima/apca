/* Note that this define is required for syscalls to work. */
#define _GNU_SOURCE

#include <unistd.h>
#include <sys/syscall.h>

#include <mpfr.h>
#include <stdlib.h>
#include <stdio.h>

#include <linux/version.h>

#if (1 || (LINUX_VERSION_CODE <= KERNEL_VERSION(3,7,0)))

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
/* Note that this define is required for syscalls to work. */
#define _GNU_SOURCE

#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>
#endif

#define RND_DIR MPFR_RNDN
#define GOLD_PREC 500

#define RUNS 1e5
#define SMALL_PREC 10
#define STEPS 1

void test_rounding(gmp_randstate_t state, int runs, int steps, int bits)
{
	int i, j, exp;
	mpfr_t grnd, gsum;
	mpfr_t rrnd, rsum;
	mpfr_t err, err_abs, err_rel, err2;

	mpfr_init2(grnd, GOLD_PREC);
	mpfr_init2(gsum, GOLD_PREC);

	mpfr_init2(rrnd, bits);
	mpfr_init2(rsum, bits);

	mpfr_init2(err, GOLD_PREC);
	mpfr_init2(err2, GOLD_PREC);
	mpfr_init2(err_abs, GOLD_PREC);
	mpfr_init2(err_rel, GOLD_PREC);


	mpfr_printf("gsum exp err err_abs err_rel err2\n");

	for(i = 0; i < runs; i++)
	{
		mpfr_set_d(err2, 0.0, RND_DIR);
		mpfr_set_d(gsum, 0.0, RND_DIR);
		mpfr_set_d(rsum, 0.0, RND_DIR);

		for(j = 0; j < steps; j++)
		{
			mpfr_urandom(grnd, state, RND_DIR);
			//mpfr_set(rrnd, grnd, RND_DIR);

			mpfr_add(rsum, rsum, grnd, RND_DIR);
			mpfr_add(gsum, gsum, grnd, RND_DIR);

			mpfr_sub(err, gsum, rsum, RND_DIR);
			mpfr_mul(err, err, err, RND_DIR);
			mpfr_add(err2, err2, err, RND_DIR);
		}

		mpfr_sub(err, gsum, rsum, RND_DIR);
		mpfr_abs(err_abs, err, RND_DIR);
		mpfr_div(err_rel, err_abs, gsum, RND_DIR);

		exp = mpfr_get_exp(rsum);

		mpfr_printf("%.20Re %d %.20Re %.20Re %.20Re %.20Re\n",
			gsum, exp, err, err_abs, err_rel, err2);
	}

	mpfr_clear(grnd);
	mpfr_clear(gsum);

	mpfr_clear(rrnd);
	mpfr_clear(rsum);

	mpfr_clear(err);
	mpfr_clear(err_abs);
	mpfr_clear(err_rel);
}

void usage()
{
	printf("exp8c: Computes the accumulated error on consecutive sums.\n");
	printf("Usage: ./exp8c RUNS STEPS BITS\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	gmp_randstate_t state;
	unsigned long int seed;
	int runs, steps, bits;

	if(argc != 4) usage();

	runs = atoi(argv[1]);
	steps = atoi(argv[2]);
	bits = atoi(argv[3]);


	gmp_randinit_default(state);

	while(getrandom(&seed, sizeof(seed), 0) != sizeof(seed));

	gmp_randseed_ui(state, seed);

	test_rounding(state, runs, steps, bits);

	return 0;
}
