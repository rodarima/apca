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

#define RUNS 100000
#define SMALL_PREC 10

void test_rounding(gmp_randstate_t state)
{
	int i, exp;
	mpfr_t gold, small, err_abs, err_rel;

	mpfr_init2(gold, GOLD_PREC);
	mpfr_init2(small, SMALL_PREC);
	mpfr_init2(err_abs, GOLD_PREC);
	mpfr_init2(err_rel, GOLD_PREC);

	mpfr_printf("g exp err_abs err_rel\n");

	for(i=0; i<RUNS; i++)
	{
		mpfr_urandom(gold, state, RND_DIR);
		mpfr_set(small, gold, RND_DIR);
		mpfr_sub(err_abs, gold, small, RND_DIR);
		mpfr_abs(err_rel, err_abs, RND_DIR);
		mpfr_div(err_rel, err_rel, gold, RND_DIR);
		exp = mpfr_get_exp(small);
		mpfr_printf("%Re %d %Re %Re\n", gold, exp, err_abs, err_rel);
	}

	mpfr_clear(gold);
	mpfr_clear(small);
	mpfr_clear(err_abs);
	mpfr_clear(err_rel);
}

int main(int argc, char *argv[])
{
	gmp_randstate_t state;
	unsigned long int seed;

	gmp_randinit_default(state);

	while(getrandom(&seed, sizeof(seed), 0) != sizeof(seed));

	gmp_randseed_ui(state, seed);

	test_rounding(state);

	return 0;
}
