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

#define RUNS 1e6
#define SMALL_PREC 10

void test_rounding(gmp_randstate_t state)
{
	int i, exp;
	mpfr_t ga, gb, gc;
	mpfr_t ra, rb, rc;
	mpfr_t err_a, err_b, err_c, err_abs, err_rel;

	mpfr_init2(ga, GOLD_PREC);
	mpfr_init2(gb, GOLD_PREC);
	mpfr_init2(gc, GOLD_PREC);

	mpfr_init2(ra, SMALL_PREC);
	mpfr_init2(rb, SMALL_PREC);
	mpfr_init2(rc, SMALL_PREC);

	mpfr_init2(err_a, GOLD_PREC);
	mpfr_init2(err_b, GOLD_PREC);
	mpfr_init2(err_c, GOLD_PREC);
	mpfr_init2(err_abs, GOLD_PREC);
	mpfr_init2(err_rel, GOLD_PREC);

	mpfr_printf("ga gb gc exp err_a, err_b, err_c, err_abs err_rel\n");

	for(i = 0; i < RUNS; i++)
	{
		mpfr_urandom(ga, state, RND_DIR);
		//mpfr_urandom(gb, state, RND_DIR);
		mpfr_set_d(gb, 1.0, RND_DIR);
		mpfr_add(gc, ga, gb, RND_DIR);

		mpfr_set(ra, ga, RND_DIR);
		mpfr_set(rb, gb, RND_DIR);
		mpfr_add(rc, ra, rb, RND_DIR);

		mpfr_sub(err_a, ga, ra, RND_DIR);
		mpfr_sub(err_b, gb, rb, RND_DIR);
		mpfr_sub(err_c, gc, rc, RND_DIR);
		mpfr_abs(err_abs, err_c, RND_DIR);
		mpfr_div(err_rel, err_abs, gc, RND_DIR);

		exp = mpfr_get_exp(rc);

		mpfr_printf("%.20Re %.20Re %.20Re %d %.20Re %.20Re %.20Re %.20Re %.20Re\n",
			ga, gb, gc, exp, err_a, err_b, err_c, err_abs, err_rel);
	}

	mpfr_clear(ga);
	mpfr_clear(gb);
	mpfr_clear(gc);

	mpfr_clear(ra);
	mpfr_clear(rb);
	mpfr_clear(rc);

	mpfr_clear(err_a);
	mpfr_clear(err_b);
	mpfr_clear(err_c);
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
