/* Note that this define is required for syscalls to work. */
#define _GNU_SOURCE

#include <unistd.h>
#include <sys/syscall.h>

#include <mpfr.h>
#include <stdlib.h>
#include <stdio.h>

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

#define RND_DIR MPFR_RNDN
#define GOLD_PREC 500


int main(int argc, char *argv[])
{
	int i, b;
	gmp_randstate_t state;
	mpfr_t g, a;
	unsigned long int seed;

	gmp_randinit_default(state);

	while(getrandom(&seed, sizeof(seed), 0) != sizeof(seed));

	printf("%ld\n", seed);

	gmp_randseed_ui(state, seed);

	if(argc != 2)
	{
		printf("Please specify bits.\n");
		return 1;
	}

	b = atoi(argv[1]);
	mpfr_init2(a, b);
	mpfr_init2(g, GOLD_PREC);

	for(i = 0; i < 10; i++)
	{
		mpfr_urandom(g, state, RND_DIR);
		mpfr_set(a, g, RND_DIR);
		mpfr_printf("%Rf ", a);
	}

	mpfr_printf("\n");

	mpfr_clear(a);
	mpfr_clear(g);

	return 0;
}
