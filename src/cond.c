#include "pca.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <mpfr.h>
#include "reduction.h"

/* Compute the condition number from the eigenvalues.
 * The eigenvalues should be sorted starting at the lowest. */
mpfr_t mp_eigcond(mpfr_t *ev, int n, mpfr_prec_t prec, mpfr_rnd_t rnd)
{
	mpfr_t c, tmp1;
	mpfr_t *low, *big;

	low = &ev[1];
	big = &ev[n];

	mpfr_init2(c, prec);
	mpfr_div(c, *big, *low, rnd);
	mpfr_abs(c, c, rnd);

	return c;
}
