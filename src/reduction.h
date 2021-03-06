#ifndef _REDUCTION_H_
#define _REDUCTION_H_

#include <mpfr.h>

#define VARS_TRED2 6

void mp_tred2(
	mpfr_t **a, int n, mpfr_t *d, mpfr_t *e,
	mpfr_prec_t prec[VARS_TRED2], mpfr_rnd_t rnd);

void tred2(float **a, int n, float *d, float *e);

#define VARS_TQLI 10

void mp_tqli(
	mpfr_t *d, mpfr_t *e, int n, mpfr_t **z,
	mpfr_prec_t prec[VARS_TQLI], mpfr_rnd_t rnd);

void tqli(float d[], float e[], int n, float **z);

void mp_eigcond(mpfr_t c, mpfr_t *ev, int n, mpfr_prec_t prec, mpfr_rnd_t rnd);
#endif /* _REDUCTION_H_ */

