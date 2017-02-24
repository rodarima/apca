#ifndef _REDUCTION_H_
#define _REDUCTION_H_

#include <mpfr.h>

void mp_tred2(mpfr_t **a, int n, mpfr_t *d, mpfr_t *e, mpfr_prec_t prec, mpfr_rnd_t rnd);

void tred2(float **a, int n, float *d, float *e);

#endif /* _REDUCTION_H_ */

