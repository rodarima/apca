#include <stdio.h>
#include <mpfr.h>

#define PREC_BIG 500

void test_bits(int n, mpfr_t expected)
{
	mpfr_t a, b, c;
	mpfr_t err;

	mpfr_init2(a, n);
	mpfr_init2(b, n);
	mpfr_init2(c, n);
	mpfr_init2(c, PREC_BIG);

	mpfr_set_d(a, 30.0, MPFR_RNDD);
	mpfr_set_d(b,  5.0, MPFR_RNDD);

	mpfr_add(c, a, b, MPFR_RNDD);
	printf("n=%5d: ", n);
	mpfr_out_str(stdout, 10, 0, c, MPFR_RNDD);
	putchar('\n');

	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);
}

int main(int argc, char *argv[])
{
	int i;
	mpfr_t gold, a, b, c, err;
	
	mpfr_init2(gold, PREC_BIG);
	mpfr_init2(err, PREC_BIG);
	mpfr_init2(a, PREC_BIG);
	mpfr_init2(b, PREC_BIG);
	mpfr_init2(c, PREC_BIG);

	//mpfr_set_d(a, 301.614398744345975432416, MPFR_RNDD);
	mpfr_set_str(a,"302.614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416614398744345975432416", 10, MPFR_RNDD);
	mpfr_set_d(b,  50.0, MPFR_RNDD);
	mpfr_add(gold, a, b, MPFR_RNDD);
	mpfr_set_d(c,  0.0, MPFR_RNDD);

	for(i = 100; i >= 3; i--)
	{
		mpfr_prec_round(a, i, MPFR_RNDD);
		mpfr_prec_round(b, i, MPFR_RNDD);
		mpfr_set_prec(c, i);

		mpfr_add(c, a, b, MPFR_RNDD);
		mpfr_sub(err, gold, c, MPFR_RNDD);
		mpfr_abs(err, err, MPFR_RNDD);
		mpfr_div(err, err, gold, MPFR_RNDD);
		mpfr_printf("%d %.10Re\n", i, err);
	}
	
	mpfr_set_prec(a, PREC_BIG);
	mpfr_set_prec(b, PREC_BIG);
	mpfr_set_prec(c, PREC_BIG);

	mpfr_clear(gold);
	mpfr_clear(err);
	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);

	return 0;
}
