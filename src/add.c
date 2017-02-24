#include <stdio.h>
#include <mpfr.h>

#define RND MPFR_RNDZ

void test_prec(mpfr_prec_t prec)
{
	int i;
	mpfr_t gold, a, b, c, err;
	float _a, _b, _c, _err;
	
	mpfr_init2(gold, prec);
	mpfr_init2(err, prec);
	mpfr_init2(a, prec);
	mpfr_init2(b, prec);
	mpfr_init2(c, prec);

	_a = 0.123456789;
	_b = 0.876543211;
	_c = _a + _b;

	mpfr_set_flt(a, _a, RND);
	mpfr_set_flt(b, _b, RND);
	mpfr_add(c, a, b, RND);

	mpfr_printf("%3d bits --------------------\n", prec);
	mpfr_printf("a:%.10Rf _a:%.10f\n", a, _a);
	mpfr_printf("b:%.10Rf _b:%.10f\n", b, _b);
	mpfr_printf("c:%.10Rf _c:%.10f\n", c, _c);

	mpfr_clear(gold);
	mpfr_clear(err);
	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);
}

void test(mpfr_prec_t prec, mpfr_rnd_t rnd, char *rnd_str)
{
	int i;
	mpfr_t gold, a, b, c, d, err;
	float _a, _b, _c, _d, _err;
	
	mpfr_init2(gold, prec);
	mpfr_init2(err, prec);
	mpfr_init2(a, prec);
	mpfr_init2(b, prec);
	mpfr_init2(c, prec);
	mpfr_init2(d, prec);

	_a = 0.123456789;
	_b = 0.876543211;
	_c = _a + _b;
	_d = _a - _b;

	mpfr_set_flt(a, _a, rnd);
	mpfr_set_flt(b, _b, rnd);
	mpfr_add(c, a, b, rnd);
	mpfr_sub(d, a, b, rnd);

	mpfr_printf("------ %3d bits, rnd mode: %s ------\n", prec, rnd_str);
	mpfr_printf("a:%.10Rf _a:%.10f\n", a, _a);
	mpfr_printf("b:%.10Rf _b:%.10f\n", b, _b);
	mpfr_printf("c:%.10Rf _c:%.10f\n", c, _c);
	mpfr_printf("d:%.10Rf _d:%.10f\n", d, _d);

	mpfr_clear(gold);
	mpfr_clear(err);
	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);
	mpfr_clear(d);
}

struct rnd_t
{
	mpfr_rnd_t rnd;
	char *rnd_str;
};

#define N_RND_MODES 5
struct rnd_t rnd_modes[N_RND_MODES] = {
	{MPFR_RNDN, "Round to Nearest"},
	{MPFR_RNDZ, "Round toward zero"},
	{MPFR_RNDU, "Round toward plus infinity"},
	{MPFR_RNDD, "Round toward minus infinity"},
	{MPFR_RNDA, "Round away from zero"},
};


int main(int argc, char *argv[])
{
	int i,j;

	for(i = 0; i < N_RND_MODES; i++)
	{
		for(j = 24; j <= 24; j++)
		{
			test(j, rnd_modes[i].rnd, rnd_modes[i].rnd_str);
		}
	}

	return 0;
}
