#include "mpvector.h"

int main()
{
	int prec = 100;
	int n = 10, m = 20;
	mpfr_t *vec, **mat;

	vec = vector_init(n, prec);
	mat = matrix_init(n, m, prec);

	vector_free(vec, n);
	matrix_free(mat, n, m);
	return 0;
}
