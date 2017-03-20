#include "pca.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <mpfr.h>

#define RND_DIR MPFR_RNDN

#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )


/**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */

void mp_tred2(mpfr_t **a, int n, mpfr_t *d, mpfr_t *e, mpfr_prec_t *prec, mpfr_rnd_t rnd)
/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.  */
{
	int l, k, j, i;
	mpfr_t scale, hh, h, g, f, tmp1;

	mpfr_init2(scale, prec[3]);
	mpfr_init2(hh, prec[4]);
	mpfr_init2(h, prec[5]);
	mpfr_init2(g, prec[6]);
	mpfr_init2(f, prec[7]);
	mpfr_init2(tmp1, prec[8]);

	for (i = n; i >= 2; i--)
	{
		l = i - 1;
		mpfr_set_d(scale, 0.0, rnd);
		mpfr_set_d(h, 0.0, rnd);
		//h = scale = 0.0;
		if (l > 1)
		{
			for (k = 1; k <= l; k++)
				mpfr_abs(tmp1, a[i][k], rnd);
				mpfr_add(scale, scale, tmp1, rnd);
				//scale += fabs(a[i][k]);
			mpfr_fprintf(stderr, "scale = %Rf\n", scale);
			if (mpfr_zero_p(scale))
			//if (scale == 0.0)
				mpfr_set(e[i], a[i][l], rnd);
				//e[i] = a[i][l];
			else
			{
				for (k = 1; k <= l; k++)
				{
					mpfr_div(a[i][k], a[i][k], scale, rnd);
					//a[i][k] /= scale;
					mpfr_mul(tmp1, a[i][k], a[i][k], rnd);
					mpfr_add(h, h, tmp1, rnd);
					//h += a[i][k] * a[i][k];
				}
				mpfr_set(f, a[i][l], rnd);
				//f = a[i][l];
				mpfr_sqrt(tmp1, h, rnd);
				if(mpfr_sgn(f) > 0)
				{
					mpfr_neg(tmp1, tmp1, rnd);
				}
				mpfr_set(g, tmp1, rnd);
				//g = f>0 ? -sqrt(h) : sqrt(h);
				mpfr_mul(e[i], scale, g, rnd);
				//e[i] = scale * g;
				mpfr_mul(tmp1, f, g, rnd);
				mpfr_sub(h, h, tmp1, rnd);
				//h -= f * g;
				mpfr_sub(a[i][l], f, g, rnd);
				//a[i][l] = f - g;
				mpfr_set_d(f, 0.0, rnd);
				//f = 0.0;
				for (j = 1; j <= l; j++)
				{
					mpfr_div(a[j][i], a[i][j], h, rnd);
					//a[j][i] = a[i][j]/h;
					mpfr_set_d(g, 0.0, rnd);
					//g = 0.0;
					for (k = 1; k <= j; k++)
						mpfr_fma(g, a[j][k], a[i][k], g, rnd);
						//g += a[j][k] * a[i][k];
					for (k = j+1; k <= l; k++)
						mpfr_fma(g, a[k][j], a[i][k], g, rnd);
						//g += a[k][j] * a[i][k];
					mpfr_div(e[j], g, h, rnd);
					//e[j] = g / h;
					mpfr_fma(f, e[j], a[i][j], f, rnd);
					//f += e[j] * a[i][j];
				}
				mpfr_add(tmp1, h, h, rnd);
				mpfr_div(hh, f, tmp1, rnd);
				//hh = f / (h + h);
				for (j = 1; j <= l; j++)
				{
					mpfr_set(f, a[i][j], rnd);
					//f = a[i][j];
					mpfr_mul(tmp1, hh, f, rnd);
					mpfr_sub(g, e[j], tmp1, rnd);
					mpfr_set(e[j], g, rnd);
					//e[j] = g = e[j] - hh * f;
					for (k = 1; k <= j; k++)
					{
						mpfr_mul(tmp1, g, a[i][k], rnd);
						mpfr_fma(tmp1, f, e[k], tmp1, rnd);
						mpfr_sub(a[j][k], a[j][k], tmp1, rnd);
						//a[j][k] -= (f * e[k] + g * a[i][k]);
					}
				}
			}
		}
		else
			mpfr_set(e[i], a[i][l], rnd);
			//e[i] = a[i][l];
		mpfr_set(d[i], h, rnd);
		//d[i] = h;
	}
	mpfr_set_d(d[1], 0.0, rnd);
	//d[1] = 0.0;
	mpfr_set_d(e[1], 0.0, rnd);
	//e[1] = 0.0;
	for (i = 1; i <= n; i++)
	{
		l = i - 1;
		if (d[i])
		{
			for (j = 1; j <= l; j++)
			{
				mpfr_set_d(g, 0.0, rnd);
				//g = 0.0;
				for (k = 1; k <= l; k++)
				{
					mpfr_fma(g, a[i][k], a[k][j], g, rnd);
					//g += a[i][k] * a[k][j];
				}
				for (k = 1; k <= l; k++)
				{
					mpfr_mul(tmp1, g, a[k][i], rnd);
					mpfr_sub(a[k][j], a[k][j], tmp1, rnd);
					//a[k][j] -= g * a[k][i];
				}
			}
		}
		mpfr_set(d[i], a[i][i], rnd);
		//d[i] = a[i][i];
		mpfr_set_d(a[i][i], 1.0, rnd);
		//a[i][i] = 1.0;
		for (j = 1; j <= l; j++)
		{
			mpfr_set_d(a[i][j], 0.0, rnd);
			mpfr_set_d(a[j][i], 0.0, rnd);
			//a[j][i] = a[i][j] = 0.0;
		}
	}

	mpfr_clear(scale);
	mpfr_clear(hh);
	mpfr_clear(h);
	mpfr_clear(g);
	mpfr_clear(f);
	mpfr_clear(tmp1);
}

/**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */

void tred2(a, n, d, e)
float **a, *d, *e;
/* float **a, d[], e[]; */
int n;
/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.  */
{
	int l, k, j, i;
	float scale, hh, h, g, f;

	for (i = n; i >= 2; i--)
	{
		l = i - 1;
		h = scale = 0.0;
		if (l > 1)
		{
			for (k = 1; k <= l; k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i] = a[i][l];
			else
			{
				for (k = 1; k <= l; k++)
				{
					a[i][k] /= scale;
					h += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = f>0 ? -sqrt(h) : sqrt(h);
				e[i] = scale * g;
				h -= f * g;
				a[i][l] = f - g;
				f = 0.0;
				for (j = 1; j <= l; j++)
				{
					a[j][i] = a[i][j]/h;
					g = 0.0;
					for (k = 1; k <= j; k++)
						g += a[j][k] * a[i][k];
					for (k = j+1; k <= l; k++)
						g += a[k][j] * a[i][k];
					e[j] = g / h;
					f += e[j] * a[i][j];
				}
				hh = f / (h + h);
				for (j = 1; j <= l; j++)
				{
					f = a[i][j];
					e[j] = g = e[j] - hh * f;
					for (k = 1; k <= j; k++)
						a[j][k] -= (f * e[k] + g * a[i][k]);
				}
			}
		}
		else
			e[i] = a[i][l];
		d[i] = h;
	}
	d[1] = 0.0;
	e[1] = 0.0;
	for (i = 1; i <= n; i++)
	{
		l = i - 1;
		if (d[i])
		{
			for (j = 1; j <= l; j++)
			{
				g = 0.0;
				for (k = 1; k <= l; k++)
					g += a[i][k] * a[k][j];
				for (k = 1; k <= l; k++)
					a[k][j] -= g * a[k][i];
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for (j = 1; j <= l; j++)
			a[j][i] = a[i][j] = 0.0;
	}
}

/**  Tridiagonal QL algorithm -- Implicit  **********************/
#if 0
void tqli(d, e, n, z)
	float d[], e[], **z;
	int n;
{
	int m, l, iter, i, k;
	float s, r, p, g, f, dd, c, b;
	void erhand();

	for (i = 2; i <= n; i++)
		e[i-1] = e[i];
	e[n] = 0.0;
	for (l = 1; l <= n; l++)
	{
		iter = 0;
		do
		{
			for (m = l; m <= n-1; m++)
			{
				dd = fabs(d[m]) + fabs(d[m+1]);
				if (fabs(e[m]) + dd == dd) break;
			}
			if (m != l)
			{
				if (iter++ == 30) erhand("No convergence in TLQI.");
				g = (d[l+1] - d[l]) / (2.0 * e[l]);
				r = sqrt((g * g) + 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m-1; i >= l; i--)
				{
					f = s * e[i];
					b = c * e[i];
					if (fabs(f) >= fabs(g))
					{
						c = g / f;
						r = sqrt((c * c) + 1.0);
						e[i+1] = f * r;
						c *= (s = 1.0/r);
					}
					else
					{
						s = f / g;
						r = sqrt((s * s) + 1.0);
						e[i+1] = g * r;
						s *= (c = 1.0/r);
					}
					g = d[i+1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					p = s * r;
					d[i+1] = g + p;
					g = c * r - b;
					for (k = 1; k <= n; k++)
					{
						f = z[k][i+1];
						z[k][i+1] = s * z[k][i] + c * f;
						z[k][i] = c * z[k][i] - s * f;
					}
				}
				d[l] = d[l] - p;
				e[l] = g;
				e[m] = 0.0;
			}
			//printf("TLQI: Iteration %d, g=%f, m=%d, l=%d\n", iter, g, m, l);
		}
		while (m != l);
	}
}
#endif
//***************************************************************************
//    3.0    3.0    3.0    3.0    3.0    3.0   35.0   45.0
//   53.0   55.0   58.0  113.0  113.0   86.0   67.0   90.0
//    3.5    3.5    4.0    4.0    4.5    4.5   46.0   59.0
//   63.0   58.0   58.0  125.0  126.0  110.0   78.0   97.0
//    4.0    4.0    4.5    4.5    5.0    5.0   48.0   60.0
//   68.0   65.0   65.0  123.0  123.0  117.0   87.0  108.0
//    5.0    5.0    5.0    5.5    5.5    5.5   46.0   63.0
//   70.0   64.0   63.0  116.0  119.0  115.0   97.0  112.0
//    6.0    6.0    6.0    6.0    6.5    6.5   51.0   69.0
//   77.0   70.0   71.0  120.0  122.0  122.0   96.0  123.0
//   11.0   11.0   11.0   11.0   11.0   11.0   64.0   75.0
//   81.0   79.0   79.0  112.0  114.0  113.0   98.0  115.0
//   20.0   20.0   20.0   20.0   20.0   20.0   76.0   86.0
//   93.0   92.0   91.0  104.0  104.5  107.0   97.5  104.0
//   30.0   30.0   30.0   30.0   30.1   30.2   84.0   96.0
//   98.0   99.0   96.0  101.0  102.0   99.0   94.0   99.0
//   30.0   33.4   36.8   40.0   43.0   45.6  100.0  106.0
//  106.0  108.0  101.0   99.0   98.0   99.0   95.0   95.0
//   42.0   44.0   46.0   48.0   50.0   51.0  109.0  111.0
//  110.0  110.0  103.0   95.5   95.5   95.0   92.5   92.0
//   60.0   61.7   63.5   65.5   67.3   69.2  122.0  124.0
//  124.0  121.0  103.0   93.2   92.5   92.2   90.0   90.8
//   70.0   70.1   70.2   70.3   70.4   70.5  137.0  132.0
//  134.0  128.0  101.0   91.7   90.2   88.8   87.3   85.8
//   78.0   77.6   77.2   76.8   76.4   76.0  167.0  159.0
//  152.0  144.0  103.0   89.8   87.7   85.7   83.7   81.8
//   98.9   97.8   96.7   95.5   94.3   93.2  183.0  172.0
//  162.0  152.0  102.0   87.5   85.3   83.3   81.3   79.3
//  160.0  157.0  155.0  152.0  149.0  147.0  186.0  175.0
//  165.0  156.0  120.0   87.0   84.9   82.8   80.8   79.0
//  272.0  266.0  260.0  254.0  248.0  242.0  192.0  182.0
//  170.0  159.0  131.0   88.0   85.8   83.7   81.6   79.6
//  382.0  372.0  362.0  352.0  343.0  333.0  205.0  192.0
//  178.0  166.0  138.0   86.2   84.0   82.0   79.8   77.5
//  770.0  740.0  710.0  680.0  650.0  618.0  226.0  207.0
//  195.0  180.0  160.0   82.9   80.2   77.7   75.2   72.7
//***************************************************************************

