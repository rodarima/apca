/*
    This file is part of MPFR C++.

    MPFR C++: Multi-precision floating point number class for C++. 
    Based on MPFR library:    http://mpfr.org

    Project homepage:    http://www.holoborodko.com/pavel/mpfr
    Contact e-mail:      pavel@holoborodko.com

    Copyright (c) 2008-2013 Pavel Holoborodko

    Contributors:
    Dmitriy Gubanov, Konstantin Holoborodko, Brian Gladman, 
    Helmut Jarausch, Fokko Beekhof, Ulrich Mutze, Heinz van Saanen, 
    Pere Constans, Peter van Hoof, Gael Guennebaud, Tsai Chia Cheng, 
    Alexei Zubanov, Jauhien Piatlicki, Victor Berger, John Westwood,
    Petr Aleksandrov, Orion Poplawski, Charles Karney, Arash Partow.

    Licensing:
    (A) MPFR C++ is under GNU General Public License ("GPL").
    
    (B) Non-free licenses may also be purchased from the author, for users who 
        do not want their programs protected by the GPL.

        The non-free licenses are for users that wish to use MPFR C++ in 
        their products but are unwilling to release their software 
        under the GPL (which would require them to release source code 
        and allow free redistribution).

        Such users can purchase an unlimited-use license from the author.
        Contact us for more details.
    
    GNU General Public License ("GPL") copyright permissions statement:
    **************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <mpreal.h>


/**  Allocation of vector storage  ***********************************/

mpreal vector_init(int n)
{

    float *v;

    v = (float *) malloc ((unsigned) n*sizeof(float));
    if (!v) erhand("Allocation failure in vector().");
    return v-1;

}

/**  Allocation of float matrix storage  *****************************/

float **matrix_init(n,m)
int n, m;
/* Allocate a float matrix with range [1..n][1..m]. */
{
    int i;
    float **mat;

    /* Allocate pointers to rows. */
    mat = (float **) malloc((unsigned) (n)*sizeof(float*));
    if (!mat) erhand("Allocation failure 1 in matrix().");
    mat -= 1;

    /* Allocate rows and set pointers to them. */
    for (i = 1; i <= n; i++)
        {
        mat[i] = (float *) malloc((unsigned) (m)*sizeof(float));
        if (!mat[i]) erhand("Allocation failure 2 in matrix().");
        mat[i] -= 1;
        }

     /* Return pointer to array of pointers to rows. */
     return mat;

}

/**  Deallocate vector storage  *********************************/

void vector_free(v,n)
float *v;
int n;
/* Free a float vector allocated by vector(). */
{
   free((char*) (v+1));
}

/**  Deallocate float matrix storage  ***************************/

void matrix_free(mat,n,m)
float **mat;
int n, m;
/* Free a float matrix allocated by matrix(). */
{
   int i;

   for (i = n; i >= 1; i--)
       {
       free ((char*) (mat[i]+1));
       }
   free ((char*) (mat+1));
}

void matrix_print(float **mat, int n, int m)
{
	int i, j;

	for(i=1; i<=n; i++)
	{
		for(j=1; j<=m; j++)
		{
			printf("%15.12f ", mat[i][j]);
		}
		printf("\n");
	}
}

void vector_print(float *vec, int n)
{
	int i;

	for(i=1; i<=n; i++)
	{
		printf("%.15f ", vec[i]);
	}
	printf("\n");
}

int main(int argc, char* argv[])
{
    using mpfr::mpreal;    
    using std::cout;
    using std::endl;
    
    // Required precision of computations in decimal digits
    // Play with it to check different precisions
    const int digits = 50; 

    // Setup default precision for all subsequent computations
    // MPFR accepts precision in bits - so we do the conversion 
    mpreal::set_default_prec(mpfr::digits2bits(digits));

    // Compute all the vital characteristics of mpreal (in current precision)
    // Analogous to lamch from LAPACK
    const mpreal one         =    1.0;
    const mpreal zero        =    0.0;
    const mpreal eps         =    std::numeric_limits<mpreal>::epsilon();
    const int    base        =    std::numeric_limits<mpreal>::radix;
    const mpreal prec        =    eps * base;
    const int bindigits      =    std::numeric_limits<mpreal>::digits(); // eqv. to mpfr::mpreal::get_default_prec();
    const mpreal rnd         =    std::numeric_limits<mpreal>::round_error();
    const mpreal maxval      =    std::numeric_limits<mpreal>::max();
    const mpreal minval      =    std::numeric_limits<mpreal>::min();
    const mpreal small       =    one / maxval;
    const mpreal sfmin       =    (small > minval) ? small * (one + eps) : minval;
    const mpreal round       =    std::numeric_limits<mpreal>::round_style();
    const int    min_exp     =    std::numeric_limits<mpreal>::min_exponent;
    const mpreal underflow   =    std::numeric_limits<mpreal>::min();
    const int    max_exp     =    std::numeric_limits<mpreal>::max_exponent;
    const mpreal overflow    =    std::numeric_limits<mpreal>::max();
    
    // Additionally compute pi with required accuracy - just for fun :)
    const mpreal pi          =    mpfr::const_pi();
        
    cout.precision(digits);    // Show all the digits
    cout << "pi         =    "<<    pi          << endl;    
    cout << "eps        =    "<<    eps         << endl;
    cout << "base       =    "<<    base        << endl;
    cout << "prec       =    "<<    prec        << endl;
    cout << "b.digits   =    "<<    bindigits   << endl;
    cout << "rnd        =    "<<    rnd         << endl;
    cout << "maxval     =    "<<    maxval      << endl;    
    cout << "minval     =    "<<    minval      << endl;    
    cout << "small      =    "<<    small       << endl;    
    cout << "sfmin      =    "<<    sfmin       << endl;    
    cout << "1/sfmin    =    "<<    1 / sfmin   << endl;    
    cout << "round      =    "<<    round       << endl;    
    cout << "max_exp    =    "<<    max_exp     << endl;
    cout << "min_exp    =    "<<    min_exp     << endl;
    cout << "underflow  =    "<<    underflow   << endl;    
    cout << "overflow   =    "<<    overflow    << endl;

    mpreal a, b, c;

    a = 2;
    b = 3;
    c = a + b;

    cout << c << endl;

    return 0;
}
