#ifndef _PCA_H_
#define _PCA_H_


float *vector(int n);
float **matrix(int n, int m);

void free_vector(float *v, int n);
void free_matrix(float **mat, int n, int m);

void print_matrix(float **mat, int n, int m);
void print_vector(float *vec, int n);

void tred2(float **a, int n, float *d, float *e);
void tqli(float *d, float *e, int n, float **z);

#endif /* _PCA_H_ */
