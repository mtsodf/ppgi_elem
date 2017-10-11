#ifndef OPERATIONSH
#define OPERATIONSH

real norm(int n, real *x);
real dot(int n, real *x, real *y);
void daxpy(int n, float alpha, real* x, real* y);
void matmul(int n, real * matrix, real* vector, real* result);
void daxpy(int n, float alpha, real* matrix, real* x, real beta, real* y);
void copy(int n, real *a, real **b);
void zero(int n, real**x);

#endif
