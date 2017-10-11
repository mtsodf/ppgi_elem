#ifndef JACOBIH
#define JACOBIH

#include "../definitions.h"

void Jacobi(int n, real * matrix, real *F, real * initialGuess, real tolerance, int maxiters);

void cg(int n, real *matrix, real *b, real *x);

#endif
