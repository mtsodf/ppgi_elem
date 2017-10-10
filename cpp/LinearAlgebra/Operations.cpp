#include <stdio.h>
#include "../definitions.h"
#include <cmath>
#include "Operations.h"


real Norm(int n, real *x){
    real v = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        v += x[i]*x[i];
    }

    return sqrt(v);
}


void daxpy(int n, float alpha, real* x, real* y){
    for (size_t i = 0; i < n; i++)
    {
        y[i] = alpha*x[i] + y[i];
    }
}
