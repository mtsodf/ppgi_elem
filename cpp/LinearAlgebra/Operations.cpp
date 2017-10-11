#include <stdio.h>
#include "../definitions.h"
#include <cmath>
#include "Operations.h"
#include <stdlib.h>


real norm(int n, real *x){
    return sqrt(dot(n,x,x));
}

real dot(int n, real *x, real *y){
    real v = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        v += x[i]*y[i];
    }    
    return v;
}

void daxpy(int n, float alpha, real* x, real* y){
    for (size_t i = 0; i < n; i++)
    {
        y[i] = alpha*x[i] + y[i];
    }
}

void matmul(int n, real *matrix, real *vector, real *result){

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            result[i] += matrix[DIM(i,j,n)]*vector[j];
        }
    }
}

void daxpy(int n, float alpha, real* matrix, real* x, real beta, real*y){
    for (size_t i = 0; i < n; i++)
    {
        real aux=0.0;
        for (size_t j = 0; j < n; j++)
        {
            aux += matrix[DIM(i,j,n)]*x[j];
        }
        y[i] = beta*y[i] + alpha*aux;
    }    
}

void copy(int n, real *a, real **b){
    *b = (real*) malloc(n*sizeof(real));
    
    for (size_t i = 0; i < n; i++)
    {
        (*b)[i] = a[i];
    }
}

void zero(int n, real**x){
    *x = (real*) malloc(n*sizeof(real));
    
    for (size_t i = 0; i < n; i++)
    {
        (*x)[i] = 0.0;
    }    
}


