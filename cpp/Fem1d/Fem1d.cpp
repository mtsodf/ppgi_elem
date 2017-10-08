#include "Fem1d.h"
#include "../definitions.h"
#include <stdlib.h>
#include <cmath>
using namespace std;

class UndefinedFuncForm: public exception
{
  virtual const char* what() const throw()
  {
    return "Undefined func form";
  }
};

real FuncForm(real qsi, char func){
    if(func == 0) return (1.0 - qsi)/2.0;
    if(func == 1) return (1.0 + qsi)/2.0;
    UndefinedFuncForm ex;
    throw ex;
}

real DFuncForm(real qsi, char func){
    if(func == 0) return -1.0/2.0;
    if(func == 1) return  1.0/2.0;
    UndefinedFuncForm ex;
    throw ex;
}

void LocalMatrix(real alpha, real beta, real he, real* lm){

    real w;

    w = sqrt(3.0)/3.0;

    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            lm[2*i+j]  = (he/2)*beta*(FuncForm(-w,i)*FuncForm(-w,j)) + (2/he)*alpha*(DFuncForm(-w,i)*DFuncForm(-w,j));
            lm[2*i+j] += (he/2)*beta*(FuncForm( w,i)*FuncForm( w,j)) + (2/he)*alpha*(DFuncForm( w,i)*DFuncForm( w,j));
        }
    }
}

void GlobalMatrix(int n, real alpha, real beta, real *K){
    real * hs;

    real h = 1.0/n;
    hs = (real*) malloc(sizeof(real)*n);

    for (size_t i = 0; i < n; i++)
    {
        hs[i] = h;
    }

    GlobalMatrix(n, alpha, beta, hs, K);
}

void GlobalMatrix(int n, real alpha, real beta, real* hs, real *K){

    real lm[4];
    int ndofs = n - 1;
    for (size_t i = 0; i < n; i++)
    {
        LocalMatrix(alpha, beta, hs[i], lm);

        if(i > 0 && i < ndofs){
            K[DIM(i-1,i-1,ndofs)] += lm[0];
            K[DIM(i-1,i+0,ndofs)] += lm[1];
            K[DIM(i+0,i-1,ndofs)] += lm[2];
            K[DIM(i+0,i+0,ndofs)] += lm[3];
        }

        if(i == 0){
            K[DIM(i,i,ndofs)] += lm[3];
        }

        if(i == ndofs){
            K[DIM(i-1,i-1,ndofs)] += lm[0];
        }
    }


}

