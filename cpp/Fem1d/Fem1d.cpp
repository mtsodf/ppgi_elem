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


void RhsLocal(real he, real f1, real f2, real *Fe){

    real w = sqrt(3)/3.0;
    Fe[0]  = f1*(he/2)*(FuncForm(-w,0) * FuncForm(-w,0) + FuncForm(w,0) * FuncForm(w,0));
    Fe[0] += f2*(he/2)*(FuncForm(-w,0) * FuncForm(-w,1) + FuncForm(w,0) * FuncForm(w,1));


    Fe[1]  = f1*(he/2)*(FuncForm(-w,1) * FuncForm(-w,0) + FuncForm(w,1) * FuncForm(w,0));
    Fe[1] += f2*(he/2)*(FuncForm(-w,1) * FuncForm(-w,1) + FuncForm(w,1) * FuncForm(w,1));
}

void RhsGlobal(int n, real h, real *fs, real *F){
    real *hs = (real*) malloc(n*sizeof(real));

    for (size_t i = 0; i < n; i++)
    {
        hs[i] = h;
    }

    RhsGlobal(n, hs, fs, F);

    free(hs);
}

void RhsGlobal(int n, real *hs, real *fs, real *F){
    int ndofs = n - 1;

    for (size_t i = 0; i < ndofs; i++)
    {
        F[i] = 0.0;
    }

    real Fe[2];
    for (size_t i = 0; i < n; i++)
    {
        if(i==0)  {
            RhsLocal(hs[i], 0.0, fs[i],Fe);
        } else if(i==n-1) 
        {
            RhsLocal(hs[i],fs[i-1], 0.0, Fe);
        } else{
            RhsLocal(hs[i],fs[i-1], fs[i], Fe);
        }
        
        if((i-1)>=0) F[i-1] += Fe[0];
        if(i<ndofs) F[i]  += Fe[1];

    }
}