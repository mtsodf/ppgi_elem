#include "Fem1d.h"
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
