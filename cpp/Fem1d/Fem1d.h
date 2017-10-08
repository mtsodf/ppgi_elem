#ifndef FEM1D_H
#define FEM1D_H

#include <exception>
#include <stdio.h>
#include "../definitions.h"


real FuncForm(real qsi, char func);
real DFuncForm(real qsi, char func);
void LocalMatrix(real alpha, real beta, real he, real* lm);
void GlobalMatrix(int n, real alpha, real beta, real* h, real *K);
void GlobalMatrix(int n, real alpha, real beta, real *K);

#endif
