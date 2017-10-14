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


void RhsLocal(real he, real f1, real f2, real *Fe);
void RhsGlobal(int n, real h, real *fs, real *F);
void RhsGlobal(int n, real *hs, real *fs, real *F);

#endif
