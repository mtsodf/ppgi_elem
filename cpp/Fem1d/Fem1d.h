#ifndef FEM1D_H
#define FEM1D_H

#include <exception>
#include <stdio.h>
#include "../definitions.h"
#define DIRICHLET 0
#define NEUMANN   1
#define DIRICHLET_NEUMANN 2
#define NEUMANN_DIRICHLET 3

real FuncForm(real qsi, char func);
real DFuncForm(real qsi, char func);
void LocalMatrix(real alpha, real beta, real he, real* lm);
void GlobalMatrix(int n, real alpha, real beta, real* h, real *K, int boundary);
void GlobalMatrix(int n, real alpha, real beta, real *K, int boundary);


void RhsLocal(real he, real f1, real f2, real *Fe);
void RhsGlobal(int n, real h, real *fs, real *F, real p, real q, real alpha, real beta, int boundary);
void RhsGlobal(int n, real *hs, real *fs, real *F, real p, real q, real alpha, real beta, int boundary);

#endif
