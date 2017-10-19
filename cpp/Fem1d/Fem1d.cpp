#include "Fem1d.h"
#include "../definitions.h"
#include <stdlib.h>
#include "../LinearAlgebra/Jacobi.h"
#include "../LinearAlgebra/Operations.h"
#include <cmath>
#include "../Utils/Utils.h"


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

void GlobalMatrix(int n, real alpha, real beta, real *K, int boundary){
    real * hs;

    real h = 1.0/n;
    hs = (real*) malloc(sizeof(real)*n);

    for (size_t i = 0; i < n; i++)
    {
        hs[i] = h;
    }

    GlobalMatrix(n, alpha, beta, hs, K, boundary);
}

void GlobalMatrix(int n, real alpha, real beta, real* hs, real *K, int boundary){

    real lm[4];
    int ndofs = n + 1;
    for (size_t i = 0; i < n; i++)
    {
        LocalMatrix(alpha, beta, hs[i], lm);
        K[DIM(i,i,ndofs)] += lm[0];
        K[DIM(i,i+1,ndofs)] += lm[1];
        K[DIM(i+1,i,ndofs)] += lm[2];
        K[DIM(i+1,i+1,ndofs)] += lm[3];
        
    }

    if(boundary == DIRICHLET || boundary == DIRICHLET_NEUMANN){
        K[DIM(0,0,ndofs)] = 1.0;
        K[DIM(0,1,ndofs)] = 0.0;
        K[DIM(1,0,ndofs)] = 0.0;
    }

    if(boundary == DIRICHLET || boundary == NEUMANN_DIRICHLET){
        K[DIM(ndofs-1,ndofs-2,ndofs)] = 0.0;
        K[DIM(ndofs-2,ndofs-1,ndofs)] = 0.0;
        K[DIM(ndofs-1,ndofs-1,ndofs)] = 1.0;
    }

}


void RhsLocal(real he, real f1, real f2, real *Fe){

    real w = sqrt(3)/3.0;
    Fe[0]  = f1*(he/2)*(FuncForm(-w,0) * FuncForm(-w,0) + FuncForm(w,0) * FuncForm(w,0));
    Fe[0] += f2*(he/2)*(FuncForm(-w,0) * FuncForm(-w,1) + FuncForm(w,0) * FuncForm(w,1));


    Fe[1]  = f1*(he/2)*(FuncForm(-w,1) * FuncForm(-w,0) + FuncForm(w,1) * FuncForm(w,0));
    Fe[1] += f2*(he/2)*(FuncForm(-w,1) * FuncForm(-w,1) + FuncForm(w,1) * FuncForm(w,1));
}

void RhsGlobal(int n, real h, real *fs, real *F, real p, real q, real alpha, real beta, int boundary){
    real *hs = (real*) malloc(n*sizeof(real));

    for (size_t i = 0; i < n; i++)
    {
        hs[i] = h;
    }

    RhsGlobal(n, hs, fs, F, p, q, alpha, beta, boundary);

    free(hs);
}

void RhsGlobal(int n, real *hs, real *fs, real *F, real p, real q, real alpha, real beta, int boundary){
    int ndofs = n + 1;

    for (size_t i = 0; i < ndofs; i++)
    {
        F[i] = 0.0;
    }

    real Fe[2];
    
    for (size_t i = 0; i < n; i++)
    {

        RhsLocal(hs[i], fs[i], fs[i+1], Fe);

        F[i]    += Fe[0];
        F[i+1]  += Fe[1];

    }

    if(boundary == DIRICHLET || boundary == DIRICHLET_NEUMANN){
        F[0] = p;
        real w = sqrt(3)/3.0; real he = hs[0];
        F[1] -= p*((he/2)*beta*(FuncForm(-w,0)*FuncForm(-w,1)) + (2/he)*alpha*(DFuncForm(-w,0)*DFuncForm(-w,1)));
        F[1] -= p*((he/2)*beta*(FuncForm( w,0)*FuncForm( w,1)) + (2/he)*alpha*(DFuncForm( w,0)*DFuncForm( w,1)));
        
    }

    if(boundary == DIRICHLET || boundary == NEUMANN_DIRICHLET){
        real w = sqrt(3)/3.0; real he = hs[n-1];
        F[ndofs - 1] = q;
        F[ndofs-2] -= q*((he/2)*beta*(FuncForm(-w,0)*FuncForm(-w,1)) + (2/he)*alpha*(DFuncForm(-w,0)*DFuncForm(-w,1)));
        F[ndofs-2] -= q*((he/2)*beta*(FuncForm( w,0)*FuncForm( w,1)) + (2/he)*alpha*(DFuncForm( w,0)*DFuncForm( w,1)));
    }

    if(boundary == NEUMANN || boundary == NEUMANN_DIRICHLET){
        F[0] -= alpha*p;
    }

    if(boundary == NEUMANN || boundary == DIRICHLET_NEUMANN){
        F[ndofs-1] += alpha*q;
    }
}



extern "C"{

    real * Fem1dTest(int n, int entrada){

            real *x, *F, *fs, *sol;
            real *K;
            real alpha, beta, p, q;
            int boundary;

            int unknowns = n + 1;

            real h = 1.0/n;

            zero(unknowns, &x);
            zero(unknowns, &F);
            zero(unknowns, &fs);
            zero(unknowns, &sol);
            zero(unknowns*unknowns, &K);

            x[0] = 0.0;
            
            for (size_t i = 1; i < unknowns; i++)
            {
                x[i] = x[i-1] + h;
            }

            /*
            ######################################################################
            # Selecao da entrada
            ######################################################################
            */

            for (size_t i = 0; i < unknowns; i++)
            {
                switch (entrada)
                {
                    case 0:
                        alpha = 1.0;
                        beta = 1.0;
                        fs[i] = 4*M_PI*M_PI*sin(2*M_PI*x[i]) + sin(2*M_PI*x[i]);
                        sol[i] = sin(2*M_PI*x[i]);
                        boundary = DIRICHLET;
                        p = 0.0;
                        q = 0.0;
                        break;
                    case 1:
                        alpha = 1.0;
                        beta = 0.0;
                        fs[i] = 2*alpha;
                        boundary = DIRICHLET;
                        p = 0.0;
                        q = 0.0;
                        break;
                    case 2:
                        alpha = 1.0;
                        beta = 0.5;
                        boundary = DIRICHLET;
                        fs[i] = 0.5*x[i];
                        p = 0.0;
                        q = 1.0;
                        break;
                    case 3:
                        alpha = 0.0;
                        beta = 1.0;
                        boundary = DIRICHLET;
                        fs[i] = beta*x[i]*(1-x[i]);
                        p = 0.0;
                        q = 0.0;
                        break;
                    case 4:
                        alpha = 2.0;
                        beta = 1.0;
                        boundary = DIRICHLET;
                        fs[i] = -7*exp(2*x[i]);
                        p = 1.0;
                        q = exp(2.0);
                        break;
                    case 5:
                        alpha = 2.0;
                        beta = 1.0;
                        boundary = NEUMANN;
                        fs[i] = -7*exp(2*x[i]);
                        p = 2.0;
                        q = 2*exp(2.0);
                        break;
                    case 6:
                        alpha = 2.0;
                        beta = 1.0;
                        boundary = NEUMANN_DIRICHLET;
                        fs[i] = -7*exp(2*x[i]);
                        p = 2.0;
                        q = exp(2.0);
                        break;
                    case 7:
                        alpha = 2.0;
                        beta = 1.0;
                        boundary = DIRICHLET_NEUMANN;
                        fs[i] = -7*exp(2*x[i]);
                        p = 1.0;
                        q = 2*exp(2.0);
                        break;
                    case 8:
                        alpha = 1.0;
                        beta = 1.0;
                        boundary = NEUMANN;
                        fs[i] = 4*M_PI*M_PI*sin(2*M_PI*x[i]) + sin(2*M_PI*x[i]);
                        p =  2*M_PI;
                        q =  2*M_PI;
                        break;
                    default:
                        break;
                }
            }
            
            /*
            ######################################################################
            # Criando Matriz Global
            ######################################################################
            */

            GlobalMatrix(n, alpha, beta, K, boundary);

            /*
            ######################################################################
            # Criando Lado Direito
            ######################################################################
            */
            RhsGlobal(n, h, fs, F, p, q, alpha, beta, boundary);


            /*
            ######################################################################
            # Solucao do Sistema Linear
            ######################################################################
            */            
            real * calc;
            zero(unknowns, &calc);
            cg(unknowns, K, F, calc);

            free(x);
            free(F);
            free(fs);
            free(sol);
            free(K);

            return calc;
    }

}
