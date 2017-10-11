#include <stdlib.h>
#include <stdio.h>
#include "../definitions.h"
#include "Jacobi.h"
#include "Operations.h"



void JacobiStep(int n, real *F, real* matrix, real *v, real *vnext, real *residue){
    for (size_t i = 0; i < n; i++)
    {
        vnext[i] = F[i];
        for (size_t j = 0; j < n; j++)
        {
            if(j!=i){
                vnext[i] -= matrix[DIM(i,j,n)]*v[j];
            }
        }
        vnext[i] /= matrix[DIM(i,i,n)];
    }

    daxpy(n, -1.0, vnext, v);
    *residue = norm(n, v);
}

void Jacobi(int n, real * matrix, real *F, real* initialGuess, real tolerance, int maxiters){
    int iters = 0;

    real* aux = (real*) malloc(sizeof(real)*n);
    real residue;

    do {
        JacobiStep(n, F, matrix, initialGuess, aux, &residue);
        iters++;
        JacobiStep(n, F, matrix, aux, initialGuess, &residue);
        iters++;
    }while(iters < maxiters && residue > tolerance);

}



void cg(int n, real *matrix, real *b, real *x){
    real phi, phiant, beta, alpha;
    real *r, *p, *q;
    daxpy(n, -1.0, matrix, x, 1.0, b);
    
    copy(n, b, &r);
    zero(n, q);
    int iters = 0;

    do{ 
        phiant = phi;
        phi = dot(n, r, r);
        if(iters == 0){
            copy(n, r, &p);
        } else{
            beta = phi/phiant;

            // p = r + beta p
            daxpy(n, 1.0, r, beta, p);
        }
        daxpy(n, 1.0, matrix, p, 0.0, q);

        alpha = phi/dot(n, p, q);

        // x = x + alpha * q
        daxpy(n,  alpha, p, x);

        //r = r - alpha * q
        daxpy(n, -alpha, q, r);
        iters++;
    }while(iters < 1000 && norm(n, r) > 1e-6);

    free(q);
    free(p);
    free(r)
}



