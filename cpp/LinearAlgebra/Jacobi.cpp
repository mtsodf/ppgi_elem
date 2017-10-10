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
    *residue = Norm(n, v);
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





