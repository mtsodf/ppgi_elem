#include "Utils.h"
#include "../definitions.h"
#include <stdio.h>

void PrintMatrix(int n, real * matrix){
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++){
            printf("%8.2f", matrix[DIM(i,j,n)]);
        }
        printf("\n");
    }
}

void PrintVec(int n, real* vec){
    for (size_t i = 0; i < n; i++)
    {
        printf("%8.4f\n", vec[i]);
    }
}