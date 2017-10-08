#include "Utils.h"
#include "../definitions.h"
#include <stdio.h>

void PrintMatrix(int n, real * matrix){
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++){
            printf("%6.2f", matrix[DIM(i,j,n)]);
        }
        printf("\n");
    }
}