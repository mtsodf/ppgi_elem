#include <stdio.h>
#include <stdlib.h>

void saxpy(int n, float a, float *x, float *y){
    for (int i = 0; i < n; ++i) y[i] = a*x[i] + y[i];
}

int main(int argc, char* argv[]){

    int n;

    if(argc < 2){
        printf("Entre com a quantidade de pontos.\n");
        return -1;
    }

    n = atoi(argv[1]);


    printf("Problema de tamanho %d\n", n);



    return 0;
}
