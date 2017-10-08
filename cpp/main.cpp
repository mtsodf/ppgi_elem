#include "Fem1d/Fem1d.h"
#include <stdlib.h>


int main(){

    printf("FuncForma %f %f\n", FuncForm(-1,1), FuncForm(1,1));

    real lm[4];

    LocalMatrix(1.0, 1.0, 0.2, lm);

    printf("%f %f %f %f\n", lm[0], lm[1], lm[2], lm[3]);
}
