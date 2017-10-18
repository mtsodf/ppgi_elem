#define CATCH_CONFIG_MAIN

#include <stdio.h>
#include "../Fem1d.h"
#include "../../Utils/Utils.h"
#include "../../ThirdParty/catch.hpp"
#include "../../definitions.h"
#include "../../LinearAlgebra/Jacobi.h"
#include "../../LinearAlgebra/Operations.h"
#include <cmath>

TEST_CASE("Test Local Matrix", "[localmatrix]"){
    real lm[4];

    LocalMatrix(1.0, 1.0, 0.1, lm);

    REQUIRE(20.06666667/2 == Approx(lm[0]).epsilon(0.0001));
    REQUIRE(-9.98333333   == Approx(lm[1]).epsilon(0.0001));
    REQUIRE(-9.98333333   == Approx(lm[2]).epsilon(0.0001));
    REQUIRE(20.06666667/2 == Approx(lm[3]).epsilon(0.0001));

}


TEST_CASE("Test Global Matrix", "[globalmatrix]"){
    real K[9];

    for (size_t i = 0; i < 9; i++)
    {
        K[i] = 0.0;
    }

    GlobalMatrix(4, 1, 1, K);

    PrintMatrix(3, K);

    REQUIRE( 8.16666667 == Approx(K[0]).epsilon(0.00001));
    REQUIRE(-3.95833333 == Approx(K[1]).epsilon(0.00001));
    REQUIRE( 0.00000000 == Approx(K[2]).epsilon(0.00001));
    REQUIRE(-3.95833333 == Approx(K[3]).epsilon(0.00001));
    REQUIRE( 8.16666667 == Approx(K[4]).epsilon(0.00001));
    REQUIRE(-3.95833333 == Approx(K[5]).epsilon(0.00001));
    REQUIRE( 0.00000000 == Approx(K[6]).epsilon(0.00001));
    REQUIRE(-3.95833333 == Approx(K[7]).epsilon(0.00001));
    REQUIRE( 8.16666667 == Approx(K[8]).epsilon(0.00001));
}


TEST_CASE("Test RHS", "[rhsvector]"){
    real Fe[2], a[2], F[3], x[3];
    printf("Ang = %f\t sin = %f\n", 2.0*M_PI*0.25, sin(2.0*M_PI*0.25));

    real hs[4], fs[3];
    hs[0] = 0.25; hs[1] = 0.25; hs[2] = 0.25; hs[3] = 0.25;
    x[0] = 0.25; x[1] = 0.5; x[2] = 0.75;
    for (size_t i = 0; i < 3; i++)
    {
        fs[i] = 4*M_PI*M_PI*sin(2*M_PI*x[i]) + sin(2*M_PI*x[i]);
    }


    RhsLocal(0.25, 0.0, fs[0], Fe);
    RhsLocal(0.25, fs[0], fs[1], a);
    REQUIRE( 6.74640293 == Approx(Fe[1]+a[0]).epsilon(0.00001));

    RhsGlobal(4, hs, fs, F);

    REQUIRE( 6.74640293    == Approx(F[0]).epsilon(0.00001));
    REQUIRE(8.88178420e-16 == Approx(F[1]).epsilon(0.00001));
    REQUIRE(-6.74640293    == Approx(F[2]).epsilon(0.00001));   

}

TEST_CASE("Test Elem Matrix", "[globalmatrix]"){

    int n = 100;
    int unknows = n - 1;

    real *x, *F, *fs, *sol;
    real *K;

    zero(unknows, &x);
    zero(unknows, &F);
    zero(unknows, &fs);
    zero(unknows, &sol);
    zero(unknows*unknows, &K);


    for (size_t i = 0; i < unknows; i++)
    {
        for (size_t j = 0; j < unknows; j++)
        {
            K[DIM(i, j, unknows)] = 0.0;
        }
    }


    real h = 1.0/n;

    x[0] = h;

    for (size_t i = 1; i < unknows; i++)
    {
        x[i] = x[i-1] + h;
    }

    REQUIRE(x[0] == Approx(h));
    REQUIRE(x[unknows-1] == Approx(1-h));

    for (size_t i = 0; i < unknows; i++)
    {
        fs[i] = 4*M_PI*M_PI*sin(2*M_PI*x[i]) + sin(2*M_PI*x[i]);
        sol[i] = sin(2*M_PI*x[i]);
    }



    GlobalMatrix(n, 1, 1, K);
    RhsGlobal(n, h, fs, F);

    printf("Matriz de Rigidez:\n");
    PrintMatrix(unknows, K);

    printf("Calc Lado Direito:\n");
    PrintVec(n-1, F);

    real * calc;

    zero(unknows, &calc);

    printf("Calc Chute Inicial:\n");
    PrintVec(n-1, calc);

    cg(unknows, K, F, calc);

    //Jacobi(unknows, K, F, calc, 1e-6, 1000);

    printf("Calc Solucao:\n");
    PrintVec(n-1, calc);

    real normSol = norm(unknows, sol);
    
    daxpy(n-1, -1.0, calc, sol);

    printf("Norma |calc - sol|/|sol| = %f\n", norm(n-1, sol)/normSol);

    free(x);
    free(F);
    free(fs);
    free(sol);
    free(K);

}   

