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
    int n = 4;
    int unknowns = n + 1;
    real K[unknowns*unknowns];

    for (size_t i = 0; i < unknowns*unknowns; i++)
    {
        K[i] = 0.0;
    }

    GlobalMatrix(4, 1, 1, K, DIRICHLET);

    REQUIRE( 1.00000000 == Approx(K[DIM(0,0,unknowns)]).epsilon(0.00001));
    REQUIRE( 0.00000000 == Approx(K[DIM(0,1,unknowns)]).epsilon(0.00001));

    REQUIRE( 0.00000000 == Approx(K[DIM(1,0,unknowns)]).epsilon(0.00001));
    REQUIRE( 8.16666667 == Approx(K[DIM(1,1,unknowns)]).epsilon(0.00001));
    REQUIRE(-3.95833333 == Approx(K[DIM(1,2,unknowns)]).epsilon(0.00001));

    REQUIRE(-3.95833333 == Approx(K[DIM(2,1,unknowns)]).epsilon(0.00001));
    REQUIRE( 8.16666667 == Approx(K[DIM(2,2,unknowns)]).epsilon(0.00001));
    REQUIRE(-3.95833333 == Approx(K[DIM(2,3,unknowns)]).epsilon(0.00001));

    REQUIRE(-3.95833333 == Approx(K[DIM(3,2,unknowns)]).epsilon(0.00001));
    REQUIRE( 8.16666667 == Approx(K[DIM(3,3,unknowns)]).epsilon(0.00001));
    REQUIRE( 0.00000000 == Approx(K[DIM(3,4,unknowns)]).epsilon(0.00001));

    REQUIRE( 0.00000000 == Approx(K[DIM(4,3,unknowns)]).epsilon(0.00001));
    REQUIRE( 1.00000000 == Approx(K[DIM(4,4,unknowns)]).epsilon(0.00001));
}


TEST_CASE("Test RHS", "[rhsvector]"){
    real Fe[2], a[2], F[5], x[5];

    real hs[4], fs[5];
    hs[0] = 0.25; hs[1] = 0.25; hs[2] = 0.25; hs[3] = 0.25;
    x[0]  = 0.0 ; x[1] = 0.25; x[2] = 0.5; x[3] = 0.75; x[4] = 1.0;
    for (size_t i = 0; i < 5; i++)
    {
        fs[i] = 4*M_PI*M_PI*sin(2*M_PI*x[i]) + sin(2*M_PI*x[i]);
    }

    RhsLocal(0.25, fs[0], fs[1], Fe);
    RhsLocal(0.25, fs[1], fs[2], a);

    REQUIRE( 6.74640293 == Approx(Fe[1]+a[0]).epsilon(0.00001));


    RhsGlobal(4, hs, fs, F, 0.0, 0.0, 1.0, 1.0, DIRICHLET);
    

    REQUIRE( 0.00000000    == Approx(F[0]).epsilon(0.00001));
    REQUIRE( 6.74640293    == Approx(F[1]).epsilon(0.00001));
    REQUIRE(8.88178420e-16 == Approx(F[2]).epsilon(0.00001));
    REQUIRE(-6.74640293    == Approx(F[3]).epsilon(0.00001));  
    REQUIRE( 0.00000000    == Approx(F[4]).epsilon(0.00001)); 

}

TEST_CASE("Test Elem Matrix", "[globalmatrix]"){

    int n = 100;
    int unknowns = n + 1;

    real *x, *F, *fs, *sol;
    real *K;

    zero(unknowns, &x);
    zero(unknowns, &F);
    zero(unknowns, &fs);
    zero(unknowns, &sol);
    zero(unknowns*unknowns, &K);


    for (size_t i = 0; i < unknowns; i++)
    {
        for (size_t j = 0; j < unknowns; j++)
        {
            K[DIM(i, j, unknowns)] = 0.0;
        }
    }


    real h = 1.0/n;

    x[0] = 0.0;

    for (size_t i = 1; i < unknowns; i++)
    {
        x[i] = x[i-1] + h;
    }

    REQUIRE(x[0] == 0.0);
    REQUIRE(x[1] == Approx(h));
    REQUIRE(x[unknowns-2] == Approx(1-h));
    REQUIRE(x[unknowns-1] == Approx(1));

    for (size_t i = 0; i < unknowns; i++)
    {
        fs[i] = 4*M_PI*M_PI*sin(2*M_PI*x[i]) + sin(2*M_PI*x[i]);
        sol[i] = sin(2*M_PI*x[i]);
    }



    GlobalMatrix(n, 1, 1, K, DIRICHLET);
    RhsGlobal(n, h, fs, F, 0.0, 0.0, 1.0, 1.0, DIRICHLET);


    real * calc;

    zero(unknowns, &calc);


    cg(unknowns, K, F, calc);

    real normSol = norm(unknowns, sol);
    
    daxpy(n-1, -1.0, calc, sol);

    printf("Norma |calc - sol|/|sol| = %f\n", norm(unknowns, sol)/normSol);

    free(x);
    free(F);
    free(fs);
    free(sol);
    free(K);

}   

