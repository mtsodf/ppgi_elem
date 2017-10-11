#define CATCH_CONFIG_MAIN

#include <stdio.h>
#include <stdlib.h>
#include "../../Utils/Utils.h"
#include "../../definitions.h"
#include "../Jacobi.h"
#include "../Operations.h"
#include "../../ThirdParty/catch.hpp"

TEST_CASE("Jacobi test", "[jacobi]"){


    //Laplacian 1D
    real matrix[16];

    for (size_t i = 0; i < 16; i++)
    {
        matrix[i] = 0.0;
    }

    for (size_t i = 0; i < 4; i++) matrix[DIM(i,i,4)] =  2.0;
    for (size_t i = 0; i < 3; i++) matrix[DIM(i,i+1,4)] = -1.0;
    for (size_t i = 0; i < 3; i++) matrix[DIM(i+1,i,4)] = -1.0;


    real F[4];
    F[0]=  0.0;
    F[1]=  0.0;
    F[2]=  0.0;
    F[3]=  5.0;

    real x[4];
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;


    Jacobi(4, matrix, F, x, 1.0e-6, 2000);

    REQUIRE(1.0 == Approx(x[0]).epsilon(1e-5));
    REQUIRE(2.0 == Approx(x[1]).epsilon(1e-5));
    REQUIRE(3.0 == Approx(x[2]).epsilon(1e-5));
    REQUIRE(4.0 == Approx(x[3]).epsilon(1e-5));
}


TEST_CASE("Matmul test", "[matmul]"){
    
        //Laplacian 1D
        real matrix[16];
    
        for (size_t i = 0; i < 16; i++)
        {
            matrix[i] = 0.0;
        }
    
        for (size_t i = 0; i < 4; i++) matrix[DIM(i,i,4)] =  2.0;
        for (size_t i = 0; i < 3; i++) matrix[DIM(i,i+1,4)] = -1.0;
        for (size_t i = 0; i < 3; i++) matrix[DIM(i+1,i,4)] = -1.0;
    
        real F[4];
        F[0]=  0.0;
        F[1]=  0.0;
        F[2]=  0.0;
        F[3]=  5.0;
    
        real x[4];
        x[0] = 1.0;
        x[1] = 2.0;
        x[2] = 3.0;
        x[3] = 4.0;
    
        real result[4];

        result[0] = 0.0;
        result[1] = 0.0;
        result[2] = 0.0;
        result[3] = 0.0;
        
        matmul(4, matrix, x, result);

        REQUIRE(0.0 == Approx(result[0]).epsilon(1e-5));
        REQUIRE(0.0 == Approx(result[1]).epsilon(1e-5));
        REQUIRE(0.0 == Approx(result[2]).epsilon(1e-5));
        REQUIRE(5.0 == Approx(result[3]).epsilon(1e-5));


    }