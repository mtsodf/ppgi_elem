#define CATCH_CONFIG_MAIN

#include <stdio.h>
#include "../Fem1d.h"
#include "../../Utils/Utils.h"
#include "../../ThirdParty/catch.hpp"

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


