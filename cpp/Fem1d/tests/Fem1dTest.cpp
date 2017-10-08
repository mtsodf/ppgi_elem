#define CATCH_CONFIG_MAIN

#include <stdio.h>
#include "../Fem1d.h"
#include "../../ThirdParty/catch.hpp"

TEST_CASE("Test Local Matrix", "[localmatrix]"){
    real lm[4];

    LocalMatrix(1.0, 1.0, 0.1, lm);

    REQUIRE(20.06666667/2 == Approx(lm[0]).epsilon(0.0001));
    REQUIRE(-9.98333333   == Approx(lm[1]).epsilon(0.0001));
    REQUIRE(-9.98333333   == Approx(lm[2]).epsilon(0.0001));
    REQUIRE(20.06666667/2 == Approx(lm[3]).epsilon(0.0001));

}
