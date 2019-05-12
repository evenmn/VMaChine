#include <gtest/gtest.h>
#include "../src/Basis/basis.h"

TEST(factorial, output test) {
    ASSERT_EQ(1,                  factorial(0));
    ASSERT_EQ(1,                  factorial(1));
    ASSERT_EQ(2,                  factorial(2));
    ASSERT_EQ(120,                factorial(5));
    ASSERT_EQ(243290200817664000, factorial(20));
}

TEST(factorialDifference, output test) {
    ASSERT_EQ(1,                  factorialDifference(1,0));
    ASSERT_EQ(2,                  factorialDifference(2,0));
    ASSERT_EQ(90,                 factorialDifference(10,8));
    ASSERT_EQ(210,                factorialDifference(15,13));
}

Test(binomial, output test) {
    ASSERT_EQ(2,                  binomial(1,1));
    ASSERT_EQ(3,                  binomial(1,2));
    ASSERT_EQ(6,                  binomial(2,2));
    ASSERT_EQ(4,                  binomial(1,3));
    ASSERT_EQ(10,                 binomial(2,3));
    ASSERT_EQ(136,                binomial(15,2));
    ASSERT_EQ(816,                binomial(15,3));
}



int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
