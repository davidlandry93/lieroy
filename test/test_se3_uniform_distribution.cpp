
#include <iostream>

#include <gtest/gtest.h>
#include <Eigen/Core>

#include "lieroy/se3.hpp"
#include "lieroy/se3_uniform_distribution.hpp"

using namespace lieroy;
using namespace std;

TEST(SE3UniformDistribution, ZeroRadiusReturnTheMean) {
    std::array<float, 6> v{{1, 2, 3, 0.1, 0.2, 0.3}};
    AlgebraSE3<float> aMean(v);
    SE3UniformDistribution<float> distribution(aMean.exp(), 0.0, 0.0);

    auto sample = distribution.sample().log();

    std::cout << sample.as_vector() << std::endl;

    for (auto i = 0; i < 6; i++) {
        ASSERT_FLOAT_EQ(aMean.as_vector()[i], sample.as_vector()[i]);
    }
}

TEST(SE3UniformDistribution, NoSampleIsLargerThanRadius) {
    std::array<float, 6> v{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    AlgebraSE3<float> aMean(v);
    const float aRadius = 0.2;
    SE3UniformDistribution<float> distribution(aMean.exp(), aRadius, aRadius);

    for (int i = 0; i < 10000; ++i) {
        auto sample = distribution.sample().log();

        for (auto i = 0; i < 6; i++) {
            ASSERT_LT(fabs(aMean.as_vector()[i]), aRadius);
        }
    }
}
