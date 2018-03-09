
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
    AlgebraSE3<float> a_mean(v);
    const float a_radius = 0.2;
    SE3UniformDistribution<float> distribution(a_mean.exp(), a_radius, a_radius);

    for (int i = 0; i < 10000; ++i) {
        auto sample = distribution.sample().log();

        for (auto i = 0; i < 6; i++) {
            ASSERT_LT(fabs(a_mean.as_vector()[i]), a_radius);
        }
    }
}

TEST(SE3UniformDistribution, PureTranslationSampling) {
    SE3UniformDistribution<float> distribution(SE3<float>::identity(), 1.0, 0.0);

    for (int i = 0; i < 10; i++) {
        SE3<float> sample = distribution.sample();

        std::cout << sample << std::endl;

        auto id = SE3<float>::identity();

        for(auto j = 0; j < 3; j++) {
            for(auto k = 0; k < 3; k++) {
                ASSERT_FLOAT_EQ(id.as_matrix()(j,k), sample.as_matrix()(j,k));
            }
        }
    }
}


TEST(SE3UniformDistribution, PureRotationSampling) {
    SE3UniformDistribution<float> distribution(SE3<float>::identity(), 0.0, 1.0);

    for (int i = 0; i < 10; i++) {
        auto sample = distribution.sample();

        for(auto j = 0; j < 3; j++) {
            ASSERT_FLOAT_EQ(0.0, sample.as_matrix()(j,3));
        }
    }
}
