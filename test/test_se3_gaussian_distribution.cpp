
#include <iostream>

#include <Eigen/Core>
#include <gtest/gtest.h>

#include "lieroy/se3.hpp"
#include "lieroy/se3_gaussian_distribution.hpp"

using namespace lieroy;
using namespace std;

TEST(SE3GaussianDistribution, SamplingAlongTranslation) {
    Eigen::Matrix<float,6,6> covariance;
    covariance.setZero();
    covariance(0,0) = 1.0;
    covariance(1,1) = 1.0;
    covariance(2,2) = 1.0;

    SE3GaussianDistribution<float> distribution(SE3<float>::identity(), covariance);
    SE3<float> sample = distribution.sample();

    std::cout << sample.as_matrix() << std::endl;

    for(auto i=0; i<4; i++) {
        ASSERT_FLOAT_EQ(1.0, sample.as_matrix()(i,i));
    }

    for(auto i=3; i < 6; i++) {
      ASSERT_FLOAT_EQ(0.0, distribution.covariance(i,i));
    }
}

TEST(SE3GaussianDistribution, SamplingAlongRotation) {
    Eigen::Matrix<float,6,6> covariance;
    covariance.setZero();
    covariance(3,3) = 1.0;
    covariance(4,4) = 1.0;
    covariance(5,5) = 1.0;

    SE3GaussianDistribution<float> distribution(SE3<float>::identity(), covariance);
    SE3<float> sample = distribution.sample();

    std::cout << sample.as_matrix() << std::endl;

    for(auto i=0; i<2; i++) {
        ASSERT_FLOAT_EQ(0.0, sample.as_matrix()(i,3));
    }
}

