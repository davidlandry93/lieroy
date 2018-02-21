#include <vector>

#include <gtest/gtest.h>

#include "lieroy/normal_random_variable.hpp"
#include "lieroy/se3.hpp"
#include "lieroy/se3_gaussian_distribution.hpp"

using namespace lieroy;

TEST(DistributionFromSample, PureTranslation) {
    Eigen::Matrix3d covariance;
    covariance.setIdentity();
    NormalRandomVariable<double,3> random_variable(covariance);

    int n_samples = 1000;
    std::vector<SE3<double>> se3_samples(n_samples);
    for(int i = 0; i < n_samples; i++) {
        Eigen::Matrix<double,4,4> m;
        m.setIdentity();
        m.block(0,3,3,1) = random_variable();

        se3_samples[i] = m;
    }

    auto distribution = SE3GaussianDistribution<double>::from_sample(se3_samples);

    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            if(i == j) {
                // Diagonal is all 1.0.
                ASSERT_FLOAT_EQ(1.0, distribution.mean(i,j));
            } else if (j != 3) {
                // Non-translation terms are 0.
                ASSERT_FLOAT_EQ(0.0, distribution.mean(i,j));
            }
        }
    }

    for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) {
            if(i >= 3 || j >= 3) {
                // Non-translation terms are 0.
                ASSERT_FLOAT_EQ(0.0, distribution.covariance(i,j));
            }
        }
    }
}

TEST(DistributionFromSample, AroundZ) {
    Eigen::Matrix<double,6,6> cov;
    cov.setZero();
    cov(5,5) = 1.0;
    NormalRandomVariable<double,6> random_variable(cov);

    int n_samples = 1000;
    std::vector<SE3<double>> se3_samples(n_samples);
    for(int i = 0; i < n_samples; i++) {
        Eigen::Matrix<double,6,1> rpy_vector = random_variable();

        Eigen::Matrix<double,4,4> m;
        m.setIdentity();
        m.block(0,0,3,3) = SO3<double>::from_rpy(rpy_vector(3), rpy_vector(4), rpy_vector(5)).as_matrix();

        se3_samples[i] = SE3<double>(m);
    }

    auto distribution = SE3GaussianDistribution<double>::from_sample(se3_samples);

    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            // Only the upper left part of the matrix is non identity.
            if(i == j && j >= 2) {
                ASSERT_FLOAT_EQ(1.0, distribution.mean(i,j));
            } else if (j >= 2 && i >= 2) {
                ASSERT_FLOAT_EQ(0.0, distribution.mean(i,j));
            }
        }
    }

    for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) {
            if(i != 5 && j != 5) {
                ASSERT_FLOAT_EQ(0.0, distribution.covariance(i,j));
            }
        }
    }
}

TEST(CovarianceConversion, PureTranslation) {
  Eigen::Matrix<double,6,6> covariance;
  covariance.setZero();

  covariance.block(0,0,3,3) = Eigen::Matrix<double,3,3>::Identity();

  auto new_cov = SE3GaussianDistribution<double>::covariance_rpy_to_lie(covariance);

  for(auto i = 3; i < 6; i++) {
    for(auto j = 3; j < 6; j++) {
      ASSERT_FLOAT_EQ(0.0, new_cov(i,j));
    }
  }
}

TEST(CovarianceConversion, PureRotation) {
  Eigen::Matrix<double,6,6> covariance;
  covariance.setZero();

  covariance.block(3,3,3,3) = Eigen::Matrix<double,3,3>::Identity();

  auto new_cov = SE3GaussianDistribution<double>::covariance_rpy_to_lie(covariance);

  for(auto i = 0; i < 3; i++) {
    for(auto j = 0; j < 3; j++) {
      ASSERT_FLOAT_EQ(0.0, new_cov(i,j));
    }
  }
}
