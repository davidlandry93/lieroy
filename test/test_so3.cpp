
#include <cmath>

#include <Eigen/Core>
#include <gtest/gtest.h>

#include "lieroy/so3.hpp"
#include "lieroy/so3_uniform_distribution.hpp"

namespace lieroy {

TEST(SO3, LogMapZeroAngleTest) {
  AlgebraSO3<float> lie(0, 0, 0);

  SO3<float> exp_map = lie.exp();
  Eigen::Matrix<float, 3, 3> exp_map_m = exp_map.as_matrix();

  Eigen::Matrix<float, 3, 3> identity = Eigen::Matrix3f::Identity();

  for (auto i = 0; i < 3; i++) {
    for (auto j = 0; j < 3; j++) {
      ASSERT_FLOAT_EQ(identity(i, j), exp_map_m(i, j));
    }
  }
}

TEST(SO3, ExpMapZeroAngleTest) {
  Eigen::Matrix<float, 3, 3> m;
  m.setIdentity();

  SO3<float> identity_transformation(m);

  AlgebraSO3<float> log_map = identity_transformation.log();
  Eigen::Matrix<float, 3, 1> v = log_map.as_vector();

  for (auto i = 0; i < 3; i++) {
    ASSERT_FLOAT_EQ(0.0, v(i));
  }
}

TEST(SO3, AroundZ) {
  Eigen::Matrix<float, 3, 3> m;
  m.setIdentity();

  float rotation_angle = M_PI / 3;
  m(0, 0) = cos(rotation_angle);
  m(1, 1) = cos(rotation_angle);
  m(0, 1) = -sin(rotation_angle);
  m(1, 0) = sin(rotation_angle);

  SO3<float> rotation(m);
  auto log_map = rotation.log();
  auto v = log_map.as_vector();

  SO3<float> rebuilt_rotation = log_map.exp();

  ASSERT_FLOAT_EQ(cos(rotation_angle), rebuilt_rotation.as_matrix()(0, 0));

  ASSERT_FLOAT_EQ(0.0, v(0));
  ASSERT_FLOAT_EQ(0.0, v(1));
  ASSERT_FLOAT_EQ(rotation_angle, v(2));
}

TEST(SO3, CloseToPiOr2Pi) {
  Eigen::Matrix<float, 3, 3> m;
  m << -3.201e-02, 9.994e-01, 1.315e-02, 9.99400000e-01, 3.18300000e-02,
      1.34200000e-02, 1.29900000e-02, 1.357e-02, -9.9982e-01;

  SO3<float> rotation(m);
  ASSERT_NO_THROW(rotation.log());
}

TEST(SO3, Orthogonal) {
    SO3UniformDistribution<double> dist(10.0);

    for(auto k = 0; k < 1000; ++k) {
        SO3<double> sample = dist.sample();

        auto sample_matrix = sample.as_matrix();

        std::cout << "Det: " << sample_matrix.determinant() << '\n';

        std::cout << sample_matrix << '\n';
        sample_matrix *= sample_matrix.transpose();
        auto identity = Eigen::Matrix<float,3,3>::Identity();

        std::cout << sample_matrix;

        for(auto i = 0; i < 3; i++) {
            for(auto j = 0; j < 3; j++) {
                ASSERT_FLOAT_EQ(sample_matrix(i,j), identity(i,j));
            }
        }
    }
}



}
