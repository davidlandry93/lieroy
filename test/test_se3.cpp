
#include <cmath>

#include <Eigen/Core>
#include <gtest/gtest.h>

#include "lieroy/normal_random_variable.hpp"
#include "lieroy/se3.hpp"
#include "lieroy/se3_gaussian_distribution.hpp"

using namespace lieroy;

TEST(SE3, LogMapZeroAngleTest) {
    Eigen::Matrix<float,4,4> m;
    m.setIdentity();

    SE3<float> identity_transformation(m);

    AlgebraSE3<float> log_map = identity_transformation.log();
    Eigen::Matrix<float,6,1> v = log_map.as_vector();

    for(auto i = 0; i < 6; i++) {
        ASSERT_FLOAT_EQ(0.0, v(i));
    }
}

TEST(SE3, TranslationTest) {
    Eigen::Matrix<float,4,4> m;
    m << 1.0, 0.0, 0.0, 1.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;

    SE3<float> translation(m);
    auto log_map = translation.log();
    auto v = log_map.as_vector();

    ASSERT_FLOAT_EQ(1.0, v(0));

    for(auto i=1; i < 6; i++) {
        ASSERT_FLOAT_EQ(0.0, v(i));
    }
}

TEST(SE3, ProductTest) {
    Eigen::Matrix<float,4,4> m;
    float theta = 0.0;

    m << cos(theta), -sin(theta), 0.0, 1.0,
        sin(theta), cos(theta), 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;

    SE3<float> t1(m);
    SE3<float> t2(m);
    SE3<float> product = t1 * t2;

    ASSERT_FLOAT_EQ(2.0, product(0,3));
}

TEST(SE3, InverseTest) {
    Eigen::Matrix<float,4,4> m;
    float theta = M_PI/2;
    m << cos(theta), -sin(theta), 0.0, 1.0,
        sin(theta), cos(theta), 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;

    SE3<float> t(m);
    SE3<float> inverse = t.inv();
    SE3<float> product = t * inverse;

    auto m_product = product.as_matrix();

    for(int i=0; i < 4; i++) {
        for(int j=0; j < 4; j++) {
            if(i == j) ASSERT_FLOAT_EQ(1.0, m_product(i,j));
            else ASSERT_FLOAT_EQ(0.0, m_product(i,j));
        }
    }
}

TEST(SE3, PureRotationLnTest) {
    SO3<double> rot = SO3<double>::from_rpy(0,0,1.0);
    SE3<double> t(Eigen::Matrix<double,3,1>::Zero(), rot);

    auto log_map = t.log();
    auto vector = log_map.as_vector();

    for(auto i=0; i < 5; i++) {
        ASSERT_FLOAT_EQ(0.0, vector(i));
    }

    ASSERT_FLOAT_EQ(1.0, vector(5));

    std::cout << vector << std::endl;
}

TEST(SE3, PerturbationTest) {
    Eigen::Matrix<float,4,4> m = Eigen::Matrix<float,4,4>::Identity();

    SE3<float> t(m);
    SE3<float> perturbated = t.perturbate(1.0);

    std::cout << perturbated << std::endl;

    ASSERT_NEAR(1.0, perturbated.as_matrix().determinant(), 1e-4);
}



