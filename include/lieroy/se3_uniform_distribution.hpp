#ifndef LIE_SE3_UNIFORM_DISTRIBUTION_HPP
#define LIE_SE3_UNIFORM_DISTRIBUTION_HPP

#include <cmath>
#include <random>

#include "algebra_se3.hpp"

#include "se3_uniform_distribution.h"

namespace lieroy {
template <typename T>
SE3UniformDistribution<T>::SE3UniformDistribution(const T& translation_radius,
                                                  const T& rotation_radius)
    : radius_translation(translation_radius),
      radius_rotation(rotation_radius) {}

template <typename T>
SE3UniformDistribution<T>::SE3UniformDistribution(const SE3<T>& mean,
                                                  const T& translation_radius,
                                                  const T& rotation_radius)
    : mean(mean),
      radius_translation(translation_radius),
      radius_rotation(rotation_radius) {}

template <typename T>
std::unique_ptr<SE3Distribution<T>> SE3UniformDistribution<T>::copy() const {
    return std::unique_ptr<SE3Distribution<T>>(
        new SE3UniformDistribution(*this));
}

template <typename T>
Eigen::Matrix<T, 3, 1> SE3UniformDistribution<T>::sample_from_sphere(
    T radius) const {
    static std::mt19937 gen{std::random_device{}()};
    static std::normal_distribution<T> dist;
    static std::uniform_real_distribution<T> uniform_dist(0.0, radius);

    Eigen::Matrix<T, 3, 1> p;

    for (auto i = 0; i < 3; i++) {
        p(i) = dist(gen);
    }

    double r = cbrt(uniform_dist(gen));
    p = r * p / p.norm();

    return p;
}

template <typename T>
std::tuple<SE3<T>, SE3<T>> SE3UniformDistribution<T>::sample_with_perturbation()
    const {
    auto translation_sample = sample_from_sphere(radius_translation);
    auto rotation_sample = sample_from_sphere(radius_rotation);

    auto lie_vector = AlgebraSE3<T>(translation_sample, rotation_sample);

    return std::make_tuple(lie_vector.exp() * mean, lie_vector.exp());
}
}

#endif
