#ifndef LIEROY_SE3_UNIFORM_DISTRIBUTION_HPP
#define LIEROY_SE3_UNIFORM_DISTRIBUTION_HPP

#include <cmath>
#include <random>

#include "algebra_se3.hpp"
#include "se3_uniform_distribution.h"
#include "so3_uniform_distribution.hpp"
#include "util.hpp"

namespace lieroy {
template <typename T>
SE3UniformDistribution<T>::SE3UniformDistribution(const T& translation_radius,
                                                  const T& rotation_radius)
    : translation_radius(translation_radius),
      rotation_radius(rotation_radius) {}

template <typename T>
SE3UniformDistribution<T>::SE3UniformDistribution(const SE3<T>& mean,
                                                  const T& translation_radius,
                                                  const T& rotation_radius)
    : mean(mean),
      translation_radius(translation_radius),
      rotation_radius(rotation_radius) {}

template <typename T>
std::unique_ptr<SE3Distribution<T>> SE3UniformDistribution<T>::copy() const {
    return std::unique_ptr<SE3Distribution<T>>(
        new SE3UniformDistribution(*this));
}

template <typename T>
std::tuple<SE3<T>, SE3<T>> SE3UniformDistribution<T>::sample_with_perturbation()
    const {
    auto translation_sample = sample_uniform_from_sphere<T>(translation_radius);
    auto rotation_sample = sample_uniform_from_sphere<T>(rotation_radius);

    auto lie_vector = AlgebraSE3<T>(translation_sample, rotation_sample);

    return std::make_tuple(lie_vector.exp() * mean, lie_vector.exp());
}
}

#endif
