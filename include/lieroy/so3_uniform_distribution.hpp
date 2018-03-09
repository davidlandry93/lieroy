#ifndef LIEROY_SO3_UNIFORM_DISTRIBUTION_HPP
#define LIEROY_SO3_UNIFORM_DISTRIBUTION_HPP

#include "util.hpp"

#include "so3_uniform_distribution.h"

namespace lieroy {
template <typename T>
SO3UniformDistribution<T>::SO3UniformDistribution(const T& radius)
    : mean(SO3<T>::identity()), radius(radius) {}

template <typename T>
SO3UniformDistribution<T>::SO3UniformDistribution(const SO3<T>& mean,
                                                  const T& radius)
    : mean(mean), radius(radius) {}

template <typename T>
std::tuple<SO3<T>, AlgebraSO3<T>> SO3UniformDistribution<T>::sample_with_perturbation() const {
    auto rotation_vector = sample_uniform_from_sphere(radius);
    auto lie_rotation = AlgebraSO3<T>(rotation_vector);

    return std::make_tuple(lie_rotation.exp() * mean, lie_rotation);
}
}

#endif
