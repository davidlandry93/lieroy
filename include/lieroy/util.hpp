#ifndef LIEROY_UTIL_HPP
#define LIEROY_UTIL_HPP

#include <cmath>
#include <random>

#include <Eigen/Core>

template <typename T>
Eigen::Matrix<T, 3, 1> sample_uniform_from_sphere(T distribution_radius) {
    static std::default_random_engine gen;
    std::normal_distribution<T> normal_distribution;
    std::uniform_real_distribution<T> uniform_dist(0.0, 1.0);

    Eigen::Matrix<T, 3, 1> sample;

    for (auto i = 0; i < 3; i++) {
        sample(i) = normal_distribution(gen);
    }

    double rescaled_radius = cbrt(uniform_dist(gen));
    sample = distribution_radius * rescaled_radius * sample / sample.norm();

    return sample;
}

#endif
