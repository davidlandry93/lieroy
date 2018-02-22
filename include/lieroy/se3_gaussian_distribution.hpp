#ifndef LIEROY_SE3_GAUSSIAN_DISTRIBUTION_HPP
#define LIEROY_SE3_GAUSSIAN_DISTRIBUTION_HPP

#include <limits>

#include "normal_random_variable.hpp"
#include "se3_gaussian_distribution.h"

namespace lieroy {
template <typename T>
SE3GaussianDistribution<T>::SE3GaussianDistribution()
    : mean(SE3<T>::identity()),
      covariance(Eigen::Matrix<T, 6, 6>::Identity()) {}

template <typename T>
SE3GaussianDistribution<T>::SE3GaussianDistribution(
    const SE3<T>& mean, const Eigen::Matrix<T, 6, 6>& covariance)
    : mean(mean), covariance(covariance) {}

template <typename T>
SE3GaussianDistribution<T>::SE3GaussianDistribution(
    const SE3GaussianDistribution<T>& other)
    : mean(other.mean), covariance(other.covariance) {}

template <typename T>
std::unique_ptr<SE3Distribution<T>> SE3GaussianDistribution<T>::copy() const {
    return std::unique_ptr<SE3GaussianDistribution<T>>(
        new SE3GaussianDistribution<T>(*this));
}

template <typename T>
SE3GaussianDistribution<T> SE3GaussianDistribution<T>::from_sample(
    const std::vector<SE3<T>>& sample) {
    return from_sample(sample, std::vector<double>(sample.size(), 1.0));
}

template <typename T>
SE3GaussianDistribution<T> SE3GaussianDistribution<T>::from_sample(
    const std::vector<SE3<T>>& sample, const std::vector<double>& weights) {
    if (sample.empty()) {
        throw std::runtime_error("Empty distribution");
    }

    double sum_of_weights = sum(weights);

    SE3<T> mean_transformation = sample[0];       // The average of the samples.
    Eigen::Matrix<T, 6, 1> average_perturbation;  // The average perturbation.
    Eigen::Matrix<T, 6, 1> previous_average_perturbation;
    Eigen::Matrix<T, 6, 6> covariance;

    // TODO: Replace n_iterations by a real convergence checker.
    int n_iterations = 0;
    auto inf = std::numeric_limits<T>::infinity();
    previous_average_perturbation << inf, inf, inf, inf, inf, inf;
    bool converged = false;

    while (n_iterations < 8 and !converged) {
        average_perturbation.setZero();
        covariance.setZero();

        auto inv_of_mean = mean_transformation.inv();

#pragma omp parallel for schedule(dynamic, 32)
        for (auto i = 0; i < sample.size(); i++) {
            SE3<T> product = sample[i] * inv_of_mean;

            AlgebraSE3<T> delta = product.log();
            auto vector = delta.as_vector();

            auto weighted_vector = weights[i] * vector;
            auto weighted_covariance = weights[i] * vector * vector.transpose();

            average_perturbation += weighted_vector;
            covariance += weighted_covariance;
        }

        average_perturbation /= sum_of_weights;
        covariance /= (std::max(sum_of_weights - 1.0, 1.0));

        auto perturbation = AlgebraSE3<T>(average_perturbation).exp();
        mean_transformation = perturbation * mean_transformation;

        ++n_iterations;
        auto delta =
            (previous_average_perturbation - average_perturbation).norm();
        converged = delta < 1e-5;
        previous_average_perturbation = average_perturbation;
    }

    return SE3GaussianDistribution<T>(mean_transformation, covariance);
}

template <typename T>
T SE3GaussianDistribution<T>::sum(const std::vector<T>& scalars) {
    T sum = 0.0;
    for (auto i = 0; i < scalars.size(); i++) {
        sum += scalars[i];
    }
    return sum;
}

template <typename T>
Eigen::Matrix<T, 6, 6> SE3GaussianDistribution<T>::covariance_rpy_to_lie(
    const Eigen::Matrix<T, 6, 6>& covariance_rpy) {
    const int sample_size = 5000;
    NormalRandomVariable<T, 6> random_distribution(covariance_rpy);

    std::vector<SE3<double>> se3_sample(sample_size);
    for (int i = 0; i < sample_size; i++) {
        auto s = random_distribution();

        Eigen::Matrix<double, 6, 1> rpy_vector(s);

        SO3<double> rotation = SO3<double>::from_rpy(
            rpy_vector(3, 0), rpy_vector(4, 0), rpy_vector(5, 0));
        se3_sample[i] = SE3<double>(rpy_vector.block(0, 0, 3, 1), rotation);
    }

    auto se3_distribution = SE3GaussianDistribution<T>::from_sample(se3_sample);

    return se3_distribution.covariance;
}

template <typename T>
std::tuple<SE3<T>, SE3<T>> SE3GaussianDistribution<T>::sample_with_perturbation() const {
    NormalRandomVariable<T, 6> normal_variable(covariance);
    Eigen::Matrix<T, 6, 1> sampled_vector = normal_variable();
    AlgebraSE3<T> perturbation_algebra(sampled_vector);

    return std::make_tuple(perturbation_algebra.exp() * mean,
                           perturbation_algebra.exp());
}
}

#endif
