#ifndef LIEROY_SE3_GAUSSIAN_DISTRIBUTION_H
#define LIEROY_SE3_GAUSSIAN_DISTRIBUTION_H

#include <Eigen/Core>

#include "lieroy/algebra_se3.hpp"
#include "lieroy/se3.hpp"
#include "lieroy/se3_distribution.hpp"

namespace lieroy {
template <typename T>
class SE3GaussianDistribution : public SE3Distribution<T> {
  public:
    SE3<T> mean;
    Eigen::Matrix<T, 6, 6> covariance;

    SE3GaussianDistribution();
    SE3GaussianDistribution(const SE3GaussianDistribution& other);
    SE3GaussianDistribution(const SE3<T>& mean,
                            const Eigen::Matrix<T, 6, 6>& covariance);
    std::unique_ptr<SE3Distribution<T>> copy() const override;
    static SE3GaussianDistribution from_sample(
        const std::vector<SE3<T>>& sample);
    static SE3GaussianDistribution from_sample(const SE3Vector& sample);
    static SE3GaussianDistribution from_sample(
        const std::vector<SE3<T>>& sample, const std::vector<double>& weights);
    static Eigen::Matrix<T, 6, 6> covariance_rpy_to_lie(
        const Eigen::Matrix<T, 6, 6>& covariance_rpy);

    std::tuple<SE3<T>, SE3<T>> sample_with_perturbation() const;

  private:
    static T sum(const std::vector<T>& scalars);
};
}

#endif
