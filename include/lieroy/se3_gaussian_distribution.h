#ifndef LIEROY_SE3_GAUSSIAN_DISTRIBUTION_H
#define LIEROY_SE3_GAUSSIAN_DISTRIBUTION_H

#include <Eigen/Core>

#include "lieroy/se3.hpp"
#include "lieroy/se3_distribution.h"

namespace lieroy {
template <typename T>
class SE3GaussianDistribution : public SE3Distribution<T> {
public:
  SE3<T> mean;
  Eigen::Matrix<T, 6, 6> covariance;

  SE3GaussianDistribution();
  SE3GaussianDistribution(const SE3GaussianDistribution &other);
  SE3GaussianDistribution(const SE3<T> &mean,
                          const Eigen::Matrix<T, 6, 6> &covariance);
  std::unique_ptr<SE3Distribution<T>> copy() const override;
  SE3<T> sample() const override;
  static SE3GaussianDistribution
  from_sample(const std::vector<SE3<T>> &sample);
  static SE3GaussianDistribution
  from_sample(const std::vector<SE3<T>> &sample,
              const std::vector<double> &weights);
  static Eigen::Matrix<T, 6, 6>
  covariance_rpy_to_lie(const Eigen::Matrix<T, 6, 6> &covariance_rpy);

private:
  static T sum(const std::vector<T> &scalars);
};
}

#endif
