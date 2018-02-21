#ifndef LIE_SE3_UNIFORM_DISTRIBUTION_H
#define LIE_SE3_UNIFORM_DISTRIBUTION_H

#include <memory>

#include <Eigen/Core>

#include "lieroy/se3.hpp"
#include "lieroy/se3_distribution.h"

namespace lie {
  template <typename T>
  class SE3UniformDistribution : public SE3Distribution<T> {
  public:
    SE3UniformDistribution(const T& translation_radius, const T& rotation_radius);
    std::unique_ptr<SE3Distribution<T>> copy() const override;
    SE3<T> sample() const override;
  private:
    T radius_translation=1.0;
    T radius_rotation=1.0;

    Eigen::Matrix<T,3,1> sample_from_sphere(T radius) const;
  };
}

#endif
