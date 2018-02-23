#ifndef LIEROY_SE3_UNIFORM_DISTRIBUTION_H
#define LIEROY_SE3_UNIFORM_DISTRIBUTION_H

#include <memory>

#include <Eigen/Core>

#include "lieroy/se3.hpp"
#include "lieroy/se3_distribution.hpp"

namespace lieroy {
template <typename T>
class SE3UniformDistribution : public SE3Distribution<T> {
  public:
    SE3UniformDistribution(const T& translation_radius, const T& rotation_radius);
    SE3UniformDistribution(const SE3<T>& mean, const T& translation_radius, const T& rotation_radius);
    std::unique_ptr<SE3Distribution<T>> copy() const override;
    std::tuple<SE3<T>, SE3<T>> sample_with_perturbation() const;

  private:
    SE3<T> mean;
    T translation_radius = 1.0;
    T rotation_radius = 1.0;

};
}

#endif
