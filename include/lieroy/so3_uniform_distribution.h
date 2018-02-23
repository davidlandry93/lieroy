#ifndef LIEROY_SO3_UNIFORM_DISTRIBUTION_H
#define LIEROY_SO3_UNIFORM_DISTRIBUTION_H

#include "lieroy/so3.hpp"
#include "lieroy/so3_distribution.hpp"

namespace lieroy {

template <typename T>
class SO3UniformDistribution : public SO3Distribution<T> {
  public:
    SO3UniformDistribution(const T& radius);
    SO3UniformDistribution(const SO3<T>& mean, const T& radius);
    std::tuple<SO3<T>, AlgebraSO3<T>> sample_with_perturbation() const override;

  private:
    SO3<T> mean;
    T radius;
};
}

#endif
