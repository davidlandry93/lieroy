#ifndef LIEROY_SO3_DISTRIBUTION_HPP
#define LIEROY_SO3_DISTRIBUTION_HPP

#include "lieroy/so3_distribution.h"

namespace lieroy {
template <typename T>
SO3<T> SO3Distribution<T>::sample() const {
    SO3<T> sample;
    AlgebraSO3<T> perturbation;
    std::tie(sample, perturbation) = sample_with_perturbation();
    return sample;
}
}
#endif
