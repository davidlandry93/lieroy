#ifndef LIEROY_SE3_DISTRIBUTION_HPP
#define LIEROY_SE3_DISTRIBUTION_HPP

#include <tuple>

#include "se3_distribution.h"

namespace lieroy {

template <typename T>
SE3<T> SE3Distribution<T>::sample() const {
    SE3<T> sample, perturbation;
    std::tie(sample, perturbation) = sample_with_perturbation();
    return sample;
}
}

#endif
