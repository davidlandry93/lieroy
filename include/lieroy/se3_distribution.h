#ifndef LIEROY_SE3_DISTRIBUTION_H
#define LIEROY_SE3_DISTRIBUTION_H

#include <memory>

#include "lieroy/se3.hpp"

namespace lieroy {
template <typename T>
class SE3Distribution {
  public:
    SE3<T> sample() const;
    virtual std::tuple<SE3<T>, SE3<T>> sample_with_perturbation() const = 0;
    virtual std::unique_ptr<SE3Distribution<T>> copy() const = 0;

  private:
};
}

#endif
