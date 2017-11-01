#ifndef PYLIE_SE3_DISTRIBUTION_H
#define PYLIE_SE3_DISTRIBUTION_H

#include <memory>

#include "pylie/se3.hpp"

namespace pylie {
  template<typename T>
  class SE3Distribution {
  public:
    virtual SE3<T> sample() const=0;
    virtual std::unique_ptr<SE3Distribution<T>> copy() const=0;
  private:
  };
}

#endif
