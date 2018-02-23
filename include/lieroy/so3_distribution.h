#ifndef LIEROY_SO3_DISTRIBUTION_H
#define LIEROY_SO3_DISTRIBUTION_H

namespace lieroy {

template <typename T>
class SO3Distribution {
  public:
    SO3<T> sample() const;
    virtual std::tuple<SO3<T>, AlgebraSO3<T>> sample_with_perturbation() const = 0;

  private:
};

}

#endif
