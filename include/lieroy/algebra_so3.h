#ifndef LIEROY_ALGEBRA_SO3_H
#define LIEROY_ALGEBRA_SO3_H

#include <array>

#include <Eigen/Core>

namespace lieroy {
    template <class T>
    class SO3;

    template <class T>
    class AlgebraSO3 {
        const T SMALL_ANGLE_THRESHOLD = 1e-6;
    public:
        AlgebraSO3();
        AlgebraSO3(T w1, T w2, T w3);
        AlgebraSO3(const std::array<T, 3>& v);
        AlgebraSO3(const Eigen::Matrix<T,3,1>& v);
        Eigen::Matrix<T,3,1> as_vector() const;
        Eigen::Matrix<T,3,3> as_skewsym_matrix() const;
        SO3<T> exp() const;
        SO3<T> exp_rodrigues() const;

    private:
        Eigen::Matrix<T,3,1> vector;
    };
}

#endif
