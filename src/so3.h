#ifndef PYLIE_SO3_H
#define PYLIE_SO3_H

#include <array>

#include <Eigen/Core>

#include "algebra_so3.hpp"

namespace pylie {
    template <class T>
    class SO3 {
        const T SMALL_ANGLE_THRESHOLD = 1e-2;
    public:
        SO3();
        SO3(const Eigen::Matrix<T,3,3>& m);
        SO3(const std::array<T, 9>& m);
        AlgebraSO3<T> log() const;
        AlgebraSO3<T> log_rodrigues() const;
        Eigen::Matrix<T,3,3> as_matrix() const;
        static SO3<T> from_rpy(T roll, T pitch, T yaw);
    private:
        Eigen::Matrix<T,3,3> matrix;
    };
}

#endif
