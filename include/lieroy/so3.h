#ifndef LIEROY_SO3_H
#define LIEROY_SO3_H

#include <array>

#include <Eigen/Core>

#include "lieroy/algebra_so3.hpp"

namespace lieroy {
    template <class T>
    class SO3 {
        const T SMALL_ANGLE_THRESHOLD = 1e-2;
    public:
        SO3();
        SO3(const Eigen::Matrix<T,3,3>& m);
        SO3(const std::array<T, 9>& m);
        SO3(const SO3<T>& other);
        static SO3<T> identity();
        AlgebraSO3<T> log() const;
        AlgebraSO3<T> log_rodrigues() const;
        Eigen::Matrix<T,3,3> as_matrix() const;
        static SO3<T> from_rpy(T roll, T pitch, T yaw);

        Eigen::Matrix<T,3,3> adjoint() const;

        SO3<T>& operator=(const SO3<T>& rhs);
        SO3<T>& operator=(SO3<T>&& rhs);
        SO3<T>& operator*(const SO3<T>& rhs);
    private:
        Eigen::Matrix<T,3,3> matrix;
    };
}

#endif
