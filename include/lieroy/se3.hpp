#ifndef LIEROY_SE3_HPP
#define LIEROY_SE3_HPP

#include <algorithm>
#include <iomanip>
#include <random>
#include <stdexcept>

#include "lieroy/algebra_so3.hpp"
#include "lieroy/normal_random_variable.hpp"
#include "lieroy/se3.h"
#include "lieroy/util.hpp"

namespace lieroy {
template <class T>
SE3<T>::SE3() : values {
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0} {}

template <class T>
SE3<T>::SE3(const std::array<T, 16>& p_values) : values{} {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            values[i][j] = p_values[i * 4 + j];
        }
    }
}

template <class T>
SE3<T>::SE3(const Eigen::Matrix<T, 4, 4>& m) : values{} {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            values[i][j] = m(i,j);
        }
    }
}

template <class T>
SE3<T>::SE3(const SE3<T>& other) : values{} {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            values[i][j] = other.values[i][j];
        }
    }
}

template <class T>
SE3<T>::SE3(const Eigen::Matrix<T, 3, 1>& translation, const SO3<T>& p_rotation)
        : values{
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0} {

    for(auto i = 0; i < 3; ++i) {
        for(auto j = 0; j < 3; ++j) {
            values[i][j] = p_rotation.as_matrix()(i,j);
        }

        values[i][3] = translation(i);
    }
}

template <class T>
T& SE3<T>::operator()(int i, int j) {
    return values[i][j];
}

template <class T>
const T& SE3<T>::operator()(int i, int j) const {
    return values[i][j];
}

template <class T>
void SE3<T>::stream_to(std::ostream& os) const {
    os << as_matrix();
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const SE3<T>& t) {
    t.stream_to(os);
    return os;
}

template <class T>
SO3<T> SE3<T>::rotation_part() const {
    auto m = as_matrix();
    return SO3<T>(m.block(0, 0, 3, 3));
}

template <class T>
AlgebraSE3<T> SE3<T>::log() const {
    SO3<T> rotation = rotation_part();
    AlgebraSO3<T> omega = rotation.log();

    T theta = sqrt(omega.as_vector().transpose() * omega.as_vector());

    T coeff2 = 0.0;
    if (theta < SMALL_ANGLE_THRESHOLD) {
        // Use Taylor series approximations to avoid numerical instability.
        coeff2 = 1 - (theta * theta) / 4;
    } else {
        T a = sin(theta) / theta;
        T b = (1 - cos(theta)) / (theta * theta);

        coeff2 = (1 / (theta * theta)) * (1 - (a / (2 * b)));
    }

    Eigen::Matrix<T, 3, 3> v_inv =
        Eigen::Matrix<T, 3, 3>::Identity() - (0.5 * omega.as_skewsym_matrix()) +
        coeff2 * omega.as_skewsym_matrix() * omega.as_skewsym_matrix();

    Eigen::Matrix<T, 3, 1> rho = v_inv * translation_part();

    return AlgebraSE3<T>(rho, omega.as_vector());
}

template <class T>
Eigen::Matrix<T, 3, 1> SE3<T>::translation_part() const {
    return as_matrix().block(0, 3, 3, 1);
}

template <class T>
SE3<T> SE3<T>::perturbate(T iso_covariance) const {
    return perturbate(iso_covariance, iso_covariance);
}

template <class T>
SE3<T> SE3<T>::perturbate(T iso_covariance, T iso_covariance_rot) const {
    AlgebraSE3<T> perturbation =
        AlgebraSE3<T>::random(iso_covariance, iso_covariance_rot);
    return perturbation.exp() * *this;
}

template <class T>
SE3<T> SE3<T>::perturbate(const Eigen::Matrix<T, 6, 6>& covariance) const {
    AlgebraSE3<T> perturbation = AlgebraSE3<T>::random(covariance);
    return perturbation.exp() * *this;
}

template <class T>
Eigen::Matrix<T, 6, 6> SE3<T>::adjoint() const {
    Eigen::Matrix<T,6,6> adjoint = Eigen::Matrix<T,6,6>::Zero();
    Eigen::Matrix<T,3,3> rotation = rotation_part().as_matrix();
    adjoint.block(0,0,3,3) = rotation;
    adjoint.block(3,3,3,3) = rotation;

    Eigen::Matrix<T,3,1> trans = translation_part();
    adjoint.block(0,3,3,3) = skew_sym(trans) * rotation;

    return adjoint;
}

template <class T>
SE3<T> SE3<T>::inv() const {
    Eigen::Matrix<T, 4, 4> inverse;
    inverse.setIdentity();

    Eigen::Matrix<T, 3, 3> rotation_matrix = rotation_part().as_matrix();
    inverse.block(0, 0, 3, 3) = rotation_matrix.transpose();

    Eigen::Matrix<T, 3, 1> t = translation_part();
    inverse.block(0, 3, 3, 1) = -rotation_matrix.transpose() * t;

    return SE3<T>(inverse);
}

template <class T>
SE3<T>& SE3<T>::operator*=(const SE3<T>& rhs) {
    auto m = as_matrix() * rhs.as_matrix();

    for(auto i = 0; i < 4; ++i) {
        for(auto j = 0; j < 4; ++j) {
            values[i][j] = m(i,j);
        }
    }

    return *this;
}

template <class T>
SE3<T>& SE3<T>::operator=(const SE3<T>& rhs) {
    for(auto i = 0; i < 4; ++i) {
        for(auto j = 0; j < 4; ++j) {
            values[i][j] = rhs.values[i][j];
        }
    }
    return *this;
}

template <class T>
SE3<T>& SE3<T>::operator=(SE3<T>&& rhs) {
    for(auto i = 0; i < 4; ++i) {
        for(auto j = 0; j < 4; ++j) {
            values[i][j] = rhs.values[i][j];
        }
    }
    return *this;
}

template <class T>
Eigen::Matrix<T, 4, 4> SE3<T>::as_matrix() const {
    Eigen::Matrix<T,4,4> m;

    for(auto i = 0; i < 4; ++i) {
        for(auto j = 0; j < 4; ++j) {
            m(i,j) = values[i][j];
        }
    }

    return m;
}

template <class T>
Eigen::Transform<T, 3, Eigen::Affine> SE3<T>::as_transform() const {
    Eigen::Transform<T, 3, Eigen::Affine, Eigen::ColMajor> t;
    t.matrix() = as_matrix();
    return t;
}

}

#endif
