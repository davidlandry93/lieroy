#ifndef LIEROY_SE3_HPP
#define LIEROY_SE3_HPP

#include <algorithm>
#include <iomanip>
#include <random>
#include <stdexcept>

#include "lieroy/algebra_so3.hpp"
#include "lieroy/normal_random_variable.hpp"
#include "lieroy/se3.h"

namespace lieroy {
    template <class T>
    SE3<T>::SE3() :
        matrix(){
        matrix.setIdentity();
    }

    template <class T>
    SE3<T>::SE3(const std::array<T, 16>& values) :
        matrix() {
        for(int i=0; i < 4; i++) {
            for(int j=0; j < 4; j++) {
                matrix(i,j) = values[i*4 + j];
            }
        }
    }

    template <class T>
    SE3<T>::SE3(const Eigen::Matrix<T,4,4>& m) :
        matrix(m) {}

    template <class T>
    SE3<T>::SE3(const SE3<T>& other) :
        matrix(other.matrix){}

    template <class T>
    SE3<T>::SE3(const Eigen::Matrix<T,3,1>& translation, const SO3<T>& rotation) :
        matrix() {
        matrix.setIdentity();
        matrix.block(0,0,3,3) = rotation.as_matrix();
        matrix.block(0,3,3,1) = translation;
    }

    template <class T>
    T& SE3<T>::operator()(int i, int j) {
        return matrix(i,j);
    }

    template <class T>
    const T& SE3<T>::operator()(int i, int j) const {
        return matrix(i,j);
    }

    template <class T>
    void SE3<T>::stream_to(std::ostream& os) const {
        os << matrix;
    }

    template <class T>
    inline std::ostream& operator<<(std::ostream& os, const SE3<T>& t) {
        t.stream_to(os);
        return os;
    }

    template <class T>
    SO3<T> SE3<T>::rotation_part() const {
        return SO3<T>(matrix.block(0,0,3,3));
    }

    template <class T>
    AlgebraSE3<T> SE3<T>::log() const {
        SO3<T> rotation(rotation_part());
        AlgebraSO3<T> omega = rotation.log();

        T theta = sqrt(omega.as_vector().transpose() * omega.as_vector());

        T coeff2 = 0.0;
        if(theta < SMALL_ANGLE_THRESHOLD) {
            // Use Taylor series approximations to avoid numerical instability.
            coeff2 = 1 - (theta*theta)/4;
        } else {
            T a = sin(theta) / theta;
            T b = (1 - cos(theta)) / (theta * theta);

            coeff2 = (1 / (theta * theta)) * (1 - (a / (2*b)));
        }

        Eigen::Matrix<T,3,3> v_inv =
            Eigen::Matrix<T,3,3>::Identity() -
            (0.5 * omega.as_skewsym_matrix()) +
            coeff2 * omega.as_skewsym_matrix() * omega.as_skewsym_matrix();

        Eigen::Matrix<T,3,1> rho = v_inv * translation_part();

        return AlgebraSE3<T>(rho, omega.as_vector());
    }

    template <class T>
    Eigen::Matrix<T,3,1> SE3<T>::translation_part() const {
        return matrix.block(0,3,3,1);
    }

    template <class T>
    SE3<T> SE3<T>::perturbate(T iso_covariance) const {
        return perturbate(iso_covariance, iso_covariance);
    }

    template <class T>
    SE3<T> SE3<T>::perturbate(T iso_covariance,
                              T iso_covariance_rot) const {
        AlgebraSE3<T> perturbation = AlgebraSE3<T>::random(iso_covariance, iso_covariance_rot);
        return perturbation.exp() * *this;
    }


    template <class T>
    SE3<T> SE3<T>::perturbate(const Eigen::Matrix<T,6,6>& covariance) const {
        AlgebraSE3<T> perturbation = AlgebraSE3<T>::random(covariance);
        return perturbation.exp() * *this;
    }

    template <class T>
    SE3<T> SE3<T>::inv() const {
        Eigen::Matrix<T,4,4> inverse;
        inverse.setIdentity();

        Eigen::Matrix<T,3,3> rotation_matrix = rotation_part().as_matrix();
        inverse.block(0,0,3,3) = rotation_matrix.transpose();

        Eigen::Matrix<T,3,1> t = matrix.block(0,3,3,1);
        inverse.block(0,3,3,1) = -rotation_matrix.transpose() * t;

        return SE3<T>(inverse);
    }

    template <class T>
    SE3<T>& SE3<T>::operator*=(const SE3<T>& rhs) {
        matrix = matrix * rhs.matrix;
        return *this;
    }

    template <class T>
    SE3<T>& SE3<T>::operator=(const SE3<T>& rhs) {
        matrix = rhs.matrix;
    }

    template <class T>
    SE3<T>& SE3<T>::operator=(SE3<T>&& rhs) {
        matrix = rhs.matrix;
    }

    template <class T>
    Eigen::Matrix<T,4,4> SE3<T>::as_matrix() const {
        return matrix;
    }

    template <class T>
    Eigen::Transform<T,3,Eigen::Affine> SE3<T>::as_transform() const {
        Eigen::Transform<T,3,Eigen::Affine,Eigen::ColMajor> t;
        t.matrix() = matrix;
        return t;
    }
}

#endif
