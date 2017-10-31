#ifndef PYLIE_ALGEBRA_SE3_HPP
#define PYLIE_ALGEBRA_SE3_HPP

#include <chrono>
#include <random>

#include "algebra_se3.h"
#include "normal_random_variable.hpp"

namespace pylie {
  template <class T>
  AlgebraSE3<T>::AlgebraSE3() :
    vector() {}

  template <class T>
  AlgebraSE3<T>::AlgebraSE3(const std::array<T, 6>& v) :
    vector () {
    for(int i = 0; i < 6; i++) {
      vector(i) = v[i];
    }
  }

  template <class T>
  AlgebraSE3<T>::AlgebraSE3(const Eigen::Matrix<T,6,1>& v) :
    vector(v) {}

  template <class T>
  AlgebraSE3<T>::AlgebraSE3(const Eigen::Matrix<T,1,6>& v) :
    vector(v.transpose()) {}

  template <class T>
  AlgebraSE3<T>::AlgebraSE3(const Eigen::Matrix<T,3,1>& rho, const Eigen::Matrix<T,3,1>& omega) :
    vector() {
    vector << rho(0), rho(1), rho(2), omega(0), omega(1), omega(2);
  }

  template <class T>
  void AlgebraSE3<T>::stream_to(std::ostream& os) const {
    os << vector;
  }

  template <class T>
  Eigen::Matrix<T,6,1> AlgebraSE3<T>::as_vector() const {
    return vector;
  }

  template <class T>
  Eigen::Matrix<T,6,1>& AlgebraSE3<T>::as_vector() {
    return vector;
  }

  template <class T>
  AlgebraSE3<T>& AlgebraSE3<T>::operator=(const AlgebraSE3<T>& rhs) {
    vector = rhs.vector;
    return *this;
  }

  template <class T>
  inline std::ostream& operator<<(std::ostream& os, const AlgebraSE3<T>& t) {
    t.stream_to(os);
    return os;
  }

  template <class T>
  SE3<T> AlgebraSE3<T>::exp() const {
    auto so3_part = rotation_part();
    auto omega = so3_part.as_vector();
    auto translation_part = vector.block(0,0,3,1);

    T theta = sqrt(omega.transpose() * omega);

    T b = 0.0, c = 0.0;
    if(theta < SMALL_ANGLE_THRESHOLD) {
      // Use Taylor series to avoid numerical instability.
      b = 0.5 - (theta*theta)/24.0;
      c = 1.0 / 6.0 - (theta*theta)/120.0;
    } else {
      b = (1 - cos(theta)) / (theta*theta);
      c = (1 - sin(theta) / theta) / (theta*theta);
    }

    auto skewsym = so3_part.as_skewsym_matrix();
    Eigen::Matrix<T,3,3> v;
    v.setIdentity();
    v += b * skewsym + c * skewsym * skewsym;

    Eigen::Matrix<T,4,4> m;
    m.setIdentity();
    m.block(0,0,3,3) = so3_part.exp().as_matrix();
    m.block(0,3,3,1) = v * translation_part;

    return SE3<T>(m);
  }

  template <class T>
  AlgebraSO3<T> AlgebraSE3<T>::rotation_part() const {
    return AlgebraSO3<T>(vector.block(3,0,3,1));
  }

  template <class T>
  Eigen::Matrix<T,3,1> AlgebraSE3<T>::translation_part() const {
    return vector.block(0,0,3,1);
  }

  template <class T>
  Eigen::Matrix<T,6,6> AlgebraSE3<T>::covariance_matrix_from_variance(T trans_var, T rot_var) {
    Eigen::Matrix<T,6,6> covariance = Eigen::Matrix<T,6,6>::Identity();

    for(int i = 0; i < 3; i++) {
      covariance(i,i) = trans_var;
      covariance(i+3, i+3) = rot_var;
    }

    return covariance;
  }

  template <class T>
  AlgebraSE3<T> AlgebraSE3<T>::random(T isotropic_covariance) {
    return random(isotropic_covariance, isotropic_covariance);
  }

  template <class T>
  AlgebraSE3<T> AlgebraSE3<T>::random(T isotropic_covariance,
                                      T isotropic_covariance_rot) {
    Eigen::Matrix<T,6,6> covariance = covariance_matrix_from_variance(isotropic_covariance,
                                                                      isotropic_covariance_rot);

    return random(covariance);
  }

  template <class T>
  AlgebraSE3<T> AlgebraSE3<T>::random(const Eigen::Matrix<T,6,6>& covariance) {
    NormalRandomVariable<T, 6> normal_variable(covariance);
    auto vec = normal_variable();

    return AlgebraSE3<T>(vec);
  }
}
#endif
