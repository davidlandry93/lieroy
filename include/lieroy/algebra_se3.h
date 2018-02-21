#ifndef LIEROY_ALGEBRA_SE3_H
#define LIEROY_ALGEBRA_SE3_H

#include <array>

#include <Eigen/Core>

#include "lieroy/algebra_so3.h"

namespace lieroy {
  template<class T>
  class SE3;

  template<class T>
  class AlgebraSE3 {
    const T SMALL_ANGLE_THRESHOLD = 1e-6;
  public:
    AlgebraSE3();
    AlgebraSE3(const Eigen::Matrix<T,3,1>& rho, const Eigen::Matrix<T,3,1>& omega);
    AlgebraSE3(const std::array<T,6>& v);
    AlgebraSE3(const Eigen::Matrix<T,6,1>& v);
    AlgebraSE3(const Eigen::Matrix<T,1,6>& v);

    Eigen::Matrix<T,3,1> translation_part() const;
    AlgebraSO3<T> rotation_part() const;
    Eigen::Matrix<T,6,1> as_vector() const;
    Eigen::Matrix<T,6,1>& as_vector();
    void stream_to(std::ostream& os) const;
    SE3<T> exp() const;
    AlgebraSE3<T>& operator=(const AlgebraSE3<T>& rhs);

    static Eigen::Matrix<T,6,6> covariance_matrix_from_variance(T trans_var, T rot_var);
    static AlgebraSE3<T> random(T isotropic_covariance);
    static AlgebraSE3<T> random(T isotropic_covariance_translation,
                                T isotropic_covariance_rot);
    static AlgebraSE3<T> random(const Eigen::Matrix<T,6,6>& covariance);
  private:
    Eigen::Matrix<T,6,1> vector;
  };

}

#endif
