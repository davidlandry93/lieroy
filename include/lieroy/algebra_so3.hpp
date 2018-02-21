#ifndef LIEROY_ALGEBRA_SO3_HPP
#define LIEROY_ALGEBRA_SO3_HPP

#include "algebra_so3.h"

namespace lieroy {
    template <class T>
    Eigen::Matrix<T,3,3> skew_symetric_matrix(const Eigen::Matrix<T,3,1>& vector) {
        Eigen::Matrix<T,3,3> m;
        m << 0.0, -vector(2), vector(1),
            vector(2), 0.0, -vector(0),
            -vector(1), vector(0), 0.0;

        return m;
    }

    template <class T>
    AlgebraSO3<T>::AlgebraSO3() : vector() {}

    template <class T>
    AlgebraSO3<T>::AlgebraSO3(T w1, T w2, T w3) :
        vector() {
        vector(0,0) = w1;
        vector(1,0) = w2;
        vector(2,0) = w3;
    }

    template <class T>
    AlgebraSO3<T>::AlgebraSO3(const std::array<T, 3>& v) :
        vector() {
        for(int i = 0; i < 3; i++) {
            vector(i,0) = v[i];
        }
    }

    template <class T>
    AlgebraSO3<T>::AlgebraSO3(const Eigen::Matrix<T,3,1>& v) :
        vector(v) {}

    template <class T>
    Eigen::Matrix<T,3,1> AlgebraSO3<T>::as_vector() const {
        return vector;
    }


    template <class T>
    Eigen::Matrix<T,3,3> AlgebraSO3<T>::as_skewsym_matrix() const {
        return skew_symetric_matrix(vector);
    }
    template <class T>
    SO3<T> AlgebraSO3<T>::exp() const {
        T theta = sqrt(vector.transpose() * vector);

        if(theta < 1e-12) {
            return SO3<T>(Eigen::Matrix<T,3,3>::Identity());
        }

        // Analytic solution from lgmath
        Eigen::Matrix<T,3,1> axis = vector / theta;
        const T sin_theta = sin(theta);
        const T cos_theta = cos(theta);

        const Eigen::Matrix<T,3,3> result = cos_theta * Eigen::Matrix<T,3,3>::Identity() +
                                            (1.0 - cos_theta) * axis * axis.transpose() +
                                            sin_theta * skew_symetric_matrix(axis);
        return SO3<T>(result);
    }

    template <class T>
    SO3<T> AlgebraSO3<T>::exp_rodrigues() const {
        T coeff1 = 0.0, coeff2 = 0.0;

        T theta = sqrt(vector.transpose() * vector);
        if (theta == 0.0) {
            coeff1 = 1.0;
            coeff2 = 0.5;
        } else {
            coeff1 = sin(theta) / theta;
            coeff2 = (1 - cos(theta)) / (theta * theta);
        }

        Eigen::Matrix<T,3,3> m;
        m.setIdentity();

        auto skewsym = as_skewsym_matrix();
        return SO3<T>(m + coeff1 * skewsym + coeff2 * skewsym * skewsym);
    }
}
#endif
