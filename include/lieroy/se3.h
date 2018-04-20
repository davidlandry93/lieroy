#ifndef LIEROY_SE3_H
#define LIEROY_SE3_H

#include <array>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "so3.hpp"
#include "algebra_se3.hpp"

namespace lieroy {
    template<class T>
    class SE3 {
        const T SMALL_ANGLE_THRESHOLD = 1e-2;
    public:
        typedef Eigen::Matrix<T, 6, 6> Covariance;

        SE3();
        SE3(const std::array<T, 16>& m);
        SE3(const Eigen::Matrix<T,4,4>& m);
        SE3(const SE3<T>& other);
        SE3(const Eigen::Matrix<T,3,1>& translation, const SO3<T>& rotation);

        static SE3<T> identity() {
            return SE3<T>();
        }

        Eigen::Matrix<T,4,4> as_matrix() const;
        Eigen::Transform<T,3,Eigen::Affine> as_transform() const;

        SO3<T> rotation_part() const;
        Eigen::Matrix<T,3,1> translation_part() const;
        AlgebraSE3<T> log() const;
        SE3<T> inv() const;
        SE3<T> perturbate(T iso_covariance) const;
        SE3<T> perturbate(T iso_covariance,
                          T iso_covariance_rot) const;
        SE3<T> perturbate(const Eigen::Matrix<T,6,6>& covariance) const;

        Eigen::Matrix<T,6,6> adjoint() const;

        T& operator()(int i, int j);
        const T& operator()(int i, int j) const;
        SE3<T>& operator*=(const SE3<T>& rhs);
        SE3<T>& operator=(SE3<T>&& rhs);
        SE3<T>& operator=(const SE3<T>& rhs);

        void stream_to(std::ostream& os) const;

    private:
        Eigen::Matrix<T,4,4> matrix;
    };

    template <typename T>
    SE3<T> operator*(SE3<T> lhs, const SE3<T>& rhs) {
        return lhs *= rhs;
    }
}

#endif
