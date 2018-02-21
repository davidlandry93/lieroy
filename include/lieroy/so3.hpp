#ifndef LIEROY_SO3_HPP
#define LIEROY_SO3_HPP

#include <cmath>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include "lieroy/so3.h"

namespace lieroy {
    template <class T>
    SO3<T>::SO3() :
        matrix() {
        matrix.setIdentity();
    }

    template <class T>
    SO3<T>::SO3(const std::array<T,9>& m) : matrix(m) {
        for(int i=0; i < 3; i++) {
            for(int j=0; j < 3; j++) {
                matrix(i,j) = m[i*3+j];
            }
        }
    }

    template <class T>
    SO3<T>::SO3(const Eigen::Matrix<T,3,3>& m) :
        matrix(m) {
    }

    template <class T>
    AlgebraSO3<T> SO3<T>::log() const {
        // Get angle
        const double phi_ba = acos(0.5*(matrix.trace()-1.0));
        const double sinphi_ba = sin(phi_ba);

        if (fabs(sinphi_ba) > 1e-9) {

            // General case, angle is NOT near 0, pi, or 2*pi
            Eigen::Matrix<T,3,1> axis;
            axis << matrix(2,1) - matrix(1,2),
                    matrix(0,2) - matrix(2,0),
                    matrix(1,0) - matrix(0,1);
            return AlgebraSO3<T>((0.5*phi_ba/sinphi_ba)*axis);

        } else if (fabs(phi_ba) > 1e-9) {

            // Angle is near pi or 2*pi
            // ** Note with this method we do not know the sign of 'phi', however since we know phi is
            //    close to pi or 2*pi, the sign is unimportant..

            // Find the eigenvalues and eigenvectors
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,3,3> > eigenSolver(matrix);

            // Try each eigenvalue
            for (int i = 0; i < 3; i++) {
                // Check if eigen value is near +1.0
                if ( fabs(eigenSolver.eigenvalues()[i] - 1.0) < 1e-4 ) {

                    // Get corresponding angle-axis
                    Eigen::Matrix<T,3,1> aaxis_ba = phi_ba*eigenSolver.eigenvectors().col(i);
                    return AlgebraSO3<T>(aaxis_ba);
                }
            }

            // Runtime error
            throw std::runtime_error("so3 logarithmic map failed to find an axis-angle, "
                                     "angle was near pi, or 2*pi, but no eigenvalues were near 1");
        } else {

            // Angle is near zero
            return AlgebraSO3<T>(Eigen::Matrix<T,3,1>::Zero());
        }
    }

    template <class T>
    AlgebraSO3<T> SO3<T>::log_rodrigues() const {
        T cos_theta = (matrix.trace() - 1.0) * 0.5;
        Eigen::Matrix<T,3,1> omega;
        omega << matrix(2,1), matrix(0,2), matrix(1,0);

        T theta_2 = omega.transpose() * omega;

        T factor = 0.0;
        if(cos_theta > (T) 0.999856) {
            // Theta was small. Use taylor approximation.
            factor = (T) 1.0 +
                theta_2 * ((T) 1/6.0) +
                theta_2 * ((T) 3/40.0) +
                theta_2 * ((T) 5/112.0);
        } else if(cos_theta > (T) -0.99) {
            // Regular angles.
            T theta = acos(cos_theta);
            T st = sqrt(theta_2);
            factor = theta / st;
        } else {
            // Angle is near pi.
            // Use Taylor approximations.
            T theta = acos(cos_theta);
            T st = sqrt(theta_2);
            factor = theta / st;
        }

        omega = factor * omega;
        return AlgebraSO3<T>(omega(0), omega(1), omega(2));
    }

    template <class T>
    Eigen::Matrix<T,3,3> SO3<T>::as_matrix() const {
        return matrix;
    }

    template <class T>
    SO3<T> SO3<T>::from_rpy(T roll, T pitch, T yaw) {
        Eigen::Transform<T, 3, Eigen::Affine> transform;
        transform.setIdentity();

        transform.rotate(Eigen::AngleAxis<T>(roll, Eigen::Vector3d::UnitX().cast<T>()));
        transform.rotate(Eigen::AngleAxis<T>(pitch, Eigen::Vector3d::UnitY().cast<T>()));
        transform.rotate(Eigen::AngleAxis<T>(yaw, Eigen::Vector3d::UnitZ().cast<T>()));

        return SO3<T>(transform.rotation().matrix());
    }
}

#endif
