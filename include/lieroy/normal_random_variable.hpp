#ifndef LIEROY_NORMAL_RANDOM_VARIABLE_H
#define LIEROY_NORMAL_RANDOM_VARIABLE_H

#include <iostream>
#include <random>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace lieroy {

    // http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
    template <typename T, int N>
    struct NormalRandomVariable
    {
        NormalRandomVariable(Eigen::Matrix<T, N, N> const& covar)
            : NormalRandomVariable(Eigen::Matrix<T, N, 1>::Zero(), covar) {}

        NormalRandomVariable(Eigen::Matrix<T,N,1> const& mean, Eigen::Matrix<T,N,N> const& covar)
            : mean(mean) {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,N,N>> eigenSolver(covar);

            Eigen::Matrix<T,N,1> eigenvalues = eigenSolver.eigenvalues();

            auto to_decompose = covar;
            if(eigenvalues(0) < 0.0) {
                // Some eigenvalues are negative.
                // The matrix was not exactly PSD.
                // Correct it before we continue.
                auto corrected_matrix = corrected_covariance_matrix(covar, eigenvalues);
                eigenSolver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,N,N>>(corrected_matrix);
            }

            transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
        }

        Eigen::Matrix<T,N,1> mean;
        Eigen::Matrix<T,N,N> transform;

        Eigen::Matrix<T,N,1> operator()() const {
            static std::mt19937 gen{ std::random_device{}() };
            static std::normal_distribution<T> dist;
            return mean + transform * Eigen::Matrix<T,N,1>{ N }.unaryExpr([&](T x) { return dist(gen); });
        }

    private:
        Eigen::Matrix<T,N,N> corrected_covariance_matrix(Eigen::Matrix<T,N,N> const& covariance, Eigen::Matrix<T,N,1> const& eigenvalues) {
            // Correct covariance matrices that are not exactly positive semidefinite but close due to numerical errors.
            //TODO: Use a less evil method to correct a covariance matrix.
            // There are more principled ways to do so but I did not find an implementation in c++.
            // The current solution is drawn form https://www.value-at-risk.net/non-positive-definite-covariance-matrices/
            // Some better solutions are quoted here https://quant.stackexchange.com/questions/2074/what-is-the-best-way-to-fix-a-covariance-matrix-that-is-not-positive-semi-defi

            Eigen::Matrix<T,N,N> corrected;
            corrected.setIdentity();
            corrected *= -2 * eigenvalues(0); // The smallest eigenvalue is the first one.
            corrected += covariance;

            std::cerr << "Covariance matrix was not positive semi definite. Correcting matrix with a " << eigenvalues(0) << " factor.";

            return corrected;
        }
    };
}

#endif
