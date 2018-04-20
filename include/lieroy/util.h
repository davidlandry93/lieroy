#ifndef LIEROY_UTIL_H
#define LIEROY_UTIL_H

#include <Eigen/Core>

template <typename T>
static Eigen::Matrix<T, 3, 1> sample_uniform_from_sphere(T distribution_radius);


template <typename T>
static Eigen::Matrix<T, 3, 3> skew_sym(const Eigen::Matrix<T,3,1>& m);

#endif
