
#include <exception>

#include <boost/format.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <Eigen/Core>

#include "pylie/algebra_se3.hpp"
#include "pylie/se3.hpp"

namespace p = boost::python;
namespace np = boost::python::numpy;

template <typename T, int I, int J>
Eigen::Matrix<T,I,J> ndarray_to_eigen_matrix(const np::ndarray& np_m) {
  Eigen::Matrix<T,I,J> eigen_m;

  for(auto i = 0; i < I; i++) {
    for(auto j = 0; j < J; j++) {
      eigen_m(i,j) = p::extract<T>(np_m[i][j]);
    }
  }

  return eigen_m;
}

template <typename T, int I, int J>
np::ndarray eigen_matrix_to_ndarray(const Eigen::Matrix<T,I,J>& eigen_m) {
  p::list matrix;

  for(auto i = 0; i < I; i++) {
    p::list row;
    for(auto j = 0; j < J; j++) {
      row.append(eigen_m(i,j));
    }
    matrix.append(row);
  }

  return np::array(matrix);
}

template <typename T, int N>
Eigen::Matrix<T,6,1> nd_array_to_eigen_vector(const np::ndarray& nd_m) {
  if(nd_m.get_nd() != 1) {
    std::runtime_error("Input has wrong shape");
  }

  Eigen::Matrix<T,6,1> eigen_m;

  for(auto i = 0; i < N; i++) {
    eigen_m(i) = p::extract<T>(nd_m[i]);
  }

  return eigen_m;
}

template <typename T, int N>
np::ndarray eigen_vector_to_array(const Eigen::Matrix<T,N,1>& eigen_m) {
  p::list array;

  for(auto i = 0; i < N; i++) {
    array.append(eigen_m(i));
  }

  return np::array(array);
}

template <typename T>
np::ndarray se3_log(const np::ndarray& m) {
  if(m.shape(0) != 4 || m.shape(1) != 4) {
    throw std::runtime_error("Input has wrong shape");
  }

  Eigen::Matrix<T,4,4> eigen_m = ndarray_to_eigen_matrix<T,4,4>(m);
  pylie::AlgebraSE3<T> log_of_lie = pylie::SE3<T>(eigen_m).log();

  return eigen_vector_to_array<T,6>(log_of_lie.as_vector());
}

template <typename T>
np::ndarray se3_exp(const np::ndarray& m) {
  Eigen::Matrix<T,6,1> eigen_m;

  if(m.get_nd() == 2 && m.shape(0) == 6 && m.shape(1) == 1) {
    eigen_m = ndarray_to_eigen_matrix<T,6,1>(m);
  } else if (m.get_nd() == 2 && m.shape(0) == 1 && m.shape(1) == 6) {
    eigen_m = ndarray_to_eigen_matrix<T,1,6>(m).transpose();
  } else if (m.get_nd() == 1 && m.shape(0) == 6) {
    eigen_m = nd_array_to_eigen_vector<T,6>(m);
  } else {
    throw std::runtime_error("Input has wrong shape");
  }

  pylie::AlgebraSE3<T> log_of_lie(eigen_m);
  return eigen_matrix_to_ndarray<T,4,4>(log_of_lie.exp().as_matrix());
}

BOOST_PYTHON_MODULE(pylie) {
  np::initialize();

  p::def("se3_log", se3_log<double>, "Compute the Lie algebra counterpart of a SE3 member. Takes a 4x4 Numpy matrix as input.");
  p::def("se3_log", se3_log<float>, "Compute the Lie algebra counterpart of a SE3 member. Takes a 4x4 Numpy matrix as input.");
  p::def("se3_exp", se3_exp<double>, "Compute the Lie Group counterpart of a member of the Lie algebra of SE3. Takes a 6-vector as input.");
  p::def("se3_exp", se3_exp<float>,  "Compute the Lie Group counterpart of a member of the Lie algebra of SE3. Takes a 6-vector as input.");
}
