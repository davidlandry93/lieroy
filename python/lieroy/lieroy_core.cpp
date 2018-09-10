
#include <exception>

#include <boost/format.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/tuple.hpp>
#include <Eigen/Core>

#include "lieroy/algebra_se3.hpp"
#include "lieroy/se3.hpp"
#include "lieroy/se3_gaussian_distribution.hpp"
#include "lieroy/se3_uniform_distribution.hpp"

using namespace lieroy;

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
  lieroy::AlgebraSE3<T> log_of_lie = lieroy::SE3<T>(eigen_m).log();

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

  lieroy::AlgebraSE3<T> log_of_lie(eigen_m);
  return eigen_matrix_to_ndarray<T,4,4>(log_of_lie.exp().as_matrix());
}

template <typename T>
p::tuple se3_gaussian_distribution_of_sample(const np::ndarray& m) {
  if(m.get_nd() != 3 || m.shape(1) != 4 || m.shape(2) != 4) {
    throw std::runtime_error("Input has wrong shape");
  }

  std::vector<lieroy::SE3<double>> transformations;
  for(auto i = 0; i < m.shape(0); ++i) {
    Eigen::Matrix<T,4,4> eigen_m;

    for(auto j = 0; j < m.shape(1); ++j) {
      for(auto k = 0; k < m.shape(2); ++k) {
        eigen_m(j,k) = p::extract<T>(m[i][j][k]);
      }
    }

    transformations.push_back(lieroy::SE3<double>(eigen_m));
  }

  auto distribution = lieroy::SE3GaussianDistribution<double>::from_sample(transformations);

  return p::make_tuple(eigen_matrix_to_ndarray<T,4,4>(distribution.mean.as_matrix()),
                       eigen_matrix_to_ndarray(distribution.covariance));
}


template <typename T>
p::list se3_sample_normal_distribution(const np::ndarray& mean, const np::ndarray& covariance, const int& n_samples) {
    auto eigen_mean = ndarray_to_eigen_matrix<T, 4, 4>(mean);
    auto eigen_covariance = ndarray_to_eigen_matrix<T, 6, 6>(covariance);
    auto distribution = SE3GaussianDistribution<T>(eigen_mean, eigen_covariance);

    p::list samples;

    for(auto i = 0; i < n_samples; i++) {
        SE3<T> perturbated, perturbation;
        std::tie(perturbated, perturbation) = distribution.sample_with_perturbation();

        samples.append(eigen_matrix_to_ndarray(perturbated.as_matrix()));
    }

    return samples;
}

template <typename T>
p::list se3_sample_uniform_distribution(const np::ndarray& mean, const T& range_translation, const T& range_rotation, const int& n_samples) {
    auto eigen_mean = ndarray_to_eigen_matrix<T, 4, 4>(mean);

    std::cout << "Ranges: " << range_translation << range_rotation << std::endl;

    auto distribution = SE3UniformDistribution<T>(eigen_mean, range_translation, range_rotation);

    p::list samples;
    for(auto i = 0; i < n_samples; i++) {
        auto perturbated = distribution.sample();
        samples.append(eigen_matrix_to_ndarray(perturbated.as_matrix()));
    }

    return samples;
}

template <typename T>
np::ndarray se3_adjoint(const np::ndarray& t) {
    auto eigen_t = ndarray_to_eigen_matrix<T,4,4>(t);
    SE3<double> se3_t(eigen_t);

    Eigen::Matrix<T,6,6> adjoint = se3_t.adjoint();

    return eigen_matrix_to_ndarray(adjoint);
}


template <typename T>
p::tuple combine_poses(const np::ndarray& t1, const np::ndarray& sigma1, const np::ndarray& t2, const np::ndarray& sigma2) {
    auto t1_eigen = ndarray_to_eigen_matrix<T,4,4>(t1);
    auto t2_eigen = ndarray_to_eigen_matrix<T,4,4>(t2);
    SE3<T> t1_group(t1_eigen);
    SE3<T> t2_group(t2_eigen);

    auto sigma1_eigen = ndarray_to_eigen_matrix<T,6,6>(sigma1);
    auto sigma2_eigen = ndarray_to_eigen_matrix<T,6,6>(sigma2);

    SE3<T> mean;
    Eigen::Matrix<T,6,6> covariance;
    std::tie(mean, covariance) = compound_poses(t1_group, sigma1_eigen, t2_group, sigma2_eigen);

    auto mean_np = eigen_matrix_to_ndarray<T,4,4>(mean.as_matrix());
    auto covariance_np = eigen_matrix_to_ndarray<T,6,6>(covariance);

    return p::make_tuple(mean_np, covariance_np);
}


BOOST_PYTHON_MODULE(lieroy_core) {
  np::initialize();

  p::def("se3_log", se3_log<double>, "Compute the Lie algebra counterpart of a SE3 member. Takes a 4x4 Numpy matrix as input.");
  p::def("se3_log", se3_log<float>, "Compute the Lie algebra counterpart of a SE3 member. Takes a 4x4 Numpy matrix as input.");
  p::def("se3_exp", se3_exp<double>, "Compute the Lie Group counterpart of a member of the Lie algebra of SE3. Takes a 6-vector as input.");
  p::def("se3_exp", se3_exp<float>,  "Compute the Lie Group counterpart of a member of the Lie algebra of SE3. Takes a 6-vector as input.");
  p::def("se3_gaussian_distribution_of_sample", se3_gaussian_distribution_of_sample<double>, "Compute a gaussian distribution from a collection of SE3 transforms");
  p::def("se3_sample_normal_distribution", se3_sample_normal_distribution<double>, "Create a collection of SE3 sampled from a normal distribution.");
  p::def("se3_sample_uniform_distribution", se3_sample_uniform_distribution<double>, "Create a collection of SE3 sampled from a uniform distribution.");
  p::def("se3_adjoint", se3_adjoint<double>, "Compute the adjoint representation of a SE3 group element.");
  p::def("se3_compound_poses", combine_poses<double>, "Combine se3 poses and their uncertainty.");
  p::def("se3_compound_poses", combine_poses<float>, "Combine se3 poses and their uncertainty.");
}
