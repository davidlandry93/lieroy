
#include <Eigen/Core>

#include "json.hpp"

#include "lieroy/se3.hpp"


namespace lieroy {
  template <typename T> SE3<T> read_json_transform(const nlohmann::json& input);
  template <typename T, int R, int C> nlohmann::json json_of_matrix(const Eigen::Matrix<T,R,C>& m);
}
