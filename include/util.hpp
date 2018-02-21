#ifndef LIEROY_UTIL_HPP
#define LIEROY_UTIL_HPP

#include "util.h"

using json = nlohmann::json;

namespace lieroy {
  template <typename T>
  SE3<T> read_json_transform(const json& input) {
    Eigen::Matrix<T, 4, 4> m;

    for(auto i = 0; i < 4; i++) {
      for(auto j = 0; j < 4; j++) {
        m(i,j) = input[i][j];
      }
    }

    return SE3<T>(m);
  }

  template <typename T, int R, int C>
  json json_of_matrix(const Eigen::Matrix<T,R,C>& m) {
    json output;

    for(auto i = 0; i < R; i++) {
      json row;
      for (auto j = 0; j < C; j++) {
        row.push_back(m(i,j));
      }
      output.push_back(row);
    }

    return output;
  }
}


#endif
