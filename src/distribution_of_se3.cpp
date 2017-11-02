
#include <iostream>
#include <vector>

#include "json.hpp"
#include "util.hpp"

#include "pylie/se3.hpp"
#include "pylie/se3_gaussian_distribution.hpp"

using json = nlohmann::json;
using namespace pylie;

int main(int argc, char** argv) {
  json input;
  std::cin >> input;

  std::vector<SE3<double>> sample;
  sample.reserve(input.size());
  for(const auto& elt : input) {
    SE3<double> t = read_json_transform<double>(elt);
    sample.push_back(t);

    std::cout << t << '\n';
  }

  auto distribution = SE3GaussianDistribution<double>::from_sample(sample);

  json output;
  output["mean"] = json_of_matrix(distribution.mean.as_matrix());
  output["covariance"] = json_of_matrix(distribution.covariance);

  std::cout << output << std::endl;

  return 0;
}
