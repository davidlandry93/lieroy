
#include <cmath>
#include <gtest/gtest.h>

#include "pylie/algebra_so3.hpp"
#include "pylie/so3.hpp"

using namespace pylie;

TEST(SO3, CycleConsistency) {

  const double MIN_ANGLE = -M_PI;
  const double MAX_ANGLE = +M_PI;
  const double STEP_ANGLE = M_PI / 32.0;

  const double MIN_RADIUS = 0;
  const double MAX_RADIUS = M_PI / sqrt(3.0);
  const int NB_SAMPLE = 25;
  const double STEP_RADIUS = (MAX_RADIUS - MIN_RADIUS) / (double)NB_SAMPLE;
  double max_error = 0.0;
  for (double r = MIN_RADIUS; r < MAX_RADIUS; r += STEP_RADIUS) {
    for (double theta = MIN_ANGLE; theta < MAX_ANGLE; theta += STEP_ANGLE) {
      for (double azimuth = MIN_ANGLE; azimuth < MAX_ANGLE;
           azimuth += STEP_ANGLE) {
        AlgebraSO3<double> lie(r * sin(theta) * cos(azimuth),
                               r * sin(theta) * sin(azimuth), r * cos(theta));
        AlgebraSO3<double> ln_exp_lie = lie.exp().log();

        Eigen::Vector3d error_v = ln_exp_lie.as_vector() - lie.as_vector();

        double error = error_v.transpose() * error_v;
        max_error = error > max_error ? error : max_error;

        ASSERT_NEAR(0.0, error, 1e-6);
      }
    }
  }
  std::cout << "Max error: " << max_error << std::endl;
}
