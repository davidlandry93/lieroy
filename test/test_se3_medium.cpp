
#include <cmath>

#include <gtest/gtest.h>

#include "se3.hpp"

using namespace pylie;

void test_cycle_consistency(double min_angle, double max_angle, int nb_sample) {
    const double step_angle = (max_angle - min_angle) / (double)nb_sample;

    const double MIN_RADIUS = 0;
    const double MAX_RADIUS = M_PI/sqrt(3.0);
    const int NB_SAMPLE_RADIUS = 25;
    const double STEP_RADIUS =  (MAX_RADIUS - MIN_RADIUS) / (double)NB_SAMPLE_RADIUS;
    double max_error = 0.0;
    for (double r = MIN_RADIUS; r < MAX_RADIUS; r += STEP_RADIUS) {
        for (double theta = min_angle; theta < max_angle; theta += step_angle) {
            for (double azimuth = min_angle; azimuth < max_angle; azimuth += step_angle) {
                Eigen::Matrix<double,6,1> v;
                v << 0,
                    0,
                    0,
                    r * sin(theta) * cos(azimuth),
                    r * sin(theta) * sin(azimuth),
                    r * cos(theta);
                AlgebraSE3<double> lie(v);
                AlgebraSE3<double> ln_exp_lie = lie.exp().log();

                Eigen::Matrix<double,6,1> error_v = ln_exp_lie.as_vector() - lie.as_vector();

                double error = error_v.transpose()* error_v;
                max_error = error > max_error ? error : max_error;

                ASSERT_NEAR(0.0, error, 1e-6);
            }
        }
    }
    std::cout << "Max error: " << max_error << std::endl;
}

TEST(SE3, CycleConsistency) {
    test_cycle_consistency(-M_PI, M_PI, 25);
}

TEST(SE3, CycleConsistencySmallAngle) {
    test_cycle_consistency(-1e-6, 1e-6, 25);
}
