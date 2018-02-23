
#include <iostream>

#include <gtest/gtest.h>

#include "util.hpp"

using namespace lieroy;

TEST(UtilTest, SampleFromSphere) {
    auto sample = sample_from_sphere<float>(5.0);
    std::cout << sample << std::endl;
}
