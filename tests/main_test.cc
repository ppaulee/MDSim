#include <gtest/gtest.h>

//#include "../src/utils/ArrayUtils.h"
#include "../src/Particle.h"
#include "../src/LennardJones.h"
#include "../src/Particle.cpp"
#include "../src/LennardJones.cpp"


// Testing Lennard Jones Potential against by hand calculated values
TEST(ForceCalculation, LennardJonesPotential) {
    ForceCalculation *algorithm = new LennardJones(5, 1);
    auto *p1 = new Particle({0,0,0},{1,0,0},1,0);
    auto *p2 = new Particle({0,2,0},{1,1,0},1,0);
    std::array<double, 3> vec = algorithm->calculateF(*p1, *p2);
    ASSERT_EQ(-0.908203125, vec[1]);
}