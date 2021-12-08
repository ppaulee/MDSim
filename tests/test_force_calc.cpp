#include "gtest/gtest.h"

// Testing Lennard Jones Potential against hand-calculated values
TEST(ForceCalculation, LennardJonesPotential) {
    ForceCalculation *algorithm = new LennardJones(5, 1);
    auto *p1 = new Particle({0,0,0},{1,0,0},1,0);
    auto *p2 = new Particle({0,2,0},{1,1,0},1,0);
    std::array<double, 3> vec = algorithm->calculateF(*p1, *p2);
    EXPECT_EQ(0, vec[0]);
    EXPECT_EQ(0.908203125, vec[1]);
    EXPECT_EQ(0, vec[2]);
}

// Testing Gravitation against hand-calculated values
TEST(ForceCalculation, Gravitation) {
    ForceCalculation *algorithm = new Gravitation();
    auto *p1 = new Particle({0,0,0},{1,0,0},1,0);
    auto *p2 = new Particle({0,2,0},{1,1,0},1,0);
    std::array<double, 3> vec = algorithm->calculateF(*p1, *p2);
    EXPECT_EQ(0, vec[0]);
    EXPECT_EQ(0.25, vec[1]);
    EXPECT_EQ(0, vec[2]);
}