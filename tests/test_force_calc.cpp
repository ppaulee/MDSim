#include "gtest/gtest.h"

// Testing Lennard Jones Potential against hand-calculated values
TEST(ForceCalculation, LennardJonesPotential) {
    ForceCalculation *algorithm = new LennardJones(5, 1, 0);
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

// Testing Harmonic Potential against hand-calculated values
TEST(ForceCalculation, HarmonicPotential) {
    HarmonicPotential *algorithm = new HarmonicPotential(1, 1);
    auto *p1 = new Particle({0,0,0},{1,0,0},1,0);
    auto *p2 = new Particle({0,2,0},{1,1,0},1,0);
    std::array<double, 3> vec = algorithm->calculateF(*p1, *p2);
    EXPECT_EQ(0, vec[0]);
    EXPECT_DOUBLE_EQ(1, vec[1]);
    EXPECT_EQ(0, vec[2]);

    auto *p3 = new Particle({0,0,0},{1,0,0},1,0);
    auto *p4 = new Particle({0,5,0},{1,1,0},1,0);
    std::array<double, 3> vec2 = algorithm->calculateFDiag(*p3, *p4);
    EXPECT_EQ(0, vec2[0]);
    EXPECT_DOUBLE_EQ(3.58578643762690495119831127579030192143032812462305192682332026, vec2[1]);
    EXPECT_EQ(0, vec2[2]);
}