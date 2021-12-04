#include "gtest/gtest.h"

// Testing size method of Particle Container
TEST(ParticleContainer, Size) {
    ParticleContainer c;
    std::array<double,3> x = {0,0,0};
    double m = 1;
    c.emplace_back(x,x,m);
    EXPECT_EQ(c.size(),1);
}

// Testing getVec Method of ParticleContainer
TEST(ParticleContainer, getVec) {
    ParticleContainer c;
    std::array<double,3> x = {0,0,0};
    double m = 1;
    c.emplace_back(x,x,m);
    std::vector<Particle> result;
    result.emplace_back(x,x,m);
    EXPECT_EQ(c.getVec()[0].getX(), result[0].getX());
    EXPECT_EQ(c.getVec()[0].getV(), result[0].getV());
    EXPECT_EQ(c.getVec()[0].getM(), result[0].getM());
}