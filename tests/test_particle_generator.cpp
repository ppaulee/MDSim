//
// Created by paulu on 05.12.2021.
//

#include "gtest/gtest.h"
#include "../src/ParticleGenerator.cpp"
#include <fstream>
#include <string>

// Tests the generateFromFile function with given test file
TEST(ParticleGenerator, generate) {
    auto container = new ParticleContainer();

    generateFromFile(*container, (char*) "../../../../tests/test_input_particle_generator.txt");

    // There must be exactly one particle in the container
    EXPECT_EQ(container->getVec().size(), 1);

    // The position of the only particle must be (0,0,0)
    EXPECT_EQ(container->getVec()[0].getX()[0], 0);
    EXPECT_EQ(container->getVec()[0].getX()[1], 0);
    EXPECT_EQ(container->getVec()[0].getX()[2], 0);

    // The velocity must be (1,2,3)
    EXPECT_EQ(container->getVec()[0].getV()[0], 1);
    EXPECT_EQ(container->getVec()[0].getV()[1], 2);
    EXPECT_EQ(container->getVec()[0].getV()[2], 3);
}

