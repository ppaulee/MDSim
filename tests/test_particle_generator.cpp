//
// Created by paulu on 05.12.2021.
//

#include "gtest/gtest.h"
#include "../src/ParticleGenerator.cpp"
#include <fstream>
#include <string>
#include "utils/ArrayUtils.h"

// Tests the generateFromFile function with given test file
TEST(ParticleGenerator, generate) {
    auto container = new ParticleContainer();

    generateFromFileTest(*container, (std::string) "1\nCube:\n0,0,0;1,2,3;1,1,1;1.1225;0.1;1;");

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

TEST(ParticleGenerator, generateSphere2D){
    auto container = new ParticleContainer();

    generateSphere2D({0,0,0}, {0,0,0}, 2, 1, 1, 0.1, *container);
    int s= container->getVec().size();
    EXPECT_EQ(container->getVec().size(), 11);
    bool center = false;
    std::array<double, 3> comp = {0,0,0};
    for(auto &p : container->getVec()){
        if(p.getX() == comp)
            center = true;
    }
    EXPECT_TRUE(center);

}
