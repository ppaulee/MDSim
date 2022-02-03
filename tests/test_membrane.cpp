//
// Created by paulu on 23.01.2022.
//
#include "gtest/gtest.h"

bool contains(std::vector<std::shared_ptr<Particle>> list, Particle &p) {
    for (auto &pa : list) {
        double x = pa->getX()[0];
        double x_ = p.getX()[0];
        double y = pa->getX()[1];
        double y_ = p.getX()[1];
        double z = pa->getX()[2];
        double z_ = p.getX()[2];
        if (x == x_ && y == y_ && z == z_) {
            return true;
        }
    }
    return false;
}

TEST(Membrane, numberParticles) {
    auto cells = new LinkedCells({144, 144, 144}, 4, 4);

    std::vector<std::array<int, 3>> pull;
    std::array<double, 3> center = {15,15,50};
    std::array<double, 3> v = {0,0,0};
    std::array<int, 3> dim_ = {5,5,1};
    generateMembrane(center, v, dim_, 5, pull, 1,*cells);
    EXPECT_EQ(cells->getParticles().size(), 25);
}

TEST(Membrane, neighbours) {
    auto cells = new LinkedCells({144, 144, 144}, 4, 4);

    std::vector<std::array<int, 3>> pull;
    std::array<double, 3> center = {15,15,50};
    std::array<double, 3> v = {0,0,0};
    std::array<int, 3> dim_ = {5,5,1};
    generateMembrane(center, v, dim_, 5, pull, 1,*cells);

    for (auto &p : cells->getParticles()) {
        auto coords = p.getX();

        // top row
        int num = 0;
        int numDiag = 0;

        // bottom left corner
        if (coords[1] == 15 || coords[0] == 15) {
            numDiag++;
        }
        // bottom right corner
        if (coords[1] == 15 || coords[0] == 35) {
            numDiag++;
        }
        // top left corner
        if (coords[1] == 35 || coords[0] == 15) {
            numDiag++;
        }
        // top right corner
        if (coords[1] == 35 || coords[0] == 35) {
            numDiag++;
        }
        // bottom row
        if (coords[1] == 15) {
            num++;
        }
        // top row
        if (coords[1] == 35) {
            num++;
        }
        // left row
        if (coords[0] == 15) {
            num++;
        }
        // right row
        if (coords[0] == 35) {
            num++;
        }

        EXPECT_EQ(p.getNeighbours().size(), 4-num);
        EXPECT_EQ(p.getNeighboursDiag().size(), 4-numDiag);
    }
}
