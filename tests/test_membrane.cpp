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
        bool top_bool = true;
        bool bottom_bool = true;
        bool left_bool = true;
        bool right_bool = true;
        bool top_left_bool = true;
        bool top_right_bool = true;
        bool bottom_left_bool = true;
        bool bottom_right_bool = true;

        // bottom left corner
        if (coords[1] == 15 || coords[0] == 15) {
            bottom_left_bool = false;
        }
        // bottom right corner
        if (coords[1] == 15 || coords[0] == 35) {
            bottom_right_bool = false;
        }
        // top left corner
        if (coords[1] == 35 || coords[0] == 15) {
            top_left_bool = false;
        }
        // top right corner
        if (coords[1] == 35 || coords[0] == 35) {
            top_right_bool = false;
        }
        // bottom row
        if (coords[1] == 15) {
            bottom_bool = false;
        }
        // top row
        if (coords[1] == 35) {
            top_bool = false;
        }
        // left row
        if (coords[0] == 15) {
            left_bool = false;
        }
        // right row
        if (coords[0] == 35) {
            right_bool = false;
        }

        if (top_bool) {
            std::array<double, 3> top = coords;
            top[1] += 5;
            EXPECT_TRUE(contains(p.getNeighbours(), cells->get(top).front()));
        }

        if (bottom_bool) {
            std::array<double, 3> bottom = coords;
            bottom[1] -= 5;
            EXPECT_TRUE(contains(p.getNeighbours(), cells->get(bottom).front()));
        }

        if (left_bool) {
            std::array<double, 3> left = coords;
            left[0] -= 5;
            EXPECT_TRUE(contains(p.getNeighbours(), cells->get(left).front()));
        }

        if (right_bool) {
            std::array<double, 3> right = coords;
            right[0] += 5;
            EXPECT_TRUE(contains(p.getNeighbours(), cells->get(right).front()));
        }

        if (top_left_bool) {
            std::array<double, 3> top_left = coords;
            top_left[0] -= 5;
            top_left[1] += 5;
            EXPECT_TRUE(contains(p.getNeighboursDiag(), cells->get(top_left).front()));
        }

        if (top_right_bool) {
            std::array<double, 3> top_right = coords;
            top_right[0] += 5;
            top_right[1] += 5;
            EXPECT_TRUE(contains(p.getNeighboursDiag(), cells->get(top_right).front()));
        }

        if (bottom_left_bool) {
            std::array<double, 3> bottom_left = coords;
            bottom_left[0] -= 5;
            bottom_left[1] -= 5;
            EXPECT_TRUE(contains(p.getNeighboursDiag(), cells->get(bottom_left).front()));
        }

        if (bottom_right_bool) {
            std::array<double, 3> bottom_right = coords;
            bottom_right[0] += 5;
            bottom_right[1] -= 5;
            EXPECT_TRUE(contains(p.getNeighboursDiag(), cells->get(bottom_right).front()));
        }
    }
}

