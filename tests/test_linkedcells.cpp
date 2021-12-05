#include "gtest/gtest.h"

// Tests correct construction of cell grid
TEST(LinkedCells, gridConstruction) {
    auto cells = new LinkedCells({2,2,2}, 1, 1);
    EXPECT_EQ(1, 1);
    for (int x = 0; x < 6; ++x) {
        for (int y = 0; y < 6; ++y) {
            for (int z = 0; z < 6; ++z) {
                bool halo = (x == 0 || y == 0 || z == 0 || x == 5 || y == 5 || z == 5);
                bool boundary = (x == 1 || y == 1 || z == 1 || x == 4 || y == 4 || z == 4);

                EXPECT_EQ(halo, cells->isHaloCell({x,y,z}));
                EXPECT_EQ(boundary, cells->isBoundaryCell({x,y,z}));
            }
        }
    }
}

// Create 2 particles, then the function numberParticles() must return 2
TEST(LinkedCells, numberParticles) {
    auto cells = new LinkedCells({2,2,2}, 1, 1);
    auto p1 = new Particle({0,0,0},{0,0,0},1,1);
    auto p2 = new Particle({1,1,1},{0,0,0},1,1);
    cells->insert(*p1);
    cells->insert(*p2);

    EXPECT_EQ(cells->numberParticles(), 2);
}

// Create particle then change the position and then let the datastrucutre move the particle into the right cell
TEST(LinkedCells, move) {
    auto cells = new LinkedCells({10,10,10}, 1, 1);
    auto p1 = new Particle({0,0,0},{0,0,0},1,1);
    cells->insert(*p1);
    p1->setX({0,0,0});

    cells->move();

    EXPECT_TRUE(cells->allInRightCell());
}