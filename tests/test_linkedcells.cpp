#include "gtest/gtest.h"

// Tests correct construction of cell grid
TEST(LinkedCells, gridConstruction) {
    auto cells = new LinkedCells({2, 2, 2}, 1, 5);
    EXPECT_EQ(1, 1);
    for (int x = 0; x < 6; ++x) {
        for (int y = 0; y < 6; ++y) {
            for (int z = 0; z < 6; ++z) {
                bool halo = (x == 0 || y == 0 || z == 0 || x == 5 || y == 5 || z == 5);
                bool boundary = (x == 1 || y == 1 || z == 1 || x == 4 || y == 4 || z == 4);

                EXPECT_EQ(halo, cells->isHaloCell({x, y, z}));
                EXPECT_EQ(boundary, cells->isBoundaryCell({x, y, z}));
            }
        }
    }
}

// Create 2 particles, then the function numberParticles() must return 2
TEST(LinkedCells, numberParticles) {
    auto cells = new LinkedCells({2, 2, 2}, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 0, 0}, 1, 1);
    auto p2 = new Particle({1, 1, 1}, {0, 0, 0}, 1, 1);
    cells->insert(*p1);
    cells->insert(*p2);

    EXPECT_EQ(cells->numberParticles(), 2);
}

// Create particle then change the position and then let the datastrucutre move the particle into the right cell
TEST(LinkedCells, move) {
    auto cells = new LinkedCells({10, 10, 10}, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 0, 0}, 1, 1);
    cells->insert(*p1);
    p1->setX({0, 0, 0});

    cells->move();

    EXPECT_TRUE(cells->allInRightCell());
    EXPECT_TRUE(cells->numberParticles() == 1);
}

TEST(LinkedCells, createGhost) {
    std::array<int, 3> bounds = {1, 1, 1};
    auto cells = new LinkedCells({4, 4, 0}, 3, 3, 0, bounds);
    auto p1 = new Particle({3.5, 8, 0}, {0, 0, 0}, 1, 0, 1, 1);
    auto p2 = new Particle({2.5, 8, 0}, {0, 0, 0}, 1, 0, 1, 1);
    auto calc = LennardJones(1, 1, 0);
    cells->forceInsert(*p1);
    cells->handleBoundary();
    cells->applyForceBuffer();
    std::list<Particle> cell1 = cells->indexGet(cells->coordToIndex({1, 2, 0}));
    //cells->move();
    EXPECT_TRUE(cell1.size() == 1);
    EXPECT_TRUE(cell1.front().getF() == calc.calculateF(*p1, *p2));

}

TEST(LinkedCells, clear) {
    std::vector<std::list<int>> test;
    test.push_back(std::list<int>({1}));
    EXPECT_TRUE(test[0].size() == 1);
    EXPECT_TRUE(test.size() == 1);
    test.at(0).clear();
    EXPECT_TRUE(test[0].size() == 0);
    EXPECT_TRUE(test.size() == 1);

}

TEST(LinkedCells, applyGravity) {
    std::array<int, 3> bounds = {0, 0, 0};
    auto cells = new LinkedCells({4, 4, 0}, 3, 3, -10, bounds);
    auto p1 = new Particle({4, 8, 0}, {0, 0, 0}, 1, 1);
    cells->forceInsert(*p1);
    cells->applyGravity();
    std::list<Particle> cell1 = cells->indexGet(cells->coordToIndex({1, 2, 0}));
    std::array<double, 3> comp = {0, -10, 0};
    EXPECT_TRUE(cell1.size() == 1);
    EXPECT_TRUE(cell1.front().getF() == comp);
}

TEST(LinkedCells, periodicBoundaryGhost) {
    std::array<int, 3> bounds = {2, 0, 0};
    auto cells = LinkedCells({4, 4, 0}, 3, 3, 0, bounds);
    auto p1 = Particle({3.5, 8, 0}, {0, 0, 0}, 1, 0, 1, 5);
    auto p1Ghost = Particle({21.5, 8, 0}, {0, 0, 0}, 1, 0, 1, 5);
    auto p2 = Particle({20.5, 8, 0}, {0, 0, 0}, 1, 0, 1, 5);
    cells.forceInsert(p1);
    cells.forceInsert(p2);
    cells.handleBoundary();
    ForceCalculation *calc = new LennardJones(5, 1);
    cells.calculateF(calc);
    std::vector<Particle> pars = cells.getParticles();
    for (auto &p: pars) {
        if (p.getX()[0] > 10) {
            std::array<double, 3> comp = calc->calculateF(p2,p1Ghost);
            EXPECT_TRUE(p.getF() == comp);
        }
    }
}