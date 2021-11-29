#include <gtest/gtest.h>

#include "../src/ParticleContainer.h"
#include "../src/ParticleContainer.cpp"
#include "../src/LennardJones.h"
#include "../src/StoermerVerlet.h"
#include "../src/Particle.cpp"
#include "../src/LennardJones.cpp"
#include "../src/StoermerVerlet.cpp"
#include "../src/LinkedCells.cpp"
#include "../src/outputWriter/vtk-unstructured.h"
#include "../src/outputWriter/vtk-unstructured.cpp"
#include "../src/outputWriter/VTKWriter.h"
#include "../src/outputWriter/VTKWriter.cpp"





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

// Testing Stoermer Verlet against hand-calculated values
TEST(ForceCalculation, StoermerVerlet) {
    ForceCalculation *algorithm = new StoermerVerlet();
    auto *p1 = new Particle({0,0,0},{1,0,0},1,0);
    auto *p2 = new Particle({0,2,0},{1,1,0},1,0);
    std::array<double, 3> vec = algorithm->calculateF(*p1, *p2);
    EXPECT_EQ(0, vec[0]);
    EXPECT_EQ(0.25, vec[1]);
    EXPECT_EQ(0, vec[2]);
}

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

// Insert a particle and then remove it. Then check if the list is empty
TEST(LinkedCells, insertAndRemove) {
    auto cells = new LinkedCells({1,1,1}, 2, 1);
    auto p = new Particle({0,0,0},{1,0,0},1,0);
    auto p1 = new Particle({1,0,0},{1,0,0},1,0);
    cells->insert(*p);
    cells->insert(*p1);
    cells->remove(*p);
    EXPECT_EQ(cells->get({0,0,0}).size(), 1);
}

// Tests correct construction of cell grid
TEST(LinkedCells, gridConstruction) {
    auto cells = new LinkedCells({1,1,1}, 1, 1);

    EXPECT_EQ(cells->isHaloCell({0,0,0}), true);

    EXPECT_EQ(cells->isHaloCell({1,0,0}), true);
    EXPECT_EQ(cells->isHaloCell({2,0,0}), true);
    EXPECT_EQ(cells->isHaloCell({3,0,0}), true);
    EXPECT_EQ(cells->isHaloCell({4,0,0}), true);

    EXPECT_EQ(cells->isHaloCell({1,4,0}), true);
    EXPECT_EQ(cells->isHaloCell({2,4,0}), true);
    EXPECT_EQ(cells->isHaloCell({3,4,0}), true);
    EXPECT_EQ(cells->isHaloCell({4,4,0}), true);

    EXPECT_EQ(cells->isHaloCell({1,0,4}), true);
    EXPECT_EQ(cells->isHaloCell({2,0,4}), true);
    EXPECT_EQ(cells->isHaloCell({3,0,4}), true);
    EXPECT_EQ(cells->isHaloCell({4,0,4}), true);

    EXPECT_EQ(cells->isHaloCell({0,1,0}), true);
    EXPECT_EQ(cells->isHaloCell({0,2,0}), true);
    EXPECT_EQ(cells->isHaloCell({0,3,0}), true);
    EXPECT_EQ(cells->isHaloCell({0,4,0}), true);

    EXPECT_EQ(cells->isHaloCell({0,1,4}), true);
    EXPECT_EQ(cells->isHaloCell({0,2,4}), true);
    EXPECT_EQ(cells->isHaloCell({0,3,4}), true);
    EXPECT_EQ(cells->isHaloCell({0,4,4}), true);

    EXPECT_EQ(cells->isHaloCell({4,1,0}), true);
    EXPECT_EQ(cells->isHaloCell({4,2,0}), true);
    EXPECT_EQ(cells->isHaloCell({4,3,0}), true);
    EXPECT_EQ(cells->isHaloCell({4,4,0}), true);

    EXPECT_EQ(cells->isHaloCell({0,0,1}), true);
    EXPECT_EQ(cells->isHaloCell({0,0,2}), true);
    EXPECT_EQ(cells->isHaloCell({0,0,3}), true);
    EXPECT_EQ(cells->isHaloCell({0,0,4}), true);

    // Boundary cells
    EXPECT_EQ(cells->isBoundaryCell({1,1,1}), true);

    EXPECT_EQ(cells->isBoundaryCell({2,1,1}), true);
    EXPECT_EQ(cells->isBoundaryCell({3,1,1}), true);

    EXPECT_EQ(cells->isBoundaryCell({1,2,1}), true);
    EXPECT_EQ(cells->isBoundaryCell({1,3,1}), true);

    EXPECT_EQ(cells->isBoundaryCell({1,1,2}), true);
    EXPECT_EQ(cells->isBoundaryCell({1,1,3}), true);

    EXPECT_EQ(cells->isBoundaryCell({3,1,2}), true);
    EXPECT_EQ(cells->isBoundaryCell({3,1,3}), true);

    EXPECT_EQ(cells->isBoundaryCell({3,2,1}), true);
    EXPECT_EQ(cells->isBoundaryCell({3,3,1}), true);

    EXPECT_EQ(cells->isBoundaryCell({1,2,3}), true);
    EXPECT_EQ(cells->isBoundaryCell({1,3,3}), true);

    EXPECT_EQ(cells->isBoundaryCell({3,2,3}), true);
    EXPECT_EQ(cells->isBoundaryCell({3,3,3}), true);

    EXPECT_EQ(cells->isHaloCell({2,2,2}), false);
    EXPECT_EQ(cells->isBoundaryCell({2,2,2}), false);
}

