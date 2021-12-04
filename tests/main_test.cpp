#include <gtest/gtest.h>

#include "../src/ParticleContainer.cpp"
#include "../src/Particle.cpp"
#include "../src/LennardJones.cpp"
#include "../src/StoermerVerlet.cpp"
#include "../src/LinkedCells.cpp"
#include "../src/outputWriter/vtk-unstructured.cpp"
#include "../src/outputWriter/VTKWriter.cpp"

#include "test_force_calc.cpp"
#include "test_linkedcells.cpp"
#include "test_particle_container.cpp"


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



