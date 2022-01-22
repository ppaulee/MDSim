#include <gtest/gtest.h>

#include "container/ParticleContainer.cpp"
#include "../src/utils/MaxwellBoltzmannDistribution.cpp"
#include "../src/Particle.cpp"
#include "forceCalculation/LennardJones.cpp"
#include "forceCalculation/MixedLennardJones.cpp"
#include "forceCalculation/Gravitation.cpp"
#include "forceCalculation/HarmonicPotential.cpp"
#include "container/LinkedCells.cpp"
#include "../src/outputWriter/vtk-unstructured.cpp"
#include "../src/outputWriter/VTKWriter.cpp"

#include "test_force_calc.cpp"
#include "test_linkedcells.cpp"
#include "test_particle_container.cpp"
#include "test_particle_generator.cpp"
#include "test_thermostats.cpp"


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



