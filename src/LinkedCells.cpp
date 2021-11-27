//
// Created by paulu on 21.11.2021.
//

#include "LinkedCells.h"
#include "utils/ArrayUtils.h"
#include "LennardJones.h"
#include <iostream>
#include <algorithm>
#include <outputWriter/VTKWriter.h>
#include "utils/MaxwellBoltzmannDistribution.h"


// TODO 2D Grid
LinkedCells::LinkedCells(std::array<int, 3> dimension, double mesh, double cutOff) {
    cutOffRadius = cutOff;
    // Add boundary and halo cells
    if (dimension[2] == 0) {
        // 2D case
        dimensions = {dimension[0] + 4, dimension[1] + 4, 0};
        size = coordToIndex({dimension[0] + 4, dimension[1] + 4, 0});
    } else {
        // 3D Case
        dimensions = {dimension[0] + 4, dimension[1] + 4, dimension[2] + 4};
        size = coordToIndex({dimension[0] + 4, dimension[1] + 4, dimension[2] + 4});
    }

    meshSize = mesh;
    // Initialise all vectors
    for(int i = 0; i < coordToIndex({dimensions[0], dimensions[1], dimensions[2]})+1; i++) {
        particles.push_back(std::vector<Particle>({}));
    }
}

int LinkedCells::coordToIndex(std::array<int, 3> coords) {
    if (dimensions[2] == 0) {
        return coords[0] * dimensions[0] + coords[1];
    } else {
        return coords[0] * dimensions[1] * dimensions[2] + coords[1] * dimensions[2] + coords[2];
    }

}

std::array<int, 3> LinkedCells::indexToCoords(int index) {
    int width_index = index / (dimensions[1] * dimensions[2]);
    int height_index = (index - width_index * dimensions[1] * dimensions[2]) / dimensions[2];
    int depth_index = index - width_index * dimensions[1] * dimensions[2] - height_index * dimensions[2];
    return {width_index, height_index, depth_index};
}

void LinkedCells::test() {
    Particle* p = new Particle({2.2,1.2,1.2},{1,0,0},1,0);
    insert(*p);
    Particle* p1 = new Particle({3.3,1.2,1.2},{1,0,0},1,0);
    insert(*p1);
    std::cout << "### TEST: " << coordToIndex(getCellCoords(*p)) << std::endl;
    std::cout << "### TEST: " << coordToIndex(getCellCoords(*p1)) << std::endl;;
    std::cout << "### TEST: " << getCellCoords(*p)[0] << getCellCoords(*p)[1] << getCellCoords(*p)[2] << std::endl;
    calculateF(new LennardJones(5,1));
    std::cout << "TEST";
}

void LinkedCells::insert(Particle &p) {
    particles[coordToIndex(getCellCoords(p))].push_back(p);
}

std::array<int, 3> LinkedCells::getCellCoords(Particle &p) {
    int x = ((int) (p.getX()[0] / meshSize));
    int y = ((int) (p.getX()[1] / meshSize) );
    int z = ((int) (p.getX()[2] / meshSize) );
    return {x,y,z};
}

void LinkedCells::remove(Particle& p) {
    int currentIndex = 0;
    // Search for right element
    for (auto p1 : particles[coordToIndex(getCellCoords(p))]) {
        if (p1 == p) {
            break;
        }
        currentIndex++;
    }
    // Remove particle by index
    if (currentIndex == 0) {
        particles[coordToIndex(getCellCoords(p))] = {};
    } else {
        auto s = particles[coordToIndex(getCellCoords(p))];
        s.erase(s.begin() + currentIndex);
    }

}

bool LinkedCells::isHaloCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0]-1 || coords[1] == dimensions[1]-1 || coords[2] == dimensions[2]-1 || coords[0] == 0 || coords[1] == 0 || coords[2] == 0;
}

bool LinkedCells::isBoundaryCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0]-2 || coords[1] == dimensions[1]-2 || coords[2] == dimensions[2]-2 || coords[0] == 1 || coords[1] == 1 || coords[2] == 1;
}

std::vector<Particle>& LinkedCells::get(std::array<double, 3> coords) {
    int x = ((int) (coords[0] / meshSize));
    int y = ((int) (coords[1] / meshSize) );
    int z = ((int) (coords[2] / meshSize) );
    return particles[coordToIndex({x,y,z})];
}

void LinkedCells::calculateF(ForceCalculation *algorithm) {
    for (auto &vec : particles) {
        for (auto p : vec) {
            p.setOldF(p.getF());
            p.setF({0,0,0});
        }
    }

    for (int x = 1; x < dimensions[0]-1; ++x) {
        for (int y = 1; y < dimensions[1]-1; ++y) {
            if (dimensions[2] == 0) {
                //2D case
                if (particles[coordToIndex({x,y,0})].empty())
                    continue;

                // Loop over all particles in cell
                for (auto &current_particle : particles[coordToIndex({x,y,0})]) {
                    current_particle.mark();
                    // Get neighbours
                    for (int x_diff = -1; x_diff <= 1; ++x_diff) {
                        for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                            // Coordinates including offset for neighboured cells
                            std::array<int, 3> c = {x + x_diff, y + y_diff, 0};
                            if (particles[coordToIndex(c)].empty() || isHaloCell(c) || (x_diff == y_diff) == 0)
                                continue;

                            // Particles in neighboured cells
                            for (auto &p : particles[coordToIndex(c)]) {
                                // Check if particle is inside of the cut off radius
                                if (ArrayUtils::L2Norm(p.getX()-current_particle.getX()) <= cutOffRadius && !p.isMarked()) {
                                    // Actual force calculation according to the used algorithm
                                    std::array<double, 3> force = algorithm->calculateF(current_particle, p);
                                    // Make use of Newtons third law
                                    current_particle.setF(current_particle.getF() + force);
                                    p.setF(p.getF() - force);
                                }

                            }
                        }
                    }

                }
            }
            else {
                for (int z = 1; z < dimensions[2]-1; ++z) {
                    if (particles[coordToIndex({x,y,z})].empty())
                        continue;

                    // Loop over all particles in cell
                    for (auto &current_particle : particles[coordToIndex({x,y,z})]) {
                        current_particle.mark();
                        // Get neighbours
                        for (int x_diff = -1; x_diff <= 1; ++x_diff) {
                            for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                                for (int z_diff = -1; z_diff <= 1; ++z_diff) {

                                    // Coordinates including offset for neighboured cells
                                    std::array<int, 3> c = {x + x_diff, y + y_diff, z + z_diff};
                                    if (particles[coordToIndex(c)].empty() || isHaloCell(c) || ((x_diff == y_diff) == z_diff) == 0)
                                        continue;

                                    // Particles in neighboured cells
                                    for (auto &p : particles[coordToIndex(c)]) {
                                        // Check if particle is inside of the cut off radius
                                        if (ArrayUtils::L2Norm(p.getX()-current_particle.getX()) <= cutOffRadius && !p.isMarked()) {
                                            // Actual force calculation according to the used algorithm
                                            std::array<double, 3> force = algorithm->calculateF(current_particle, p);
                                            // Make use of Newtons third law
                                            current_particle.setF(current_particle.getF() + force);
                                            p.setF(p.getF() - force);
                                        }

                                    }
                                }
                            }
                        }

                    }

                }
            }

        }
    }
}

void LinkedCells::calculateV(double delta_t) {
    for (auto &vec: particles) {
        for (auto &p : vec) {
            std::array<double, 3> res = p.getF() + p.getOldF();
            res = (1 / (2 * p.getM())) * res;
            res = delta_t * res;
            res = res + p.getV();
            p.setV(res);
        }
    }
}

void LinkedCells::calculateX(double delta_t) {
    for (auto &vec: particles) {
        for (auto &p : vec) {
            std::array<double, 3> res = (1 / (2 * p.getM())) * p.getOldF();
            res = (delta_t * delta_t) * res;

            std::array<double, 3> res2 = delta_t * p.getV();
            res = res + res2;
            p.setX(res + p.getX());
        }
    }

}

void LinkedCells::move() {
    for (int i = 0; i < size; ++i) {
        for (auto &p : particles[i]) {
            if (coordToIndex(getCellCoords(p)) != i) {
                remove(p);
                insert(p);
            }
        }
    }
}

void LinkedCells::simulate(double delta_t, ForceCalculation *algorithm) {
    // calculate new x
    calculateX(delta_t);
    // calculate new f
    calculateF(algorithm);
    // calculate new v
    calculateV(delta_t);
    // moves particles to correct cell
    move();
}

void LinkedCells::plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(numberParticles());
    for (auto &vec : particles) {
        for (auto &p : vec) {
            writer.plotParticle(p);
        }

    }
    writer.writeFile(out_name, iteration);
}

void LinkedCells::addBrownianMotion(double averageV, int dimension) {
    for (auto &vec : particles) {
        for (auto &p : vec) {
            p.setV(p.getV() + maxwellBoltzmannDistributedVelocity(averageV, dimension));
        }
    }
}

int LinkedCells::numberParticles() {
    int result = 0;
    for (auto &vec : particles) {
        result += vec.size();
    }
    return result;
}


