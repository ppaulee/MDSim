//
// Created by paulu on 21.11.2021.
//
#include "LinkedCells.h"
#include "utils/ArrayUtils.h"
#include <outputWriter/VTKWriter.h>
#include <stdexcept>
#include <iostream>


LinkedCells::LinkedCells(std::array<int, 3> dimension, double mesh, double cutOff, double sigma) {
    if (dimension[0] % 2 != 0 || dimension[1] % 2 != 0 || dimension[2] % 2 != 0) {
        throw std::invalid_argument("Dimensions must be even");
    }

    reflectionDistance = cbrt(sqrt(2)) * sigma;
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
    for (int i = 0; i < coordToIndex({dimensions[0], dimensions[1], dimensions[2]}) + 1; i++) {
        particles.push_back(std::list<Particle>({}));
    }
}

int LinkedCells::coordToIndex(std::array<int, 3> coords) {
    if (dimensions[2] == 0) {
        return (coords[0]) * dimensions[0] + (coords[1]);
    } else {
        return coords[0] * dimensions[1] * dimensions[2] + coords[1] * dimensions[2] + coords[2];
    }
}

// Only 3D
std::array<int, 3> LinkedCells::indexToCoords(int index) {
    int width_index = index / (dimensions[1] * dimensions[2]);
    int height_index = (index - width_index * dimensions[1] * dimensions[2]) / dimensions[2];
    int depth_index = index - width_index * dimensions[1] * dimensions[2] - height_index * dimensions[2];
    return {width_index, height_index, depth_index};
}

void LinkedCells::insert(Particle &p) {
    // Transfer the coordinates into a system where (0,0,0) is the smallest coordinate
    if (dimensions[2] == 0) {
        std::array<double, 3> tmp = {(double) dimensions[0] / 2, (double) dimensions[1] / 2, 0};
        p.setX(tmp + p.getX());
    } else {
        std::array<double, 3> tmp = {(double) dimensions[0] / 2, (double) dimensions[1] / 2,
                                     (double) dimensions[2] / 2};
        p.setX(tmp + p.getX());
    }
    particles[coordToIndex(getCellCoords(p))].push_back(p);
}

void LinkedCells::forceInsert(Particle &p) {

    particles[coordToIndex(getCellCoords(p))].push_back(p);
}

std::array<int, 3> LinkedCells::getCellCoords(Particle &p) {
    int x = ((int) (p.getX()[0] / meshSize));
    int y = ((int) (p.getX()[1] / meshSize));
    int z = ((int) (p.getX()[2] / meshSize));
    return {x, y, z};
}

void LinkedCells::remove(Particle &p, int index = -1) {
    int i = index;
    if (index == -1) {
        i = coordToIndex(getCellCoords(p));
    }
    particles[i].remove(p);
}

bool LinkedCells::isHaloCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0] - 1 || coords[1] == dimensions[1] - 1 || coords[2] == dimensions[2] - 1 ||
           coords[0] == 0 || coords[1] == 0 || (coords[2] == 0 && dimensions[2] != 0);
}

bool LinkedCells::isBoundaryCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0] - 2 || coords[1] == dimensions[1] - 2 || coords[2] == dimensions[2] - 2 ||
           coords[0] == 1 || coords[1] == 1 || coords[2] == 1;
}

std::list<Particle> &LinkedCells::get(std::array<double, 3> coords) {
    int x = ((int) (coords[0] / meshSize));
    int y = ((int) (coords[1] / meshSize));
    int z = ((int) (coords[2] / meshSize));
    return particles[coordToIndex({x, y, z})];
}

std::list<Particle> &LinkedCells::indexGet(int i) {
    return particles[i];
}

void LinkedCells::calculateF(ForceCalculation *algorithm) {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setOldF(p.getF());
            p.setF({0, 0, 0});
            p.unmark();
        }
    }

    for (int x = 0; x < dimensions[0] - 1; ++x) {
        for (int y = 0; y < dimensions[1] - 1; ++y) {
            if (dimensions[2] == 0) {
                //2D case
                // Loop over all particles in cell
                for (auto &current_particle: particles[coordToIndex({x, y, 0})]) {
                    current_particle.mark();
                    // Get neighbours
                    for (int x_diff = -1; x_diff <= 1; ++x_diff) {
                        for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                            // Coordinates including offset for neighboured cells
                            std::array<int, 3> c = {x + x_diff, y + y_diff, 0};

                            if (coordToIndex(c) < 0) {
                                continue;
                            }

                            // Particles in neighboured cells
                            for (auto &p: particles[coordToIndex(c)]) {
                                // Check if particle is inside the cut-off radius
                                if (ArrayUtils::L2Norm(p.getX() - current_particle.getX()) <= cutOffRadius &&
                                    !(p == current_particle) && !p.isMarked()) {
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
            } else {
                for (int z = 1; z < dimensions[2] - 1; ++z) {
                    // Loop over all particles in cell
                    for (auto &current_particle: particles[coordToIndex({x, y, z})]) {
                        current_particle.mark();
                        // Get neighbours
                        for (int x_diff = -1; x_diff <= 1; ++x_diff) {
                            for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                                for (int z_diff = -1; z_diff <= 1; ++z_diff) {

                                    // Coordinates including offset for neighboured cells
                                    std::array<int, 3> c = {x + x_diff, y + y_diff, z + z_diff};
                                    // Particles in neighboured cells
                                    for (auto &p: particles[coordToIndex(c)]) {
                                        // Check if particle is inside of the cut off radius
                                        if (ArrayUtils::L2Norm(p.getX() - current_particle.getX()) <= cutOffRadius &&
                                            !(p == current_particle) && !p.isMarked()) {
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
        for (auto &p: vec) {
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
        for (auto &p: vec) {
            std::array<double, 3> res = (1 / (2 * p.getM())) * p.getOldF();
            res = (delta_t * delta_t) * res;

            std::array<double, 3> res2 = delta_t * p.getV();
            res = res + res2;

            p.setX(res + p.getX());
        }
    }
}

void LinkedCells::move() {
    std::vector<Particle> tmp_stack_insert = {};
    std::vector<Particle> tmp_stack_remove = {};
    std::vector<int> tmp_stack_remove_index = {};
    // Find particles in wrong cells and add them to the vectors above
    for (int i = 0; i < size; ++i) {
        for (auto &p: particles[i]) {
            if (coordToIndex(getCellCoords(p)) != i) {
                tmp_stack_remove.push_back(p);
                tmp_stack_insert.push_back(p);
                tmp_stack_remove_index.push_back(i);
            }
            p.unmark();
        }
    }
    // Remove and reinsert particles in wrong cells
    for (long unsigned int i = 0; i < tmp_stack_insert.size(); ++i) {
        // remove
        auto &to_remove = tmp_stack_remove[i];
        int index = tmp_stack_remove_index[i];
        particles[index].remove(to_remove);

        // insert
        auto &p = tmp_stack_insert[i];

        if (coordToIndex(getCellCoords(p)) >= 0 &&
            (long unsigned int) coordToIndex(getCellCoords(p)) < particles.size()) {
            particles[coordToIndex(getCellCoords(p))].push_back(p);
        }

    }
}

bool LinkedCells::allInRightCell() {
    for (int i = 0; i < size; ++i) {
        for (auto &p: particles[i]) {
            if (coordToIndex(getCellCoords(p)) != i) {
                return false;
            }
        }
    }
    return true;
}

void LinkedCells::simulate(ForceCalculation *algorithm, double delta_t) {
    // calculate new x
    calculateX(delta_t);

    //deleteParticlesInHalo();
    // moves particles to correct cell
    move();
    // TODO here boundary conditions
    //reflection with ghost particle
    createGhosts();

    // calculate new f
    calculateF(algorithm);
    deleteParticlesInHalo();
    // calculate new v
    calculateV(delta_t);
}

void LinkedCells::createGhosts() {
    std::list<Particle> ghosts;
    for (auto &vec: particles) {
        if (!vec.empty()) {
            //if (isBoundaryCell(getCellCoords(vec.front()))) {
            for (auto &p: vec) {
                // check distance of the particle to every boundary, create ghost particle if it's smaller than reflectionDistance
                // add all ghost particles to be added to a list and create them after iterating over the particles
                double leftDistance = p.getX()[0] - meshSize;
                if (leftDistance <= reflectionDistance) {
                    ghosts.emplace_front(
                            Particle({(meshSize - leftDistance), p.getX()[1], p.getX()[2]}, {0, 0, 0}, p.getM(),
                                     p.getType()));
                    std::cout << "left ghost created:\n";
                }
                double bottomDistance = p.getX()[1] - meshSize;
                if (bottomDistance <= reflectionDistance) {
                    ghosts.emplace_front(
                            Particle({p.getX()[0], (meshSize - bottomDistance), p.getX()[2]}, {0, 0, 0}, p.getM(),
                                     p.getType()));
                    std::cout << "bottom ghost created:\n";
                    std::cout << p.getX()[0] << "," << (meshSize - bottomDistance) << "," << p.getX()[2] << "\n";
                    std::cout << "Original:\n";
                    std::cout << p.getX()[0] << "," << p.getX()[1] << "," << p.getX()[2] << "\n";
                    //exit(1);
                }
                double rightDistance = (dimensions[0] - 1) * meshSize - p.getX()[0];
                if (rightDistance <= reflectionDistance) {
                    ghosts.emplace_front(
                            Particle({((dimensions[0] - 1) * meshSize + rightDistance), p.getX()[1], p.getX()[2]},
                                     {0, 0, 0},
                                     p.getM(),
                                     p.getType()));
                    std::cout << "right ghost created:\n";
                }
                double topDistance = (dimensions[1] - 1) * meshSize - p.getX()[1];
                if (topDistance <= reflectionDistance) {
                    ghosts.emplace_front(
                            Particle({(p.getX()[0]), (dimensions[1] - 1) * meshSize + topDistance, p.getX()[2]},
                                     {0, 0, 0}, p.getM(),
                                     p.getType()));
                    std::cout << "top ghost created:\n";
                }

                //}
            }

        }
    }
    for (auto par: ghosts) {
        forceInsert(par);
    }
    //ghosts.clear();
}

void LinkedCells::plotParticles(int iteration) {
    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(numberParticles());
    for (auto &vec: particles) {
        for (auto &p: vec) {
            std::array<double, 3> tmp;
            if (dimensions[2] == 0) {
                // Transfer the coordinates back to original scheme
                tmp = {(double) dimensions[0] / 2, (double) dimensions[1] / 2, 0};
                p.setX((std::array<double, 3>) (p.getX() - tmp));
            } else {
                // Transfer the coordinates back to original scheme
                tmp = {(double) dimensions[0] / 2, (double) dimensions[1] / 2, (double) dimensions[2] / 2};
                p.setX((std::array<double, 3>) (p.getX() - tmp));
            }
            writer.plotParticle(p);
            // Transfer the coordinates back to calculate them correctly again
            p.setX((std::array<double, 3>) (p.getX() + tmp));
        }

    }
    writer.writeFile(out_name, iteration);
}

void LinkedCells::addBrownianMotion(double averageV, int dimension) {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setV(p.getV() + maxwellBoltzmannDistributedVelocity(averageV, dimension));
        }
    }
}

int LinkedCells::numberParticles() {
    int result = 0;
    for (auto &vec: particles) {
        result += vec.size();
    }
    return result;
}

void LinkedCells::deleteParticlesInHalo() {
    for (int x = 0; x < dimensions[0] - 1; x++) {
        for (int y = 0; y < dimensions[1] - 1; y++) {
            if (dimensions[2] == 0) {
                if (isHaloCell({x, y, 0}) && !particles.at(coordToIndex({x, y, 0})).empty()) {
                    std::list<Particle> newList;
                    particles.at(coordToIndex({x, y, 0})) = newList;
                    std::cout << "Particle in halo detected \n";
                    //exit(1);
                }
            } else {
                for (int z = 1; z < dimensions[2] - 1; ++z) {
                    if (isHaloCell({x, y, z})) {
                        particles[coordToIndex({x, y, z})].clear();
                    }
                }
            }
        }
    }
}

std::array<int, 3> LinkedCells::getDimensions() {
    return dimensions;
}




