//
// Created by paulu on 21.11.2021.
//
#include "LinkedCells.h"
#include "utils/ArrayUtils.h"
#include <outputWriter/VTKWriter.h>
#include <stdexcept>
#include <iostream>


LinkedCells::LinkedCells(std::array<int, 3> dimension, double mesh, double cutOff, double gravConst,
                         std::array<int, 3> &boundaryConds) {
    if (dimension[0] % 2 != 0 || dimension[1] % 2 != 0 || dimension[2] % 2 != 0) {
        throw std::invalid_argument("Dimensions must be even");
    }
    grav = gravConst;
    boundaryCondition = boundaryConds;
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

LinkedCells::LinkedCells(std::array<int, 3> dimension, double mesh, double cutOff) {
    if (dimension[0] % 2 != 0 || dimension[1] % 2 != 0 || dimension[2] % 2 != 0) {
        throw std::invalid_argument("Dimensions must be even");
    }
    grav = 0;
    boundaryCondition = {0, 0, 0};
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
    if (forceCalcs.size() == p.getType()) {
        forceCalcs.push_back(LennardJones(p.getEpsilon(), p.getSigma(), p.getType()));
        reflectionDistance.push_back((cbrt(sqrt(2)) * p.getSigma()) / 2);
        if (forceCalcs.size() >= 1) {
            for (auto &calc: forceCalcs) {
                if (calc.getType() != p.getType()) {
                    double tempEpsilon = sqrt(calc.getEpsilon() * p.getEpsilon());
                    double tempSigma = (calc.getSigma() + p.getSigma()) / 2;
                    auto temp = MixedLennardJones(tempEpsilon, tempSigma, 0, p.getType());
                    mixedForceCalcs.push_back(temp);
                }
            }
        }
    }
    particles[coordToIndex(getCellCoords(p))].push_back(p);
}

void LinkedCells::forceInsert(Particle &p) {
    if (forceCalcs.size() == p.getType()) {
        forceCalcs.push_back(LennardJones(p.getEpsilon(), p.getSigma(), p.getType()));
        reflectionDistance.push_back((cbrt(sqrt(2)) * p.getSigma()) / 2);
        if (forceCalcs.size() >= 1) {
            for (auto &calc: forceCalcs) {
                if (calc.getType() != p.getType()) {
                    double tempEpsilon = sqrt(calc.getEpsilon() * p.getEpsilon());
                    double tempSigma = (calc.getSigma() + p.getSigma()) / 2;
                    auto temp = MixedLennardJones(tempEpsilon, tempSigma, 0, p.getType());
                    mixedForceCalcs.push_back(temp);
                }
            }
        }
    }
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
    initCalculateF();
    for (int x = 0; x < dimensions[0] - 1; ++x) {
        for (int y = 0; y < dimensions[1] - 1; ++y) {
            if (dimensions[2] == 0) {
                TwoDcalculateF(x, y, algorithm);
            } else {
                for (int z = 0; z < dimensions[2] - 1; ++z) {
                    ThreeDcalculateF(x, y, z, algorithm);
                }
            }

        }
    }
}

void LinkedCells::initCalculateF() {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setOldF(p.getF());
            p.setF({0, 0, 0});
            p.unmark();
        }
    }
}

void LinkedCells::TwoDcalculateF(int x, int y, ForceCalculation *algorithm) {
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
                calculateNeighbouredF(c, algorithm, current_particle);
            }
        }
    }
}

void LinkedCells::ThreeDcalculateF(int x, int y, int z, ForceCalculation *algorithm) {
    // Loop over all particles in cell
    for (auto &current_particle: particles[coordToIndex({x, y, z})]) {
        current_particle.mark();
        // Get neighbours
        for (int x_diff = -1; x_diff <= 1; ++x_diff) {
            for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                for (int z_diff = -1; z_diff <= 1; ++z_diff) {

                    // Coordinates including offset for neighboured cells
                    std::array<int, 3> c = {x + x_diff, y + y_diff, z + z_diff};

                    if (coordToIndex(c) < 0) {
                        continue;
                    }
                    calculateNeighbouredF(c, algorithm, current_particle);
                }
            }
        }
    }
}

void LinkedCells::calculateNeighbouredF(std::array<int, 3> c, ForceCalculation *algorithm, Particle& current_particle) {
    // Particles in neighboured cells
    for (auto &p: particles[coordToIndex(c)]) {
        // Check if particle is inside the cut-off radius
        if (ArrayUtils::L2Norm(p.getX() - current_particle.getX()) <= cutOffRadius && !(p == current_particle) && !p.isMarked()) {
            // Actual force calculation according to the used algorithm
            std::array<double, 3> force = algorithm->calculateF(current_particle, p);

            // Make use of Newtons third law
            current_particle.setF(current_particle.getF() + force);
            p.setF(p.getF() - force);
        }
    }
    ghosts.clear();
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
    handleBoundary();


    // calculate new f
    calculateF(algorithm);
    deleteParticlesInHalo();
    if (grav != 0) {
        applyGravity();
    }
    // calculate new v
    calculateV(delta_t);
}

void LinkedCells::handleBoundary() {
    std::list<Particle> ins;
    std::list<Particle> del;
    //int c = 0;
    for (auto &vec: particles) {
        if (!vec.empty()) {
            if (isBoundaryCell(getCellCoords(vec.front())) || isHaloCell(getCellCoords(vec.front()))) {
                for (auto &p: vec) {
                    //std::cout << "par nr " << c << "\n";
                    //c++;
                    if (boundaryCondition[0] == 1) {
                        if (isBoundaryCell(getCellCoords(p))) {
                            double leftDistance = p.getX()[0] - meshSize;
                            if (leftDistance <= reflectionDistance[p.getType()]) {
                                Particle temp = Particle({(meshSize - leftDistance), p.getX()[1], p.getX()[2]},
                                                         {0, 0, 0},
                                                         p.getM(),
                                                         p.getType(), p.getSigma(), p.getEpsilon());
                                p.setF(p.getF() + forceCalcs[p.getType()].calculateF(p, temp));

                            }
                            double rightDistance = (dimensions[0] - 1) * meshSize - p.getX()[0];
                            if (rightDistance <= reflectionDistance[p.getType()]) {
                                Particle temp = Particle(
                                        {((dimensions[0] - 1) * meshSize + rightDistance), p.getX()[1], p.getX()[2]},
                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                p.setF(p.getF() + forceCalcs[p.getType()].calculateF(p, temp));
                            }
                        }
                    } else if (boundaryCondition[0] == 2) {
                        // It's in the boundary -> create ghost in halo cell on the other side
                        if (isBoundaryCell(getCellCoords(p))) {
                            std::cout << "got here1\n";
                            if (getCellCoords(p)[0] == 1) {
                                double leftDistance = p.getX()[0] - meshSize;
                                // if (leftDistance <= reflectionDistance[p.getType()]) {
                                Particle par = Particle(
                                        {((dimensions[0] - 1) * meshSize + leftDistance), p.getX()[1], p.getX()[2]},
                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                ghosts.push_back(par);
                                //}
                            } else {
                                double rightDistance = (dimensions[0] - 1) * meshSize - p.getX()[0];
                                //if (rightDistance <= reflectionDistance[p.getType()]) {
                                Particle par = Particle({(meshSize - rightDistance), p.getX()[1], p.getX()[2]},
                                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(),
                                                        p.getEpsilon());
                                ghosts.push_back(par);
                                //}
                            }
                            //std::cout << "got here3\n";
                        }
                        // It's in the halo -> move it to the boundary on the other side
                        if (isHaloCell(getCellCoords(p))) {
                            bool rem = false;
                            std::cout << "got here2\n";
                            if (getCellCoords(p)[0] == 0) {
                                double leftDistance = abs(meshSize - p.getX()[0]);
                                //if (leftDistance <= reflectionDistance[p.getType()]) {
                                Particle part1 = Particle(
                                        {((dimensions[0] - 1) * meshSize - leftDistance), p.getX()[1], p.getX()[2]},
                                        p.getV(), p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                part1.setOldF(p.getOldF());
                                part1.setF(p.getF());
                                ins.push_back(part1);
                                rem = true;
                                //}
                            } else {
                                double rightDistance = abs(p.getX()[0] - (dimensions[0] - 1) * meshSize);
                                //if (rightDistance <= reflectionDistance[p.getType()]) {
                                Particle part2 = Particle({meshSize + rightDistance, p.getX()[1], p.getX()[2]},
                                                          p.getV(),
                                                          p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                part2.setOldF(p.getOldF());
                                part2.setF(p.getF());
                                ins.push_back(part2);
                                rem = true;
                                //}
                            }
                            if (rem)
                                del.push_back(p);
                        }
                    }
                    if (boundaryCondition[1] == 1) {
                        if (isBoundaryCell(getCellCoords(p))) {
                            double bottomDistance = p.getX()[1] - meshSize;
                            if (bottomDistance <= reflectionDistance[p.getType()]) {
                                Particle temp = Particle({p.getX()[0], (meshSize - bottomDistance), p.getX()[2]},
                                                         {0, 0, 0}, p.getM(), p.getType(), p.getSigma(),
                                                         p.getEpsilon());
                                p.setF(p.getF() + forceCalcs[p.getType()].calculateF(p, temp));
                            }
                            double topDistance = (dimensions[1] - 1) * meshSize - p.getX()[1];
                            if (topDistance <= reflectionDistance[p.getType()]) {
                                Particle temp = Particle(
                                        {(p.getX()[0]), (dimensions[1] - 1) * meshSize + topDistance, p.getX()[2]},
                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                p.setF(p.getF() + forceCalcs[p.getType()].calculateF(p, temp));
                            }
                        }
                    } else if (boundaryCondition[1] == 2) {
                        // It's in the boundary -> create ghost in halo cell on the other side
                        if (isBoundaryCell(getCellCoords(p))) {
                            if (getCellCoords(p)[1] == 1) {
                                double topDistance = (dimensions[1] - 1) * meshSize - p.getX()[1];
                                Particle par = Particle({p.getX()[0], meshSize - topDistance, p.getX()[2]}, {0, 0, 0},
                                                        p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                ghosts.push_back(par);
                            } else {
                                double bottomDistance = p.getX()[1] - meshSize;
                                Particle par = Particle(
                                        {p.getX()[0], ((dimensions[1] - 1) * meshSize + bottomDistance), p.getX()[2]},
                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                            }
                        } else {
                            bool rem = false;
                            if (getCellCoords(p)[1] == 1) {
                                double topDistance = abs(p.getX()[1] - (dimensions[1] - 1) * meshSize);
                                Particle part1 = Particle({p.getX()[0], meshSize + topDistance, p.getX()[2]}, p.getV(),
                                                          p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                part1.setF(p.getF());
                                part1.setOldF(p.getF());
                                ins.push_back(part1);
                                rem = true;
                            } else {
                                double bottomDistance = meshSize - p.getX()[1];
                                Particle part2 = Particle(
                                        {p.getX()[0], ((dimensions[1] - 1) * meshSize - bottomDistance), p.getX()[2]},
                                        p.getV(), p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                part2.setF(p.getF());
                                part2.setOldF(p.getOldF());
                                ins.push_back(part2);
                                rem = true;
                            }
                            if (rem)
                                del.push_back(p);
                        }
                    }
                    if (dimensions[2] != 0) {
                        if (boundaryCondition[2] == 1) {
                            if (isBoundaryCell(getCellCoords(p))) {
                                double frontDistance = p.getX()[2] - meshSize;
                                if (frontDistance <= reflectionDistance[p.getType()]) {
                                    Particle temp = Particle({p.getX()[0], p.getX()[1], (meshSize - frontDistance)},
                                                             {0, 0, 0}, p.getM(), p.getType(), p.getSigma(),
                                                             p.getEpsilon());
                                    p.setF(p.getF() + forceCalcs[p.getType()].calculateF(p, temp));
                                }
                                double backDistance = (dimensions[2] - 1) * meshSize - p.getX()[2];
                                if (backDistance <= reflectionDistance[p.getType()]) {
                                    Particle temp = Particle(
                                            {p.getX()[0], p.getX()[1], (dimensions[2] - 1) * meshSize + backDistance},
                                            {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon());
                                    p.setF(p.getF() + forceCalcs[p.getType()].calculateF(p, temp));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    deleteParticlesInHalo();
    std::cout << "got here4\n";
    for (auto par: ins) {
        forceInsert(par);
    }
    std::cout << "got here5\n";
    for (auto par: del) {
        remove(par);
    }
    std::cout << "got here6\n";
    //std::cout << "nr particles: " << numberParticles() << "\n";
    ghosts.clear();
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
                if (isHaloCell({x, y, 0})) {
                    std::list<Particle> newList;
                    particles[(coordToIndex({x, y, 0}))].clear();
                    particles.at(coordToIndex({x, y, 0})) = newList;
                    //std::cout << "Particle in halo detected \n";
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

double LinkedCells::calcKineticEnergy() {
    double energy = 0;
    for (auto &vec: particles) {
        for (auto &p: vec) {
            double scalar_product = p.getV()[0] * p.getV()[0] + p.getV()[1] * p.getV()[1] + p.getV()[2] * p.getV()[2];
            energy += p.getM() * scalar_product / 2;
        }
    }
    return energy;
}

void LinkedCells::scaleVelocity(double scale) {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setV(scale * p.getV());
        }
    }
}

void LinkedCells::applyGravity() {
    std::array<double, 3> temp = {0, 1, 0};
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setF(p.getF() + (p.getM() * grav) * temp);
        }
    }
}




