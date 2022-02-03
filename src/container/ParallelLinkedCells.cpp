//
// Created by Jonas on 03.02.22.
//

#include "ParallelLinkedCells.h"

#include "utils/ArrayUtils.h"
#include <outputWriter/VTKWriter.h>
#include <stdexcept>
#include <iostream>
#include <omp.h>

#ifdef _OPENMP

ParallelLinkedCells::ParallelLinkedCells(std::array<int, 3> dimension, double mesh, double cutOff, double gravConst,
                                         std::array<int, 3> &boundaryConds, int strat) {
    if (dimension[0] % 2 != 0 || dimension[1] % 2 != 0 || dimension[2] % 2 != 0) {
        throw std::invalid_argument("Dimensions must be even");
    }
    grav = gravConst;
    boundaryCondition = boundaryConds;
    cutOffRadius = cutOff;
    currentId = 0;
    strategy = strat;
    // Add boundary and halo cells
    if (dimension[2] == 0) {
        // 2D case
        dimensions = {dimension[0] + 4, dimension[1] + 4, 0};
        size = (dimension[0] + 4) * (dimension[1] + 4) - 1;
    } else {
        // 3D Case
        dimensions = {dimension[0] + 4, dimension[1] + 4, dimension[2] + 4};
        size = coordToIndex({dimension[0] + 4, dimension[1] + 4, dimension[2] + 4});
    }

    meshSize = mesh;
    // Initialise all vectors
    for (int i = 0; i < coordToIndex({dimensions[0], dimensions[1], dimensions[2]}) + 1; i++) {
        particles.push_back(std::vector<Particle>({}));
    }
}

ParallelLinkedCells::ParallelLinkedCells(std::array<int, 3> dimension, double mesh, double cutOff) {
    if (dimension[0] % 2 != 0 || dimension[1] % 2 != 0 || dimension[2] % 2 != 0) {
        throw std::invalid_argument("Dimensions must be even");
    }
    grav = 0;
    boundaryCondition = {0, 0, 0};
    cutOffRadius = cutOff;
    currentId = 0;
    strategy = 1;
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
        particles.push_back(std::vector<Particle>({}));
    }
}

int ParallelLinkedCells::coordToIndex(std::array<int, 3> coords) {
    if (dimensions[2] == 0) {
        return (coords[1]) * dimensions[0] + (coords[0]);
    } else {
        //TODO fix this for 3D case
        return coords[0] * dimensions[1] * dimensions[2] + coords[1] * dimensions[2] + coords[2];
    }
}

// Only 3D
std::array<int, 3> ParallelLinkedCells::indexToCoords(int index) {
    int width_index = index / (dimensions[1] * dimensions[2]);
    int height_index = (index - width_index * dimensions[1] * dimensions[2]) / dimensions[2];
    int depth_index = index - width_index * dimensions[1] * dimensions[2] - height_index * dimensions[2];
    return {width_index, height_index, depth_index};
}

void ParallelLinkedCells::insert(Particle &p) {
    p.setID(currentId);
    currentId++;
    if (strategy == 2)
        p.initParallelBuffer(omp_get_max_threads());
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
                    auto temp = MixedLennardJones(tempEpsilon, tempSigma, calc.getType(), p.getType());
                    mixedForceCalcs.push_back(temp);
                }
            }
        }
    }
    particles.at(coordToIndex(getCellCoords(p))).push_back(p);
}

void ParallelLinkedCells::forceInsert(Particle &p) {
    p.setID(currentId);
    currentId++;
    if (strategy == 2)
        p.initParallelBuffer(omp_get_max_threads());
    if (forceCalcs.size() == p.getType()) {
        forceCalcs.push_back(LennardJones(p.getEpsilon(), p.getSigma(), p.getType()));
        reflectionDistance.push_back((cbrt(sqrt(2)) * p.getSigma()) / 2);
        if (forceCalcs.size() >= 1) {
            for (auto &calc: forceCalcs) {
                if (calc.getType() != p.getType()) {
                    double tempEpsilon = sqrt(calc.getEpsilon() * p.getEpsilon());
                    double tempSigma = (calc.getSigma() + p.getSigma()) / 2;
                    auto temp = MixedLennardJones(tempEpsilon, tempSigma, calc.getType(), p.getType());
                    mixedForceCalcs.push_back(temp);
                }
            }
        }
    }
    particles.at(coordToIndex(getCellCoords(p))).push_back(p);
}

std::array<int, 3> ParallelLinkedCells::getCellCoords(Particle &p) {
    int x = ((int) (p.getX()[0] / meshSize));
    int y = ((int) (p.getX()[1] / meshSize));
    int z = ((int) (p.getX()[2] / meshSize));
    return {x, y, z};
}

void ParallelLinkedCells::remove(Particle &p, int index = -1) {
    int i = index;
    if (index == -1) {
        i = coordToIndex(getCellCoords(p));
    }
    for (std::vector<Particle>::iterator it = particles[i].begin(); it != particles[i].end(); ++it) {
        if (p.getID() == it->getID()) {
            particles.at(i).erase(it);
            break;
        }
    }
}

bool ParallelLinkedCells::isHaloCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0] - 1 || coords[1] == dimensions[1] - 1 || coords[2] == dimensions[2] - 1 ||
           coords[0] == 0 || coords[1] == 0 || (coords[2] == 0 && dimensions[2] != 0);
}

bool ParallelLinkedCells::isBoundaryCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0] - 2 || coords[1] == dimensions[1] - 2 || coords[2] == dimensions[2] - 2 ||
           coords[0] == 1 || coords[1] == 1 || coords[2] == 1;
}

std::vector<Particle> &ParallelLinkedCells::get(std::array<double, 3> coords) {
    int x = ((int) (coords[0] / meshSize));
    int y = ((int) (coords[1] / meshSize));
    int z = ((int) (coords[2] / meshSize));
    return particles[coordToIndex({x, y, z})];
}

std::vector<Particle> &ParallelLinkedCells::indexGet(int i) {
    return particles[i];
}

std::array<double, 3> ParallelLinkedCells::calculateLJForce(Particle &p1, Particle &p2) {
    std::array<double, 3> force;
    if (p2.getType() == p1.getType()) {
        auto calc = LennardJones(p1.getEpsilon(), p1.getSigma(), p1.getType());
        force = calc.calculateF(p1, p2);
    } else {
        auto calc = MixedLennardJones(sqrt(p1.getEpsilon() * p2.getEpsilon()), (p1.getSigma() + p2.getSigma()) / 2,
                                      p1.getType(), p2.getType());
        force = calc.calculateF(p1, p2);
    }
    return force;
}

void ParallelLinkedCells::parallelCalculateF() {
    initCalculateF();
#pragma omp parallel for shared(particles) schedule(static)
    for (int x = 0; x < dimensions[0]; ++x) {
        for (int y = 0; y < dimensions[1]; ++y) {
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
                            if (c[0] < 0 || c[0] >= dimensions[0] || c[1] < 0 || c[1] >= dimensions[1]) {
                                continue;
                            }
                            if (coordToIndex(c) < 0) {
                                continue;
                            }
                            // Particles in neighboured cells
                            for (auto &p: particles[coordToIndex(c)]) {
                                // Check if particle is inside the cut-off radius
                                //p.setCalcLock();
                                //current_particle.setCalcLock();
                                //current_particle.setLock();
                                std::array<double, 3> current_X = current_particle.getX();
                                //current_particle.unSetLock();
                                //p.setLock();
                                double distance = ArrayUtils::L2Norm(p.getX() - current_X);
                                //p.unSetLock();
                                if (distance <= cutOffRadius &&
                                    !(p == current_particle) && !p.isMarked()) {
                                    //p.unSetCalcLock();
                                    // Actual force calculation according to the used algorithm
                                    std::array<double, 3> force;
                                    force = calculateLJForce(current_particle, p);

                                    // Make use of Newtons third law
                                    if (strategy == 1) {
                                        current_particle.setLock();
                                        current_particle.setF(current_particle.getF() + force);
                                        current_particle.unSetLock();
                                        p.setLock();
                                        p.setF(p.getF() - force);
                                        p.unSetLock();
                                    } else if (strategy == 2) {
                                        int threadNum = omp_get_thread_num();
                                        current_particle.setParallelForce(threadNum,
                                                                          current_particle.getParallelForce(threadNum) +
                                                                          force);
                                        p.setParallelForce(threadNum, p.getParallelForce(threadNum) - force);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                for (int z = 0; z < dimensions[2] - 1; ++z) {
                    // Loop over all particles in cell
                    for (auto &current_particle: particles[coordToIndex({x, y, z})]) {
                        current_particle.mark();
                        // Get neighbours
                        for (int x_diff = -1; x_diff <= 1; ++x_diff) {
                            for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                                for (int z_diff = -1; z_diff <= 1; ++z_diff) {

                                    // Coordinates including offset for neighboured cells
                                    std::array<int, 3> c = {x + x_diff, y + y_diff, z + z_diff};

                                    if (c[0] < 0 || c[0] >= dimensions[0] || c[1] < 0 || c[1] >= dimensions[1] ||
                                        c[2] < 0 || c[2] >= dimensions[2]) {
                                        continue;
                                    }
                                    if (coordToIndex(c) < 0) {
                                        continue;
                                    }
                                    for (auto &p: particles[coordToIndex(c)]) {
                                        // Check if particle is inside the cut-off radius
                                        if (ArrayUtils::L2Norm(p.getX() - current_particle.getX()) <= cutOffRadius &&
                                            !(p == current_particle) &&
                                            !p.isMarked()) {
                                            // Actual force calculation according to the used algorithm
                                            std::array<double, 3> force;
                                            force = calculateLJForce(current_particle, p);

                                            // Make use of Newtons third law
                                            if (strategy == 1) {
                                                current_particle.setLock();
                                                current_particle.setF(current_particle.getF() + force);
                                                current_particle.unSetLock();
                                                p.setLock();
                                                p.setF(p.getF() - force);
                                                p.unSetLock();
                                            } else if (strategy == 2) {
                                                int threadNum = omp_get_thread_num();
                                                current_particle.setParallelForce(threadNum,
                                                                                  current_particle.getParallelForce(
                                                                                          threadNum) + force);
                                                p.setParallelForce(threadNum, p.getParallelForce(threadNum) - force);
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
    // add up buffers from the threads
    if (strategy == 2) {
#pragma omp parallel for schedule(static)
        for (auto &cell: particles) {
            for (auto &p: cell) {
                for (int i = 0; i < omp_get_max_threads(); i++)
                    p.setF(p.getF() + p.getParallelForce(i));
            }
        }
    }
    ghosts.clear();
}

void ParallelLinkedCells::initCalculateF() {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setOldF(p.getF());
            p.setF({0, 0, 0});
            p.unmark();
            if (strategy == 2) {
                std::array<double, 3> temp = {0, 0, 0};
                for (int i = 0; i < omp_get_max_threads(); i++)
                    p.setParallelForce(i, temp);
            }
            // p.getCalculated().clear();
        }
    }
}

void ParallelLinkedCells::calculateV(double delta_t) {
#pragma omp parallel for
    for (auto &vec: particles) {
        for (auto &p: vec) {
            std::array<double, 3> vBefore = p.getV();
            std::array<double, 3> res = p.getF() + p.getOldF();
            res = (1 / (2 * p.getM())) * res;
            res = delta_t * res;
            res = res + p.getV();
            p.setV(res);
        }
    }
}

void ParallelLinkedCells::calculateX(double delta_t) {
#pragma omp parallel for
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

void ParallelLinkedCells::move() {
    int numberStart = numberParticles();
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
        remove(to_remove, index);
        // insert
        auto &p = tmp_stack_insert[i];
        std::array<int, 3> c = getCellCoords(p);
        if (c[0] < 0 || c[0] > dimensions[0] || c[1] < 0 || c[1] > dimensions[1]) {
            continue;
        }
        // if (coordToIndex(getCellCoords(p)) >= 0 &&
        //    (long unsigned int) coordToIndex(getCellCoords(p)) < particles.size()) {
        //if (!isOutOfScope(p))
        particles.at(coordToIndex(getCellCoords(p))).push_back(p);
        // }
    }
    //if(numberParticles() < 4)
}

bool ParallelLinkedCells::allInRightCell() {
    for (int i = 0; i < size; ++i) {
        for (auto &p: particles[i]) {
            if (coordToIndex(getCellCoords(p)) != i) {
                return false;
            }
        }
    }
    return true;
}

void ParallelLinkedCells::simulate(ForceCalculation *algorithm, double delta_t) {
    // calculate new x
    calculateX(delta_t);
    //deleteParticlesInHalo();
    // moves particles to correct cell
    move();
    handleBoundary();
    move();
    // calculate new f
    parallelCalculateF();
    applyForceBuffer();

    //deleteParticlesInHalo();

    if (grav != 0) {
        applyGravity();
    }
    // calculate new v
    calculateV(delta_t);
}

void ParallelLinkedCells::applyForceBuffer() {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setF(p.getF() + p.getForceBuffer());
            p.setForceBuffer({0, 0, 0});
        }
    }
}

void ParallelLinkedCells::calculateGhostForce(Particle &ghost) {
    std::array<int, 3> cell = getCellCoords(ghost);
    if (dimensions[2] == 0) {
        for (int x_diff = -1; x_diff <= 1; ++x_diff) {
            for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                std::array<int, 3> c = {cell[0] + x_diff, cell[1] + y_diff, 0};
                if (c[0] < 0 || c[0] >= dimensions[0] || c[1] < 0 || c[1] >= dimensions[1])
                    continue;
                for (auto &p: particles[coordToIndex(c)]) {
                    double distance = ArrayUtils::L2Norm(p.getX() - ghost.getX());
                    if (distance > 0 && distance <= cutOffRadius) {
                        std::array<double, 3> force;
                        force = calculateLJForce(p, ghost);
                        p.setForceBuffer(p.getForceBuffer() + force);
                    }
                }
            }
        }
    } else {
        for (int x_diff = -1; x_diff <= 1; ++x_diff) {
            for (int y_diff = -1; y_diff <= 1; ++y_diff) {
                for (int z_diff = -1; z_diff <= 1; ++z_diff) {
                    std::array<int, 3> c = {cell[0] + x_diff, cell[1] + y_diff, cell[2] + z_diff};
                    if (c[0] < 0 || c[0] >= dimensions[0] || c[1] < 0 || c[1] >= dimensions[1] ||
                        c[2] < 0 || c[2] >= dimensions[2])
                        continue;
                    for (auto &p: particles[coordToIndex(c)]) {
                        double distance = ArrayUtils::L2Norm(p.getX() - ghost.getX());
                        if (distance > 0 && distance <= cutOffRadius) {
                            std::array<double, 3> force;
                            force = calculateLJForce(p, ghost);
                            std::array<double, 3> fBefore = p.getForceBuffer();
                            p.setForceBuffer(p.getForceBuffer() + force);
                        }
                    }
                }
            }
        }
    }
}

void ParallelLinkedCells::handleReflectionBoundary(Particle &p) {
    //for (auto &vec: particles) {
    //if (!vec.empty()) {
    //if (isBoundaryCell(getCellCoords(vec.front())) || isHaloCell(getCellCoords(vec.front()))) {
    //for (auto &p: vec) {
    if (boundaryCondition[0] == 1) {
        if (isBoundaryCell(getCellCoords(p)) || isHaloCell(getCellCoords(p))) {
            double leftDistance = p.getX()[0] - meshSize;
            if (leftDistance <= reflectionDistance[p.getType()]) {
                Particle temp = Particle(
                        {(meshSize - leftDistance), p.getX()[1], p.getX()[2]},
                        {0, 0, 0},
                        p.getM(),
                        p.getType(), p.getSigma(), p.getEpsilon(), p.getID());
                p.setForceBuffer(
                        p.getForceBuffer() +
                        forceCalcs[p.getType()].calculateF(p, temp));
                //ghosts.push_back(temp);

            }
            double rightDistance = (dimensions[0] - 1) * meshSize - p.getX()[0];
            if (rightDistance <= reflectionDistance[p.getType()]) {
                Particle temp = Particle(
                        {((dimensions[0] - 1) * meshSize + rightDistance), p.getX()[1],
                         p.getX()[2]},
                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                        p.getID());
                p.setForceBuffer(
                        p.getForceBuffer() +
                        forceCalcs[p.getType()].calculateF(p, temp));
                //ghosts.push_back(temp);
            }
        }
    }
    if (boundaryCondition[1] == 1) {
        if (isBoundaryCell(getCellCoords(p)) || isHaloCell(getCellCoords(p))) {
            double bottomDistance = p.getX()[1] - meshSize;
            if (bottomDistance <= reflectionDistance[p.getType()]) {
                Particle temp = Particle({p.getX()[0], (meshSize - bottomDistance), p.getX()[2]},
                                         {0, 0, 0}, p.getM(), p.getType(), p.getSigma(),
                                         p.getEpsilon(), p.getID());
                p.setForceBuffer(p.getForceBuffer() + forceCalcs[p.getType()].calculateF(p, temp));
                //ghosts.push_back(temp);
            }
            double topDistance = (dimensions[1] - 1) * meshSize - p.getX()[1];
            if (topDistance <= reflectionDistance[p.getType()]) {
                Particle temp = Particle(
                        {(p.getX()[0]), (dimensions[1] - 1) * meshSize + topDistance, p.getX()[2]},
                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon(), p.getID());
                p.setForceBuffer(p.getForceBuffer() + forceCalcs[p.getType()].calculateF(p, temp));
                //ghosts.push_back(temp);
            }
        }
    }
    if (dimensions[2] != 0) {
        if (boundaryCondition[2] == 1) {
            if (isBoundaryCell(getCellCoords(p))) {
                double frontDistance = p.getX()[2] - meshSize;
                if (frontDistance <= reflectionDistance[p.getType()]) {
                    Particle temp = Particle({p.getX()[0], p.getX()[1], (meshSize - frontDistance)},
                                             {0, 0, 0}, p.getM(), p.getType(), p.getSigma(),
                                             p.getEpsilon(), p.getID());
                    p.setForceBuffer(
                            p.getForceBuffer() + forceCalcs[p.getType()].calculateF(p, temp));
                }
                double backDistance = (dimensions[2] - 1) * meshSize - p.getX()[2];
                if (backDistance <= reflectionDistance[p.getType()]) {
                    Particle temp = Particle(
                            {p.getX()[0], p.getX()[1],
                             (dimensions[2] - 1) * meshSize + backDistance},
                            {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                            p.getID());
                    p.setForceBuffer(
                            p.getForceBuffer() + forceCalcs[p.getType()].calculateF(p, temp));
                }
            }
        }
    }
}
//}
//}
//}
//}

void ParallelLinkedCells::handleBoundary() {
    int numberStart = numberParticles();
    std::vector<Particle> ins;
    std::vector<Particle> del;
    for (auto &vec: particles) {
        if (!vec.empty()) {
            if (isBoundaryCell(getCellCoords(vec.front())) || isHaloCell(getCellCoords(vec.front()))) {
                for (auto &p: vec) {
                    handleReflectionBoundary(p);
                    bool corner = false;
                    if (dimensions[2] != 0) {
                        bool rem = false;
                        if (boundaryCondition[0] == 2 && boundaryCondition[2] == 2) {
                            if (p.getX()[0] == 0 && p.getX()[2] == 0) {
                                //front left corner
                                double leftDistance = meshSize - p.getX()[0];
                                double frontDistance = meshSize - p.getX()[2];
                                Particle par = Particle(
                                        {((dimensions[0] - 1) * meshSize - leftDistance), p.getX()[1],
                                         ((dimensions[2] - 1) * meshSize - frontDistance)},
                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon(), p.getID());
                                par.setOldF(p.getOldF());
                                par.setF(p.getF());
                                par.setForceBuffer(p.getForceBuffer());
                                //part1.setCalculated(p.getCalculated());
                                ins.push_back(par);
                                rem = true;
                                corner = true;
                            } else if (p.getX()[0] == dimensions[0] - 1 && p.getX()[2] == 0) {
                                //front right corner
                                double rightDistance = p.getX()[0] - (dimensions[0] - 1) * meshSize;
                                double frontDistance = meshSize - p.getX()[2];
                                Particle par = Particle(
                                        {meshSize + rightDistance, p.getX()[1],
                                         ((dimensions[2] - 1) * meshSize - frontDistance)},
                                        p.getV(),
                                        p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                        p.getID());
                                par.setOldF(p.getOldF());
                                par.setF(p.getF());
                                par.setForceBuffer(p.getForceBuffer());
                                // part2.setCalculated(p.getCalculated());
                                ins.push_back(par);
                                rem = true;
                                corner = true;
                            } else if (p.getX()[0] == 0 && p.getX()[2] == dimensions[2] - 1) {
                                //back left corner
                                double leftDistance = meshSize - p.getX()[0];
                                double backDistance = p.getX()[2] - (dimensions[2] - 1) * meshSize;
                                Particle par = Particle({((dimensions[0] - 1) * meshSize - leftDistance), p.getX()[1],
                                                         meshSize + backDistance},
                                                        p.getV(),
                                                        p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                                        p.getID());
                                par.setOldF(p.getOldF());
                                par.setF(p.getF());
                                par.setForceBuffer(p.getForceBuffer());
                                // part2.setCalculated(p.getCalculated());
                                ins.push_back(par);
                                rem = true;
                                corner = true;
                            } else if (p.getX()[0] == dimensions[0] - 1 && p.getX()[2] == dimensions[2] - 1) {
                                //back right corner
                                double rightDistance = p.getX()[0] - (dimensions[0] - 1) * meshSize;
                                double backDistance = p.getX()[2] - (dimensions[2] - 1) * meshSize;
                                Particle par = Particle(
                                        {meshSize + rightDistance, p.getX()[1], meshSize + backDistance},
                                        p.getV(),
                                        p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                        p.getID());
                                par.setOldF(p.getOldF());
                                par.setF(p.getF());
                                par.setForceBuffer(p.getForceBuffer());
                                // part2.setCalculated(p.getCalculated());
                                ins.push_back(par);
                                rem = true;
                                corner = true;
                            }
                            if (rem)
                                del.push_back(p);
                        }
                    }
                    if (!corner && boundaryCondition[0] == 2) {
                        // It's in the boundary -> create ghost in halo cell on the other side
                        //if (isBoundaryCell(getCellCoords(p))) {
                        if (getCellCoords(p)[0] == 1) {
                            double leftDistance = p.getX()[0] - meshSize;
                            // if (leftDistance <= reflectionDistance[p.getType()]) {
                            Particle par = Particle(
                                    {((dimensions[0] - 1) * meshSize + leftDistance), p.getX()[1], p.getX()[2]},
                                    {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon(), p.getID());
                            calculateGhostForce(par);
                            //ghosts.push_back(par);
                            //}
                        } else if (getCellCoords(p)[0] == dimensions[0] - 2) {
                            double rightDistance = (dimensions[0] - 1) * meshSize - p.getX()[0];
                            //if (rightDistance <= reflectionDistance[p.getType()]) {
                            Particle par = Particle({(meshSize - rightDistance), p.getX()[1], p.getX()[2]},
                                                    {0, 0, 0}, p.getM(), p.getType(), p.getSigma(),
                                                    p.getEpsilon(), p.getID());
                            calculateGhostForce(par);
                            //ghosts.push_back(par);
                            //}
                        }
                        //}
                        // It's in the halo -> move it to the boundary on the other side
                        //else if (isHaloCell(getCellCoords(p))) {
                        bool rem = false;
                        if (getCellCoords(p)[0] == 0) {
                            double leftDistance = meshSize - p.getX()[0];
                            //if (leftDistance <= reflectionDistance[p.getType()]) {
                            Particle part1 = Particle(
                                    {((dimensions[0] - 1) * meshSize - leftDistance), p.getX()[1], p.getX()[2]},
                                    p.getV(), p.getM(), p.getType(), p.getSigma(), p.getEpsilon(), p.getID());
                            part1.setOldF(p.getOldF());
                            part1.setF(p.getF());
                            part1.setForceBuffer(p.getForceBuffer());
                            //part1.setCalculated(p.getCalculated());
                            ins.push_back(part1);
                            rem = true;
                            //}
                        } else if (getCellCoords(p)[0] == dimensions[0] - 1) {
                            double rightDistance = p.getX()[0] - (dimensions[0] - 1) * meshSize;
                            //if (rightDistance <= reflectionDistance[p.getType()]) {
                            Particle part2 = Particle({meshSize + rightDistance, p.getX()[1], p.getX()[2]},
                                                      p.getV(),
                                                      p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                                      p.getID());
                            part2.setOldF(p.getOldF());
                            part2.setF(p.getF());
                            part2.setForceBuffer(p.getForceBuffer());
                            // part2.setCalculated(p.getCalculated());
                            ins.push_back(part2);
                            rem = true;
                            //}
                        }
                        if (rem)
                            del.push_back(p);
                        //}
                    }
                    if (boundaryCondition[1] == 2) {
                        // It's in the boundary -> create ghost in halo cell on the other side
                        if (isBoundaryCell(getCellCoords(p))) {
                            if (getCellCoords(p)[1] == 1) {
                                double topDistance = (dimensions[1] - 1) * meshSize - p.getX()[1];
                                Particle par = Particle({p.getX()[0], meshSize - topDistance, p.getX()[2]},
                                                        {0, 0, 0},
                                                        p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                                        p.getID());
                                //ghosts.push_back(par);
                                calculateGhostForce(par);
                            } else {
                                double bottomDistance = p.getX()[1] - meshSize;
                                Particle par = Particle(
                                        {p.getX()[0], ((dimensions[1] - 1) * meshSize + bottomDistance),
                                         p.getX()[2]},
                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon(), p.getID());
                                //ghosts.push_back(par);
                                calculateGhostForce(par);
                            }
                        } else if (isHaloCell(getCellCoords(p))) {
                            bool rem = false;
                            if (getCellCoords(p)[1] == 0) {
                                double topDistance = abs(p.getX()[1] - (dimensions[1] - 1) * meshSize);
                                Particle part1 = Particle({p.getX()[0], meshSize + topDistance, p.getX()[2]},
                                                          p.getV(),
                                                          p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                                          p.getID());
                                part1.setF(p.getF());
                                part1.setOldF(p.getF());
                                ins.push_back(part1);
                                rem = true;
                            } else {
                                double bottomDistance = meshSize - p.getX()[1];
                                Particle part2 = Particle(
                                        {p.getX()[0], ((dimensions[1] - 1) * meshSize - bottomDistance),
                                         p.getX()[2]},
                                        p.getV(), p.getM(), p.getType(), p.getSigma(), p.getEpsilon(), p.getID());
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
                        if (!corner && boundaryCondition[2] == 2) {
                            // It's in the boundary -> create ghost in halo cell on the other side
                            //if (isBoundaryCell(getCellCoords(p))) {
                            if (getCellCoords(p)[2] == 1) {
                                double frontDistance = p.getX()[2] - meshSize;
                                // if (leftDistance <= reflectionDistance[p.getType()]) {
                                Particle par = Particle(
                                        {p.getX()[0], p.getX()[1],
                                         ((dimensions[2] - 1) * meshSize + frontDistance)},
                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                        p.getID());
                                //ghosts.push_back(par);
                                calculateGhostForce(par);
                                //}
                            } else if (getCellCoords(p)[2] == dimensions[2] - 2) {
                                double backDistance = (dimensions[2] - 1) * meshSize - p.getX()[2];
                                //if (rightDistance <= reflectionDistance[p.getType()]) {
                                Particle par = Particle({p.getX()[0], p.getX()[1], meshSize - backDistance},
                                                        {0, 0, 0}, p.getM(), p.getType(), p.getSigma(),
                                                        p.getEpsilon(), p.getID());
                                //ghosts.push_back(par);
                                calculateGhostForce(par);
                                //}
                            }
                            //}
                            // It's in the halo -> move it to the boundary on the other side
                            //else if (isHaloCell(getCellCoords(p))) {
                            bool rem = false;
                            if (getCellCoords(p)[2] == 0) {
                                double frontDistance = meshSize - p.getX()[2];
                                //if (leftDistance <= reflectionDistance[p.getType()]) {
                                Particle part1 = Particle(
                                        {p.getX()[0], p.getX()[1],
                                         ((dimensions[2] - 1) * meshSize - frontDistance)},
                                        p.getV(), p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                        p.getID());
                                part1.setOldF(p.getOldF());
                                part1.setF(p.getF());
                                part1.setForceBuffer(p.getForceBuffer());
                                //part1.setCalculated(p.getCalculated());
                                ins.push_back(part1);
                                rem = true;
                                //}
                            } else if (getCellCoords(p)[2] == dimensions[2] - 1) {
                                double backDistance = p.getX()[2] - (dimensions[2] - 1) * meshSize;
                                //if (rightDistance <= reflectionDistance[p.getType()]) {
                                Particle part2 = Particle({p.getX()[0], p.getX()[1], meshSize + backDistance},
                                                          p.getV(),
                                                          p.getM(), p.getType(), p.getSigma(), p.getEpsilon(),
                                                          p.getID());
                                part2.setOldF(p.getOldF());
                                part2.setF(p.getF());
                                part2.setForceBuffer(p.getForceBuffer());
                                // part2.setCalculated(p.getCalculated());
                                ins.push_back(part2);
                                rem = true;
                                //}
                            }
                            if (rem)
                                del.push_back(p);
                            //}
                        }
                    }
                }
            }
        }
    }
    deleteParticlesInHalo();
    for (auto par: ins) {
        forceInsertNoId(par);
    }
    move();
    ins.clear();
    del.clear();
}

void ParallelLinkedCells::plotParticles(int iteration) {
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

void ParallelLinkedCells::addBrownianMotion(double averageV, int dimension) {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setV(p.getV() + maxwellBoltzmannDistributedVelocity(averageV, dimension));
        }
    }
}

int ParallelLinkedCells::numberParticles() {
    int result = 0;
    for (auto &vec: particles) {
        result += vec.size();
    }
    return result;
}

void ParallelLinkedCells::deleteParticlesInHalo() {
    for (int x = 0; x < dimensions[0]; x++) {
        for (int y = 0; y < dimensions[1]; y++) {
            if (dimensions[2] == 0) {
                if (isHaloCell({x, y, 0}) && !particles[(coordToIndex({x, y, 0}))].empty()) {
                    std::list<Particle> newList;
                    particles[(coordToIndex({x, y, 0}))].clear();
                    //particles.at(coordToIndex({x, y, 0})) = newList;
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

bool ParallelLinkedCells::isOutOfScope(Particle &p) {
    if (p.getX()[0] < 0 || p.getX()[0] > dimensions[0] * meshSize || p.getX()[1] < 0 ||
        p.getX()[1] > dimensions[1] * meshSize)
        return true;
    if (dimensions[2] != 0) {
        //TODO 3D case
    }
    return false;
}

std::array<int, 3> ParallelLinkedCells::getDimensions() {
    return dimensions;
}

double ParallelLinkedCells::calcKineticEnergy() {
    double energy = 0;
    for (auto &vec: particles) {
        for (auto &p: vec) {
            double scalar_product =
                    p.getV()[0] * p.getV()[0] + p.getV()[1] * p.getV()[1] + p.getV()[2] * p.getV()[2];
            energy += p.getM() * scalar_product / 2;
        }
    }
    return energy;
}

void ParallelLinkedCells::scaleVelocity(double scale) {
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setV(scale * p.getV());
        }
    }
}

void ParallelLinkedCells::applyGravity() {
    std::array<double, 3> temp = {0, 1, 0};
    for (auto &vec: particles) {
        for (auto &p: vec) {
            p.setF(p.getF() + (p.getM() * grav) * temp);
        }
    }
}

std::vector<Particle> ParallelLinkedCells::getParticles() {
    std::vector<Particle> res;
    for (auto &vec: particles) {
        for (auto &p: vec) {
            res.push_back(p);
        }
    }
    return res;
}

void ParallelLinkedCells::forceInsertNoId(Particle &p) {
    /*if (forceCalcs.size() == p.getType()) {
        forceCalcs.push_back(LennardJones(p.getEpsilon(), p.getSigma(), p.getType()));
        reflectionDistance.push_back((cbrt(sqrt(2)) * p.getSigma()) / 2);
        if (forceCalcs.size() >= 1) {
            for (auto &calc: forceCalcs) {
                if (calc.getType() != p.getType()) {
                    double tempEpsilon = sqrt(calc.getEpsilon() * p.getEpsilon());
                    double tempSigma = (calc.getSigma() + p.getSigma()) / 2;
                    auto temp = MixedLennardJones(tempEpsilon, tempSigma, calc.getType(), p.getType());
                    mixedForceCalcs.push_back(temp);
                }
            }
        }
    }*/
    if (strategy == 2)
        p.initParallelBuffer(omp_get_max_threads());
    if (!isOutOfScope(p))
        particles[coordToIndex(getCellCoords(p))].push_back(p);
}

void ParallelLinkedCells::setMembraneSimulation() {
    return;

}

void ParallelLinkedCells::setNanoScaleFlowSimulation() {
    return;
}

void ParallelLinkedCells::initMembrane() {
    return;
}

void ParallelLinkedCells::simulateMembrane(double delta_t, bool pullState) {
    return;
}

#endif //_OPENMP