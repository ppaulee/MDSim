//
// Created by Jonas on 03.02.22.
//
#pragma once

#include <vector>
#include <list>
#include <memory>
#include "forceCalculation/MixedLennardJones.h"
#include "forceCalculation/LennardJones.h"
#include "Particle.h"
#include "forceCalculation/ForceCalculation.h"
#include "SimulationContainer.h"

#ifdef _OPENMP
#ifndef PSEMOLDYN_GROUPB_PARALLELLINKEDCELLS_H
#define PSEMOLDYN_GROUPB_PARALLELLINKEDCELLS_H


class ParallelLinkedCells : public SimulationContainer {

private:

    /**
     * Stores the Id of the next inserted particle
     */
    int currentId;

    /**
     * Stores used parallelisation strategy
     * 1 = blocking with locks
     * 2 = Threadbuffers
     */
    int strategy;

    /**
     * This vector stores cells. These cells store a vector of particles in this cell
     */
    std::vector<std::vector<Particle>> particles;

    /**
     * Stores LennardJones objects for force calculations between particles of the same type
     */
    std::vector<LennardJones> forceCalcs;
    /**
     * Stores LennardJones objects for force calculations between particles of different types
     */
    std::vector<MixedLennardJones> mixedForceCalcs;

    /**
     * Used for periodic boundaries, stores ghosts of particles in boundary cells that are mirrored into the halo on the other side
     */
    std::vector<Particle> ghosts;

    /**
     * Distance at which ghost particles are created
     */
    std::vector<double> reflectionDistance;

    /**
     * Stores which boundary condition is used for which boundary pair
     * First element for X-axis ("left and right")
     * Second element for Y-axis ("top and bottom")
     * Last element for Z-axis ("front and back")
     *
     * 0 = outflow
     * 1 = reflection
     * 2 = periodic
     */
    std::array<int, 3> boundaryCondition;

    /**
     * Width, Height, Depth
     */
    std::array<int, 3> dimensions;
    int size;

    /**
     * This is the length of the sides of one cube/square
     */
    double meshSize;
    double cutOffRadius;

    /**
     * Constant for gravity
     */
    double grav;

    /**
     * Translates a single index to coordinates
     *
     * @param coords Coordinates (x,y,z)
     * @return Index of the vector
     */
    std::array<int, 3> indexToCoords(int index);

    /**
     * Removes a particle from the data structure
     *
     * @param p Particle to remove
     * @param index optional: removes particle from given index
     */
    void remove(Particle &p, int index);

    /**
     * Sets all forces of particles to zero and sets oldF
     * Is called in calculateF()
     */
    void initCalculateF();

    /**
     * Does the same as forceInsert() but without giving a new ID to the particle and incrementing current_id
     *
     * @param p Particle to be inserted
     */
    void forceInsertNoId(Particle &p);


public:
    /**
     * Finds the right cell for the particle
     *
     * @param p Particle
     * @return Cell coordinates (x,y,z)
     */
    std::array<int, 3> getCellCoords(Particle &p);

    /**
     * Note that a layer of halo and boundary cells are added on the side. E.g. if we have a 1x1x1 cube then the actual dimensions are (1+4)x(1+4)x(1+4). The coordinate of the inner cell is (0+2,0+2,0+2)
     * General note: intern the smallest possible coordinate is (0,0,0) but we want to support also negative coordinates. Therefore to the position of the particle we add half the dimension in every dimension. The output will be retransferred to
     * the original coordinates (minus half the dimension in every dimension)
     *
     * @param dimension Dimensions (x,y,z) of the inner cells. All numbers must be even! -x/2,-y/2,-z/2 are the minimum coordinates and x/2,y/2,z/2 are the maximum coordinates
     * @param mesh Mesh size of the grid
     * @param cutOff cut off radius
     * @param gravConst Constant for gravitational force calculation (0 disables gravity)
     * @param boundaryConds Array setting boundary conditions
     */
    explicit ParallelLinkedCells(std::array<int, 3> dimension, double mesh, double cutOff, double gravConst,
                                 std::array<int, 3> &boundaryConds, int strat);

    /**
     * Note that a layer of halo and boundary cells are added on the side. E.g. if we have a 1x1x1 cube then the actual dimensions are (1+4)x(1+4)x(1+4). The coordinate of the inner cell is (0+2,0+2,0+2)
     * General note: intern the smallest possible coordinate is (0,0,0) but we want to support also negative coordinates. Therefore to the position of the particle we add half the dimension in every dimension. The output will be retransferred to
     * the original coordinates (minus half the dimension in every dimension). Using this constructor disables gravitation.
     *
     * @param dimension Dimensions (x,y,z) of the inner cells. All numbers must be even! -x/2,-y/2,-z/2 are the minimum coordinates and x/2,y/2,z/2 are the maximum coordinates
     * @param mesh Mesh size of the grid
     * @param cutOff cut off radius
     */
    explicit ParallelLinkedCells(std::array<int, 3> dimension, double mesh, double cutOff);

    /**
     * Translates coordinates to a single index for the vector
     *
     * @param coords Coordinates (x,y,z)
     * @return Index of the vector
     */
    int coordToIndex(std::array<int, 3> coords);

    /**
     *
     * @param coords Coordinates of the cell
     * @return Returns the content of a cell
     */
    std::vector<Particle> &get(std::array<double, 3> coords);

    /**
     * Same as get, but uses a "L" shaped coordinate system instead of a cross, meaning only indices >= 0 are allowed
     *
     * @param coords Coordinates of the cell
     * @return Returns the content of a cell
     */
    std::list<Particle> &absoluteGet(std::array<double, 3> coords);

    /**
     * Get cell using an index
     *
     * @param i index of a cell
     * @return The content of a cell
     */
    std::vector<Particle> &indexGet(int i);

    /**
     * Inserts a particle into the data structure
     *
     * @param p Particle to insert
     */
    void insert(Particle &p) override;

    /**
     * Inserts a particle into the data structure without changing its positional attribute
     *
     * @param p Particle to insert
     */
    void forceInsert(Particle &p) override;

    std::vector<Particle> getParticles() override;

    /**
     * Checks if the cell with these coordinates is a halo cell
     *
     * @param coords Coordinates
     * @return
     */
    bool isHaloCell(std::array<int, 3> coords);

    /**
     * Checks if the cell with these coordinates is a boundary cell
     *
     * @param coords Coordinates
     * @return
     */
    bool isBoundaryCell(std::array<int, 3> coords);

    /**
     * Multithreaded version of calculateF
     */
    void parallelCalculateF();

    /**
     * Checks if all particles are in the right cell
     *
     * @return
     */
    bool allInRightCell();

    /**
     * Calculates new velocity
     *
     * @param delta_t Size of timestep
     */
    void calculateV(double delta_t);

    /**
     * Calculates new position
     *
     * @param delta_t Size of timestep
     */
    void calculateX(double delta_t);

    /**
     * Moves particles which no longer belong to their cell to their corresponding cell. Will be called after in simulate()
     */
    void move();

    /**
     * Runs calculateX(), calculateF(), calculateV(), move() in order
     *
     * @param delta_a
     * @param algorithm
     */
    void simulate(ForceCalculation *algorithm, double delta_a) override;

    void applyForceBuffer();

    /**
     * Iterates over particle in boundary cells and applies boundary conditions according to boundaryCondition
     */
    void handleBoundary();

    /**
     * Plots particles
     *
     * @param iteration current Iteration
     */
    void plotParticles(int iteration) override;

    void addBrownianMotion(double averageV, int dimension) override;

    /**
     *
     * @return Number of particles inside the data structure
     */
    int numberParticles() override;

    /**
     * Deletess all particles in halo cells
     */
    void deleteParticlesInHalo();

    /**
     *
     * @return dimensions
     */
    std::array<int, 3> getDimensions();

    double calcKineticEnergy() override;

    void scaleVelocity(double scale) override;

    /**
     * Applies gravitational force to all particles
     */
    void applyGravity();

    bool isOutOfScope(Particle &p);

    std::array<double, 3> calculateLJForce(Particle &p1, Particle &p2);

    void calculateGhostForce(Particle &ghost);

    void handleReflectionBoundary(Particle &p);

    void setMembraneSimulation() override;

    void setNanoScaleFlowSimulation() override;

    void initMembrane() override;

    void simulateMembrane(double delta_t, bool pullState) override;
};


#endif //PSEMOLDYN_GROUPB_PARALLELLINKEDCELLS_H
#endif //_OPENMP