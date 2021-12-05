//
// Created by paulu on 21.11.2021.
//
#pragma once

#include <vector>
#include <list>
#include <memory>
#include "Particle.h"
#include "forceCalculation/ForceCalculation.h"
#include "SimulationContainer.h"



class LinkedCells : public SimulationContainer {

private:
    /**
     * This vector stores cells. These cells store a vector of particles in this cell
     */
    std::vector<std::list<Particle>> particles;

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
     * Finds the right cell for the particle
     *
     * @param p Particle
     * @return Cell coordinates (x,y,z)
     */
    std::array<int, 3> getCellCoords(Particle& p);

    /**
     * Translates coordinates to a single index for the vector
     *
     * @param coords Coordinates (x,y,z)
     * @return Index of the vector
     */
    int coordToIndex(std::array<int, 3> coords);

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
    void remove(Particle& p, int index);




public:

    /**
     * Note that a layer of halo and boundary cells are added on the side. E.g. if we have a 1x1x1 cube then the actual dimensions are (1+4)x(1+4)x(1+4). The coordinate of the inner cell is (0+2,0+2,0+2)
     * General note: intern the smallest possible coordinate is (0,0,0) but we want to support also negative coordinates. Therefore to the position of the particle we add half the dimension in every dimension. The output will be retransferred to
     * the original coordinates (minus half the dimension in every dimension)
     *
     * @param dimension Dimensions (x,y,z) of the inner cells. All numbers must be even! -x/2,-y/2,-z/2 are the minimum coordinates and x/2,y/2,z/2 are the maximum coordinates
     * @param mesh Mesh size of the grid
     * @param cutOff cut off radius
     */
    explicit LinkedCells(std::array<int, 3> dimension, double mesh, double cutOff);

    /**
     *
     * @param coords Coordinates of the cell
     * @return Returns the content of a cell
     */
    std::list<Particle>& get(std::array<double, 3> coords);

    /**
     * Inserts a particle into the data structure
     *
     * @param p Particle to insert
     */
    void insert(Particle& p) override;

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
     * Calculates the force between particles with a cutoff range
     *
     * @param algorithm Algorithm to compute the force
     */
    void calculateF(ForceCalculation *algorithm);


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
    int numberParticles();

    /**
     * Deletess all particles in halo cells
     */
    void deleteParticlesInHalo();

    /**
     * Deletes all particles out of scope
     */
    void deleteOutside();

    /**
     *
     * @return dimensions
     */
    std::array<int, 3> getDimensions();

    /**
     * Transfers the coordiantes to a format which we can work with
     * We add half the dimension on each dimension
     *
     * @param p Particle which should be transferred
     * @return Transferred Particle
     */
    Particle& transfer(Particle& p);

    /**
     * Transferes the particle back
     * See transform() for transformation rule
     *
     * @param p Particle which should be transferred
     * @return Retransferred Particle
     */
    Particle& retransfer(Particle& p);

};


