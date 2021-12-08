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


/**
 * \image html benchmark_cmp_algorithms.png
 * This benachmark was tested on a Ryzen 7 1700 @3.0 GHz with 16GB RAM
 *
 * Here we can see the performance difference between the different algorithms. The LinkedCells algorithm is encoded as 1 and the naive force calculation is encoded as 0.
 *
 * We can create a linear regression model with this data. In fact we can observe a quadratic relationship therefore, we use a linear model with a polynomial with degree 2. (time = a_0 + a_1*numberOfParticles + a_2*numberOfParticles^2 + e   where e ~ N(0,1)) Below the regression model for the naive approach is plotted.
 * \image html naive_regression.png
 * \image latex naive_regression.eps
 * Both of the covariates are highly significant. And R^2 is 1, which tells us that our model fits perfectly.
 * Now we can take a look at the plot where we used the LinkedCells algorithm instead of the naive approach:
 * \image html lc_regression.png
 * \image latex lc_regression.eps
 * Both of the covariates are again significant. But now the quadratic term becomes less significant (p-value < 0.019). We can interpret this that this improved approach is "more linear". R^2 is again very high (0.99).
 * Here is a plot where we can observe that the LinkedCells algorithm is a way better than the first naive approach.
 * \image html both_regression.png
 * \image latex both_regression.eps
 * Now we can extrapolate our models.
 * \image html extrapolation.png
 * \image latex extrapolation.eps
 * For one iteration with 100000 Particles the LinkedCells algorithm is 133.4s faster than the naive approach.
 *
 */
class LinkedCells : public SimulationContainer {

private:
    /**
     * This vector stores cells. These cells store a vector of particles in this cell
     */
    std::vector<std::list<Particle>> particles;
    double reflectionDistance;

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
     * @param sigma sigma parameter of Lennard Jones Potential
     */
    explicit LinkedCells(std::array<int, 3> dimension, double mesh, double cutOff, double sigma);

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
    std::list<Particle> &get(std::array<double, 3> coords);

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
    std::list<Particle> &indexGet(int i);

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
    void forceInsert(Particle &p);

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
     * Iterates over particle in boundary cells and creates ghost particles for reflection if they are near enough to the boundary
     */
    void createGhosts();

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
     *
     * @return dimensions
     */
    std::array<int, 3> getDimensions();

};



