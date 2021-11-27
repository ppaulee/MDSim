//
// Created by paulu on 21.11.2021.
//

#include <vector>
#include "Particle.h"
#include "ForceCalculation.h"


class LinkedCells {

private:
    /**
     * This vector stores cells. These cells store a vector of particles in this cell
     */
    std::vector<std::vector<Particle>> particles;

    /**
     * Width, Height, Depth
     */
    std::array<int, 3> dimensions;
    int size;
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


public:

    /**
     *
     * @param dimension Dimensions (x,y,z) of the inner cells
     * @param mesh Mesh size of the grid
     * @param cutOff cut off radius
     */
    explicit LinkedCells(std::array<int, 3> dimension, double mesh, double cutOff);

    /**
     *
     * @param coords Coordinates of the cell
     * @return Returns the content of a cell
     */
    std::vector<Particle>& get(std::array<double, 3> coords);

    /**
     * Inserts a particle into the data structure
     *
     * @param p Particle to insert
     */
    void insert(Particle& p);

    /**
     * Removes a particle from the data structure
     *
     * @param p Particle to remove
     */
    void remove(Particle& p);

    void test();

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
     * Calculates new velocity
     */
    void calculateV();

    /**
     * Calculates new position
     */
    void calculateX();

};


