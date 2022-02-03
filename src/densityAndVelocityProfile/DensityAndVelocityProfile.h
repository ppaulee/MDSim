//
// Created by ethan on 2/1/2022.
//

#pragma once

#include "./container/SimulationContainer.h"
#include <vector>
#include <list>


class DensityAndVelocityProfile {
private:
    /**
     *
     */
     int binSize;

    /**
     * width for density calculation
     */
     double width;

    /**
     * h of liquid
     */
     double h;

    /**
     * starting Point of the liquid particles for the bins
     */
     std::array<double, 3> liquidStartPoint;

    /**
     * Bins for densities
     */
     std::vector<std::list<double>> densities;

     /**
      * Bins for velocities
      */
     std::vector<std::list<double>> velocities;

public:
    /**
     *
     * @param binSize size for the densities and velocities bins according the x axis
     * @param width
     */
    explicit DensityAndVelocityProfile(int binSize, double width, double h, std::array<double, 3> liquidStartPoint);

    /**
     * Calculates densities and velocities with given data structure (particles)
     *
     * @param container Data structure with particles
     */
    void calculateDnV(SimulationContainer& container);

    /**
     * Calculates densities with given data structure (particles)
     *
     * @param container Data structure with particles
     */
    void calculateDensities(SimulationContainer& container);

    /**
     * Calculates velocities with given data structure (particles)
     *
     * @param container Data structure with particles
     */
    void calculateVelocities(SimulationContainer& container);

    std::vector<std::list<double>> getDensities();

    std::vector<std::list<double>> getVelocities();
};



