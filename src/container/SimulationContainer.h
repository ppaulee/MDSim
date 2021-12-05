//
// Created by paulu on 05.12.2021.
//
#pragma once

#include "utils/MaxwellBoltzmannDistribution.h"
#include "forceCalculation/ForceCalculation.h"

class SimulationContainer {

public:

    /**
     * Plots particles
     *
     * @param iteration current Iteration
     */
    virtual void plotParticles(int iteration) = 0;


    /**
     * Simulates the particles
     *
     * @param algorithm Force calculation algorithm used
     * @param delta_t Timestep size
     */
    virtual void simulate(ForceCalculation *algorithm, double delta_t) = 0;

    /**
     * Inserts a particles into the datastructure
     *
     * @param p
     */
    virtual void insert(Particle& p) = 0;

    /**
     * Adds a brownian motion to all particles
     *
     * @param averageV
     * @param dimension
     */
    virtual void addBrownianMotion(double averageV, int dimension) = 0;
};

