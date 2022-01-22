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
     * Simulates a membrane
     *
     * @param delta_t Timestep size
     */
    virtual void simulateMembrane(double delta_t, bool pullState) = 0;

    /**
     * Inserts a particles into the datastructure
     *
     * @param p
     */
    virtual void insert(Particle& p) = 0;

    virtual void forceInsert(Particle &p) = 0;

    /**
     * Adds a brownian motion to all particles
     *
     * @param averageV
     * @param dimension
     */
    virtual void addBrownianMotion(double averageV, int dimension) = 0;

    /**
     * Calculates the kinetic energy in the system
     *
     * @return Kinetic energy in the whole system
     */
    virtual double calcKineticEnergy() = 0;

    /**
     *
     * @return Number of particles
     */
    virtual int numberParticles() = 0;

    /**
     * Scale all velocities by a specific value
     *
     * @param scale Value to scale
     */
    virtual void scaleVelocity(double scale) = 0;

    /**
     * Get all particles in the data structure
     * @return vector of particles
     */
    virtual std::vector<Particle> getParticles() = 0;
};

