//
// Created by paulu on 12.12.2021.
//

#pragma once
#include "./container/SimulationContainer.h"

class Thermostats {
private:
    /**
     * current Temperature in Kelvin
     */
    double currentTemp;

    /**
     * Target temperature in Kelvin
     */
    double targetTemperature;

    int step;
    int currentStep;

    double maxDeltaTemp;

    int dim;

    bool holdTemp;

    double deltaT;
    double endTime;

    /**
     * Calculates the sign of the number
     *
     * @param x
     * @return -1 for <0, 0 for ==0, 1 for >0
     */
    int sign(double x);

public:
    /**
     *
     * @param delta_t Timesteps per iteration
     * @param end_time End time of the simulation
     * @param initialTemp Initial Temperature
     * @param stepSize The number of time steps after which the thermostat is periodically applied.
     * @param targetTemp Target temperature
     * @param max_delta_temp the maximal absolute value allowed for the temperature should to be changed in one scaling step.
     * @param dimensions dimensions of the simulation
     */
    explicit Thermostats(SimulationContainer& container, double delta_t, double end_time, double initialTemp, int stepSize, int dimensions, double targetTemp = -1, double max_delta_temp = -1);

    /**
     * Calculates current Temperature in the system and sets currentTemperature
     *
     * @param container Data structure with particles
     */
    void calcCurrentTemperature(SimulationContainer& container);

    /**
     * Control the temperature
     *
     * @param container Data structure with particles
     * @param current_time Current time in simulation
     */
    void adjustTemperature(SimulationContainer& container, double current_time);

    /**
     *
     * @return current temperature
     */
    double getCurrentTemperature();



};



