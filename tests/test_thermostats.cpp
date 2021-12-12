//
// Created by paulu on 12.12.2021.
//

#include "gtest/gtest.h"
#include "../src/Thermostats.h"
#include "../src/Thermostats.cpp"

/**
 * Warms up the simulation and compare the end temperature with the expected temperature
 */
TEST(Thermostats, warm_up) {
    auto cells = new LinkedCells({40, 40, 0}, 1, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 1, 0}, 1, 1);
    auto p2 = new Particle({1, 1, 0}, {1, 1, 0}, 1, 1);

    cells->insert(*p1);
    cells->insert(*p2);

    double delta_t = 1;
    double end_t = 5;
    int initialTemp = 20;
    int stepSize = 1;
    int dimensions = 2;
    int targetTemp = 40;
    double epsilon = 5;
    double sigma = 1;

    double current_time = 0;
    auto algorithm = new LennardJones(epsilon, sigma);

    auto thermo = new Thermostats(*cells, delta_t, end_t, initialTemp, stepSize, dimensions, targetTemp);
    thermo->calcCurrentTemperature(*cells);
    thermo->adjustTemperature(*cells, current_time);
    while (current_time < end_t) {
        cells->simulate(algorithm, delta_t);
        thermo->calcCurrentTemperature(*cells);
        thermo->adjustTemperature(*cells, current_time);
        thermo->calcCurrentTemperature(*cells);
        current_time += delta_t;
    }
    EXPECT_NEAR(thermo->getCurrentTemperature(), 40, 0.001);
}

/**
 * Warms up the simulation and compare the end temperature with the expected temperature. Here we use a step size of 2 and a end_time of 6. Here
 * we have a multiple of 2 as the end_time
 */
TEST(Thermostats, warm_up_steps) {
    auto cells = new LinkedCells({40, 40, 0}, 1, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 1, 0}, 1, 1);
    auto p2 = new Particle({1, 1, 0}, {1, 1, 0}, 1, 1);

    cells->insert(*p1);
    cells->insert(*p2);

    double delta_t = 1;
    double end_t = 6;
    int initialTemp = 20;
    int stepSize = 2;
    int dimensions = 2;
    int targetTemp = 40;
    double epsilon = 5;
    double sigma = 1;

    double current_time = 0;
    auto algorithm = new LennardJones(epsilon, sigma);

    auto thermo = new Thermostats(*cells, delta_t, end_t, initialTemp, stepSize, dimensions, targetTemp);
    thermo->calcCurrentTemperature(*cells);
    thermo->adjustTemperature(*cells, current_time);
    while (current_time < end_t) {
        cells->simulate(algorithm, delta_t);
        thermo->calcCurrentTemperature(*cells);
        thermo->adjustTemperature(*cells, current_time);
        thermo->calcCurrentTemperature(*cells);
        current_time += delta_t;
    }
    EXPECT_NEAR(thermo->getCurrentTemperature(), 40, 0.001);
}

/**
 * Warms up the simulation and compare the end temperature with the expected temperature. Here we use a step size of 2 and a end_time of 6. Here
 * we do not have a multiple of 2 as the end_time
 */
TEST(Thermostats, warm_up_steps_off) {
    auto cells = new LinkedCells({40, 40, 0}, 1, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 1, 0}, 1, 1);
    auto p2 = new Particle({1, 1, 0}, {1, 1, 0}, 1, 1);

    cells->insert(*p1);
    cells->insert(*p2);

    double delta_t = 1;
    double end_t = 5;
    int initialTemp = 20;
    int stepSize = 2;
    int dimensions = 2;
    int targetTemp = 40;
    double epsilon = 5;
    double sigma = 1;

    double current_time = 0;
    auto algorithm = new LennardJones(epsilon, sigma);

    auto thermo = new Thermostats(*cells, delta_t, end_t, initialTemp, stepSize, dimensions, targetTemp);
    thermo->calcCurrentTemperature(*cells);
    thermo->adjustTemperature(*cells, current_time);
    while (current_time < end_t) {
        cells->simulate(algorithm, delta_t);
        thermo->calcCurrentTemperature(*cells);
        thermo->adjustTemperature(*cells, current_time);
        thermo->calcCurrentTemperature(*cells);
        current_time += delta_t;
    }
    EXPECT_NEAR(thermo->getCurrentTemperature(), 40, 0.001);
}

/**
 * Cools down the simulation and compare the end temperature with the expected temperature
 */
TEST(Thermostats, cool_down) {
    auto cells = new LinkedCells({40, 40, 0}, 1, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 1, 0}, 1, 1);
    auto p2 = new Particle({1, 1, 0}, {1, 1, 0}, 1, 1);

    cells->insert(*p1);
    cells->insert(*p2);

    double delta_t = 1;
    double end_t = 5;
    int initialTemp = 40;
    int stepSize = 1;
    int dimensions = 2;
    int targetTemp = 20;
    double epsilon = 5;
    double sigma = 1;

    double current_time = 0;
    auto algorithm = new LennardJones(epsilon, sigma);

    auto thermo = new Thermostats(*cells, delta_t, end_t, initialTemp, stepSize, dimensions, targetTemp);
    thermo->calcCurrentTemperature(*cells);
    thermo->adjustTemperature(*cells, current_time);
    while (current_time < end_t) {
        cells->simulate(algorithm, delta_t);
        thermo->calcCurrentTemperature(*cells);
        thermo->adjustTemperature(*cells, current_time);
        thermo->calcCurrentTemperature(*cells);
        current_time += delta_t;
    }
    EXPECT_NEAR(thermo->getCurrentTemperature(), 20, 0.001);
}

/**
 * Cools down the simulation and compare the end temperature with the expected temperature. Here we use a step size of 2 and a end_time of 6. Here
 * we have a multiple of 2 as the end_time
 */
TEST(Thermostats, cool_down_steps) {
    auto cells = new LinkedCells({40, 40, 0}, 1, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 1, 0}, 1, 1);
    auto p2 = new Particle({1, 1, 0}, {1, 1, 0}, 1, 1);

    cells->insert(*p1);
    cells->insert(*p2);

    double delta_t = 1;
    double end_t = 6;
    int initialTemp = 40;
    int stepSize = 2;
    int dimensions = 2;
    int targetTemp = 20;
    double epsilon = 5;
    double sigma = 1;

    double current_time = 0;
    auto algorithm = new LennardJones(epsilon, sigma);

    auto thermo = new Thermostats(*cells, delta_t, end_t, initialTemp, stepSize, dimensions, targetTemp);
    thermo->calcCurrentTemperature(*cells);
    thermo->adjustTemperature(*cells, current_time);
    while (current_time < end_t) {
        cells->simulate(algorithm, delta_t);
        thermo->calcCurrentTemperature(*cells);
        thermo->adjustTemperature(*cells, current_time);
        thermo->calcCurrentTemperature(*cells);
        current_time += delta_t;
    }
    EXPECT_NEAR(thermo->getCurrentTemperature(), 20, 0.001);
}

/**
 * Hold the temperature of the simulation and compares the end temperature with the initial temperature
 */
TEST(Thermostats, hold) {
    auto cells = new LinkedCells({40, 40, 0}, 1, 1, 5);
    auto p1 = new Particle({0, 0, 0}, {0, 1, 0}, 1, 1);
    auto p2 = new Particle({1, 1, 0}, {1, 1, 0}, 1, 1);

    cells->insert(*p1);
    cells->insert(*p2);

    double delta_t = 1;
    double end_t = 5;
    int initialTemp = 40;
    int stepSize = 1;
    int dimensions = 2;
    int targetTemp = 40;
    double epsilon = 5;
    double sigma = 1;

    double current_time = 0;
    auto algorithm = new LennardJones(epsilon, sigma);

    auto thermo = new Thermostats(*cells, delta_t, end_t, initialTemp, stepSize, dimensions, targetTemp);
    thermo->calcCurrentTemperature(*cells);
    thermo->adjustTemperature(*cells, current_time);
    while (current_time < end_t) {
        cells->simulate(algorithm, delta_t);
        thermo->calcCurrentTemperature(*cells);
        thermo->adjustTemperature(*cells, current_time);
        thermo->calcCurrentTemperature(*cells);
        current_time += delta_t;
    }
    EXPECT_NEAR(thermo->getCurrentTemperature(), 40, 0.001);
}