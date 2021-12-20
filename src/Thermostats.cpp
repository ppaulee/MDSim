//
// Created by paulu on 12.12.2021.
//

#include "Thermostats.h"

Thermostats::Thermostats(SimulationContainer& container, double delta_t, double end_time, double initialTemp, int stepSize, int dimensions, double targetTemp, double max_delta_temp) {
    currentTemp = initialTemp;
    maxDeltaTemp = max_delta_temp;
    step = stepSize;
    currentStep = 0;
    holdTemp = false;
    if (targetTemp == -1) {
        holdTemp = true;
        targetTemperature = currentTemp;
    } else {
        targetTemperature = targetTemp;
    }
    dim = dimensions;
    endTime = end_time;
    deltaT = delta_t;

    // Set temp to initial value
    calcCurrentTemperature(container);
    double scale = sqrt(initialTemp / currentTemp);
    container.scaleVelocity(scale);
}

void Thermostats::calcCurrentTemperature(SimulationContainer& container) {
    double energy = container.calcKineticEnergy();
    currentTemp = ((double )(2*energy)) / ((double) (container.numberParticles() * dim));
}

void Thermostats::adjustTemperature(SimulationContainer& container, double current_time) {
    if (currentStep % step == 0) {
        // Hold the current temperature
        if (holdTemp) {
            double scale = sqrt(targetTemperature / currentTemp);
            // Cap by maxDeltaTemp
            if (maxDeltaTemp == -1) {
                // Calculate difference of current Temp and target Temp
                double difference = targetTemperature- currentTemp;
                // Cap the difference to maxDeltaTemp
                double difference_capped = std::min(std::abs(difference), maxDeltaTemp);
                // Add capped temperature difference to currentTemp (note that we keep the sign)
                double target_temp = currentTemp + (sign(difference) * difference_capped);
                scale = sqrt(target_temp / currentTemp);
            }
            container.scaleVelocity(scale);
        } else {
            // Calculate the overall number of steps the simulation needs
            double numberOfSteps = endTime / deltaT;
            // Calculate the number of steps which are left
            double numberOfStepsLeft = numberOfSteps - (current_time / deltaT);
            // Calculate the difference of the temperature
            double difference = targetTemperature- currentTemp;
            // Cap the difference to maxDeltaTemp
            double difference_capped = difference;
            if (maxDeltaTemp != -1) {
                difference_capped = sign(difference) * std::min(std::abs(difference), maxDeltaTemp);
            }
            // Calculate the increase of temperature for one timestep in the current system
            double tempPerStep = (difference_capped) / numberOfStepsLeft;

            // Scale the tempPerStep with the stepSize
            double temp = step * tempPerStep;
            // If the stepSize is less than the remaining steps then we use this to scale
            if (step > numberOfStepsLeft) {
                temp = numberOfStepsLeft * tempPerStep;
            }
            // Cap again to maxDeltaTemp
            if (maxDeltaTemp != -1) {
                temp = std::min(temp, maxDeltaTemp);
            }
            container.scaleVelocity(sqrt((temp + currentTemp) / currentTemp));
        }
    }
    currentStep++;
}

int Thermostats::sign(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

double Thermostats::getCurrentTemperature() {
    return currentTemp;
}