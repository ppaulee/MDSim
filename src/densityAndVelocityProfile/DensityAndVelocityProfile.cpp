//
// Created by ethan on 2/1/2022.
//

#include "DensityAndVelocityProfile.h"
#include <cmath>
#include "../utils/ArrayUtils.h"

DensityAndVelocityProfile::DensityAndVelocityProfile(int size, double w, double hh, std::array<double, 3> startPoint) {
    binSize = size;
    width = w;
    h = hh;
    liquidStartPoint = startPoint;

    for (int i = 0; i < binSize; i++) {
        densities.push_back(std::list<double>({}));
    }

    for (int i = 0; i < binSize; i++) {
        velocities.push_back(std::list<double>({}));
    }

}

void DensityAndVelocityProfile::calculateDnV(SimulationContainer &container) {
    for(auto &d : densities){
        d.clear();
    }
    for(auto &v : velocities){
        v.clear();
    }
    calculateDensities(container);
    calculateVelocities(container);
}

void DensityAndVelocityProfile::calculateDensities(SimulationContainer &container) {
    for(int i = 0; i<densities.size();i++){
        double innerBound = 1.0 * i;
        double outerBound = innerBound + width;
        for(auto &current_p : container.getParticles()){
            int sizeN = 0;
            for(auto &p : container.getParticles()){
                if (((ArrayUtils::L2Norm(p.getX() - current_p.getX())) <=  (outerBound-innerBound)) && !(p == current_p)){
                    sizeN++;
                }
            }
            double density = sizeN*3 / (4*M_PI*(outerBound*outerBound*outerBound-innerBound*innerBound*innerBound));
            int index = (int) ((current_p.getX()[0]/h) -liquidStartPoint[0]);
            if(index<0){
                index = 0;
            }
            if(index>=velocities.size()){
                index = velocities.size()-1;
            }
            densities[index].emplace_back(density);
        }
    }
}

void DensityAndVelocityProfile::calculateVelocities(SimulationContainer &container) {
    std::vector<Particle> pars = container.getParticles();
    for(auto &current_p : pars){
        double v = ArrayUtils::L2Norm(current_p.getV());
        int index = (int) ((current_p.getX()[0]/h) -liquidStartPoint[0]);
        if(index<0){
            index = 0;
        }
        if(index>=velocities.size()){
            index = velocities.size()-1;
        }
        velocities[index].emplace_back(v);
    }
}

std::vector<std::list<double>> DensityAndVelocityProfile::getDensities() {
    return densities;
}

std::vector<std::list<double>> DensityAndVelocityProfile::getVelocities() {
    return velocities;
}