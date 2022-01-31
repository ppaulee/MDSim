//
// Created by paulu on 22.01.2022.
//
#pragma once

#include "ForceCalculation.h"

class HarmonicPotential {
public:
    // stiffness constant
    double k;
    // average bond length
    double r0;
    explicit HarmonicPotential(double k, double r0);

    std::array<double, 3> calculateF(const Particle p1, const Particle p2);

    std::array<double, 3> calculateFDiag(const Particle p1, const Particle p2);
};
