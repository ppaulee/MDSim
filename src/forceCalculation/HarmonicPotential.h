//
// Created by paulu on 22.01.2022.
//
#pragma once

#include "ForceCalculation.h"

class HarmonicPotential : public ForceCalculation {
public:
    // stiffness constant
    double k;
    // average bond length
    double r0;
    explicit HarmonicPotential(double k, double r0);

    std::array<double, 3> calculateF(Particle p1, Particle p2) override;

    std::array<double, 3> calculateFDiag(Particle p1, Particle p2);
};
