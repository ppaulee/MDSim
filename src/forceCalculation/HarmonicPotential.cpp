//
// Created by paulu on 22.01.2022.
//

#include <utils/ArrayUtils.h>
#include "HarmonicPotential.h"

std::array<double, 3> HarmonicPotential::calculateF(Particle p1, Particle p2) {
    double l2norm = ArrayUtils::L2Norm(p1.getX() - p2.getX());
    double tmp = k * (l2norm - r0);
    return tmp * ( (1 / l2norm) * (p2.getX() - p1.getX()) );
}

std::array<double, 3> HarmonicPotential::calculateFDiag(Particle p1, Particle p2) {
    double l2norm = ArrayUtils::L2Norm(p1.getX() - p2.getX());
    double tmp = k * (l2norm - sqrt(2) * r0);
    return tmp * ( (1 / l2norm) * (p2.getX() - p1.getX()) );
}

HarmonicPotential::HarmonicPotential(double k, double r0) : k(k), r0(r0) {};

