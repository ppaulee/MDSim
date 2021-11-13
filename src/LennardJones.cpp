//
// Created by Jonas on 13.11.21.
//

#include <utils/ArrayUtils.h>
#include "LennardJones.h"

std::array<double, 3> LennardJones::calculateF(Particle p1, Particle p2) {
    double norm = ArrayUtils::L2Norm(p1.getX() - p2.getX());
    double temp = pow(sigma / norm, 6);
    temp = temp - 2 * pow(sigma / norm, 12);
    temp = temp * (24 * epsilon) / (norm * norm);
    std::array<double, 3> vec = temp * (p1.getX() - p2.getX());
    return vec;
}

LennardJones::LennardJones(double epsilon, double sigma) : epsilon(epsilon), sigma(sigma) {}
