//
// Created by Jonas on 13.11.21.
//

#include "utils/ArrayUtils.h"
#include "LennardJones.h"

std::array<double, 3> LennardJones::calculateF(Particle p1, Particle p2) {
    double normNoRoot = 0;
    std::array<double, 3> difference = p1.getX() - p2.getX();
    for (int i = 0; i < 3; i++) {
        normNoRoot += difference.at(i) * difference.at(i);
    }
    double powSigmaSix = pow(sigma, 6);
    double powNormNoRootThree = pow(normNoRoot, 3);
    double temp = (powSigmaSix * powNormNoRootThree) - (2 * (powSigmaSix * powSigmaSix));
    temp = (-24 * epsilon * temp) / (normNoRoot * powNormNoRootThree * powNormNoRootThree);
    std::array<double, 3> vec = temp * difference;
    return vec;
}

LennardJones::LennardJones(double epsilon, double sigma) : epsilon(epsilon), sigma(sigma) {}
