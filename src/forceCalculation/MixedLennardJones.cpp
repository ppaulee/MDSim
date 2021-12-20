//
// Created by Jonas on 18.12.21.
//

#include "MixedLennardJones.h"
#include "utils/ArrayUtils.h"

std::array<double, 3> MixedLennardJones::calculateF(Particle p1, Particle p2) {
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

double MixedLennardJones::getEpsilon() const {
    return epsilon;
}

double MixedLennardJones::getSigma() const {
    return sigma;
}

int MixedLennardJones::getType1() const {
    return type1;
}

int MixedLennardJones::getType2() const {
    return type2;
}

MixedLennardJones::MixedLennardJones(double epsilonArg, double sigmaArg, int type1Arg, int type2Arg){
    epsilon = epsilonArg;
    sigma = sigmaArg;
    type1 = type1Arg;
    type2 = type2Arg;
}
