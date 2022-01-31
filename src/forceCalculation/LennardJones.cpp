//
// Created by Jonas on 13.11.21.
//

#include "utils/ArrayUtils.h"
#include "LennardJones.h"

std::array<double, 3> LennardJones::calculateF(Particle p1, Particle p2) {
    if (isMembrane_bool) {
        if (ArrayUtils::L2Norm(p1.getX() - p2.getX()) <= 1.1225) {
            return {0,0,0};
        } else {
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
    } else {
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
}


LennardJones::LennardJones(double epsilon, double sigma, int type) : epsilon(epsilon), sigma(sigma), type(type) {
    isMembrane_bool = false;
}

double LennardJones::getEpsilon() const {
    return epsilon;
}

double LennardJones::getSigma() const {
    return sigma;
}

int LennardJones::getType() const {
    return type;
}

void LennardJones::setMembrane() {
    isMembrane_bool = true;
}

bool LennardJones::isMembrane() {
    return isMembrane_bool;
}

LennardJones::LennardJones(double epsilon, double sigma) : epsilon(epsilon), sigma(sigma) {
    isMembrane_bool = false;
}
