//
// Created by paulu on 30.01.2022.
//

#include "TruncatedLennardJones.h"
#include "utils/ArrayUtils.h"

std::array<double, 3> TruncatedLennardJones::calculateF(Particle p1, Particle p2) {
    if (ArrayUtils::L2Norm(p1.getX() - p2.getX()) >= 1.1225*sigma) {
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
}

TruncatedLennardJones::TruncatedLennardJones(double epsilon1, double sigma1, int type, double epsilon,
                                             double sigma) : LennardJones(epsilon1, sigma1, type), epsilon(epsilon),
                                                             sigma(sigma) {}