//
// Created by Jonas on 10.11.21.
//

#include <utils/ArrayUtils.h>
#include "StoermerVerlet.h"

std::array<double, 3> StoermerVerlet::calculateF(Particle p1, Particle p2) {
    //double f = (p1.getM() * p2.getM()) / pow(ArrayUtils::L2Norm(ArrayUtils::elementWisePairOp(p1.getX(), p2.getX(), sub)), 3);
    double f = (p1.getM() * p2.getM()) / pow(ArrayUtils::L2Norm(p1.getX() - p2.getX()), 3);
    //std::array<double, 3> vec = ArrayUtils::elementWisePairOp(p2.getX(),p1.getX(), sub);
    std::array<double, 3> vec = p2.getX() - p1.getX();
    //vec = ArrayUtils::elementWiseScalarOp(f, vec, multiply);
    vec = f * vec;
    return vec;
}
