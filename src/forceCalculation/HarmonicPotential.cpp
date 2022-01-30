//
// Created by paulu on 22.01.2022.
//

#include <utils/ArrayUtils.h>
#include "HarmonicPotential.h"
#include <iostream>

std::array<double, 3> HarmonicPotential::calculateF(const Particle &p1, const Particle &p2) {
    double l2norm = ArrayUtils::L2Norm(p1.getX() - p2.getX());
    if (std::isinf(l2norm)) {
        std::cout << "NAN in calculateF early" << std::endl;
        return {0,0,0};
    }
    double tmp = k * (l2norm - r0);
    std::array<double, 3> tmp_t = tmp * ( (1 / l2norm) * (p2.getX() - p1.getX()) );
    if (std::isnan(tmp_t[0]) || std::isnan(tmp_t[1]) || std::isnan(tmp_t[2])) {
        std::cout << "NAN in calculateF" << std::endl;
        return {0,0,0};
    }
    return tmp * ( (1 / l2norm) * (p2.getX() - p1.getX()) );
}

std::array<double, 3> HarmonicPotential::calculateFDiag(const Particle &p1, const Particle &p2) {
    double l2norm = ArrayUtils::L2Norm(p1.getX() - p2.getX());
    if (std::isinf(l2norm)) {
        std::cout << "NAN in calculateFDiag early" << std::endl;
        return {0,0,0};
    }
    double tmp = k * (l2norm - sqrt(2) * r0);
    std::array<double, 3> tmp_t = tmp * ( (1 / l2norm) * (p2.getX() - p1.getX()) );
    if (std::isnan(tmp_t[0]) || std::isnan(tmp_t[1]) || std::isnan(tmp_t[2])) {
        std::cout << "NAN in calculateFDiag" << std::endl;
        return {0,0,0};
    }
    return tmp * ( (1 / l2norm) * (p2.getX() - p1.getX()) );
}

HarmonicPotential::HarmonicPotential(double k, double r0) : k(k), r0(r0) {};

