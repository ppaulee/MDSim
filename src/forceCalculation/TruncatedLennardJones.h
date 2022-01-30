//
// Created by paulu on 30.01.2022.
//
#pragma once

#include "LennardJones.h"


class TruncatedLennardJones : public LennardJones {

private:
    /**
     *  Parameter for Lennard Jones potential
     */
    double epsilon;
    /**
     *  Parameter for Lennard Jones potential
     */
    double sigma;
public:
    std::array<double, 3> calculateF(Particle p1, Particle p2);

    TruncatedLennardJones(double epsilon1, double sigma1, int type, double epsilon,
                          double sigma);
};

