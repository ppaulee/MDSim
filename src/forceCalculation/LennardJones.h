//
// Created by Jonas on 13.11.21.
//
#pragma once


#include "ForceCalculation.h"

/**
 * Uses Lennard Jones potential for force calculations
 */
class LennardJones : public ForceCalculation {
private:

    bool isMembrane_bool;

    bool isNano_bool;

    /**
     *  Parameter for Lennard Jones potential
     */
    double epsilon;
    /**
     *  Parameter for Lennard Jones potential
     */
    double sigma;
    /**
     *  This instance of Lennard Jones is used for calculations between 2 particles of this, same type
     */
     int type;
public:

    LennardJones(double epsilon, double sigma, int type);

    LennardJones(double epsilon, double sigma);

    double getEpsilon() const;

    double getSigma() const;

    int getType() const;

    std::array<double, 3> calculateF(Particle p1, Particle p2);

    void setMembrane();

    bool isMembrane();

    void setNano();

    bool isNano();
};

