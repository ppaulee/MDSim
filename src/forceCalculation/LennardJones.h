//
// Created by Jonas on 13.11.21.
//

#ifndef PSEMOLDYN_GROUPB_LENNARDJONES_H
#define PSEMOLDYN_GROUPB_LENNARDJONES_H


#include "ForceCalculation.h"

/**
 * Uses Lennard Jones potential for force calculations
 */
class LennardJones : public ForceCalculation {
private:
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
};


#endif //PSEMOLDYN_GROUPB_LENNARDJONES_H
