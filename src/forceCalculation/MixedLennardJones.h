//
// Created by Jonas on 18.12.21.
//

#ifndef PSEMOLDYN_GROUPB_MIXEDLENNARDJONES_H
#define PSEMOLDYN_GROUPB_MIXEDLENNARDJONES_H


#include "ForceCalculation.h"

/**
 * Uses Lennard Jones potential for force calculations
 * This variant also stores between which particle types the parameters are valid
 * The parameters epsilon and sigma should be set according to mixing rules
 */
class MixedLennardJones : public ForceCalculation {
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
     *  This instance of Mixed Lennard Jones is used for calculations between 2 particles of these 2 different type
     *  This is the type of the first particle
     */
    int type1;

    /**
     *  This instance of Mixed Lennard Jones is used for calculations between 2 particles of these 2 different type
     *  This is the type of the second particle
     */
    int type2;

public:

    MixedLennardJones(double epsilonArg, double sigmaArg, int type1Arg, int type2Arg);

    double getEpsilon() const;

    double getSigma() const;

    int getType1() const;

    int getType2() const;

    std::array<double, 3> calculateF(Particle p1, Particle p2);
};


#endif //PSEMOLDYN_GROUPB_MIXEDLENNARDJONES_H
