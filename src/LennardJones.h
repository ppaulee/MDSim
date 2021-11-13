//
// Created by Jonas on 13.11.21.
//

#ifndef PSEMOLDYN_GROUPB_LENNARDJONES_H
#define PSEMOLDYN_GROUPB_LENNARDJONES_H


#include "ForceCalculation.h"

class LennardJones : public ForceCalculation {
private:
    double epsilon;
    double sigma;
public:

    LennardJones(double epsilon, double sigma);

    std::array<double, 3> calculateF(Particle p1, Particle p2);
};


#endif //PSEMOLDYN_GROUPB_LENNARDJONES_H
