//
// Created by Jonas on 10.11.21.
//

#ifndef PSEMOLDYN_GROUPB_FORCECALCULATION_H
#define PSEMOLDYN_GROUPB_FORCECALCULATION_H


#include <array>
#include "Particle.h"

class ForceCalculation {
public:
    /**
     * Calculates force between 2 particles
     *
     * @return new force
     */
    virtual std::array<double, 3> calculateF(Particle p1, Particle p2) = 0;
};


#endif //PSEMOLDYN_GROUPB_FORCECALCULATION_H
