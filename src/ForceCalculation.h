//
// Created by Jonas on 10.11.21.
//

#ifndef PSEMOLDYN_GROUPB_FORCECALCULATION_H
#define PSEMOLDYN_GROUPB_FORCECALCULATION_H


#include <array>
#include "Particle.h"

/**
 * Abstract class whose children represent different types of algorithms for force calculation
 */
class ForceCalculation {
public:
    /**
     * Calculates force between 2 particles using the classes algorithm
     *
     * @param p1 First Particle
     * @param p2 Second Particle
     * @return new force
     */
    virtual std::array<double, 3> calculateF(Particle p1, Particle p2) = 0;
};


#endif //PSEMOLDYN_GROUPB_FORCECALCULATION_H
