//
// Created by schnitzel on 10.11.21.
//

#ifndef PSEMOLDYN_GROUPB_GRAVITATION_H
#define PSEMOLDYN_GROUPB_GRAVITATION_H

#include "ForceCalculation.h"

/**
 * Uses Stoermer-Verlet algorithm for force calculation
 */
class Gravitation : public ForceCalculation {
public:
    Gravitation() = default;

    std::array<double, 3> calculateF(Particle p1, Particle p2);
};

#endif //PSEMOLDYN_GROUPB_GRAVITATION_H
