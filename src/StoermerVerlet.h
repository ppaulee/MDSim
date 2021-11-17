//
// Created by schnitzel on 10.11.21.
//

#ifndef PSEMOLDYN_GROUPB_STOERMERVERLET_H
#define PSEMOLDYN_GROUPB_STOERMERVERLET_H

#include "ForceCalculation.h"

/**
 * Uses Stoermer-Verlet algorithm for force calculation
 */
class StoermerVerlet : public ForceCalculation {
public:
    StoermerVerlet() = default;

    std::array<double, 3> calculateF(Particle p1, Particle p2);
};

#endif //PSEMOLDYN_GROUPB_STOERMERVERLET_H
