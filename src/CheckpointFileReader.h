//
// Created by Jonas on 20.12.21.
//

#ifndef PSEMOLDYN_GROUPB_CHECKPOINTFILEREADER_H
#define PSEMOLDYN_GROUPB_CHECKPOINTFILEREADER_H

#include "Particle.h"
#include "container/SimulationContainer.h"
#include <vector>

class CheckpointFileReader {

public:
    CheckpointFileReader();

    virtual ~CheckpointFileReader();

    void readFile(SimulationContainer &particles, char *filename);
};

#endif //PSEMOLDYN_GROUPB_CHECKPOINTFILEREADER_H
