//
// Created by Jonas on 20.12.21.
//

#ifndef PSEMOLDYN_GROUPB_CHECKPOINTFILEWRITER_H
#define PSEMOLDYN_GROUPB_CHECKPOINTFILEWRITER_H

#include <container/SimulationContainer.h>

class CheckpointFileWriter{
public:
    CheckpointFileWriter();

    void writeFile(SimulationContainer &particles, const std::string &filename);
};

#endif //PSEMOLDYN_GROUPB_CHECKPOINTFILEWRITER_H
