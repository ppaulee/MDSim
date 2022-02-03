//
// Created by ethan on 2/1/2022.
//

#ifndef PSEMOLDYN_GROUPB_DENSITYANDVELOCITYTOCSVWRITER_H
#define PSEMOLDYN_GROUPB_DENSITYANDVELOCITYTOCSVWRITER_H

#include "DensityAndVelocityProfile.h"
#include <string>

class DensityAndVelocityToCSVWriter {
private:
    std::string fileName;
    int linesCount;

public:
    explicit DensityAndVelocityToCSVWriter(std::string s);

    void writeToCSV(DensityAndVelocityProfile &profile, int iteration);

    double averageValue(std::list<double> list);
};



#endif //PSEMOLDYN_GROUPB_DENSITYANDVELOCITYTOCSVWRITER_H
