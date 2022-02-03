//
// Created by ethan on 2/1/2022.
//

#include "DensityAndVelocityToCSVWriter.h"
#include <iomanip>
#include <sstream>
#include <fstream>


DensityAndVelocityToCSVWriter::DensityAndVelocityToCSVWriter(std::string filename){
    fileName = filename;
    linesCount = 0;
};

void DensityAndVelocityToCSVWriter::writeToCSV(DensityAndVelocityProfile &profile, int iteration) {
    std::fstream file;
    file.open(fileName, std::ios::out | (linesCount ? std::ios::app : std::ios::trunc));

    file << "Iteration ," << iteration << "\n";
    linesCount++;
    file << "density\n";
    linesCount++;
    file << "bin Number:";
    for(int  i = 1; i<=profile.getDensities().size(); i++){
        file << "," << i;
    }
    file << "\n";
    linesCount++;
    for(int  i = 0; i<profile.getDensities().size(); i++){
        file << "," << averageValue(profile.getDensities()[i]);
    }
    file << "\n\n";
    linesCount++;
    file << "velocity\n";
    linesCount++;
    file << "bin Number:";
    for(int  i = 1; i<=profile.getVelocities().size(); i++){
        file << "," << i;
    }
    file << "\n";
    linesCount++;
    for(int  i = 0; i<profile.getVelocities().size(); i++){
        file << "," << averageValue(profile.getVelocities()[i]);
    }
    file << "\n";
    linesCount++;
    file.close();
}

double DensityAndVelocityToCSVWriter::averageValue(std::list<double> list) {
    return std::accumulate(list.begin(), list.end(), 0.0) / list.size();
}