//
// Created by Jonas on 20.12.21.
//

#include "CheckpointFileWriter.h"
#include <iomanip>
#include <sstream>
#include <fstream>

void CheckpointFileWriter::writeFile(SimulationContainer &particles, const std::string &filename) {

    std::ofstream file;
    std::stringstream strstr;
    strstr << filename << ".txt";

    file.open(strstr.str().c_str());
    file << particles.numberParticles() << std::endl;

    std::vector<Particle> pars = particles.getParticles();
    for (auto &p: pars) {
        for (auto x: p.getX()) {
            file << x << " ";
        }
        for (auto v: p.getV()) {
            file << v << " ";
        }
        for (auto f: p.getF()) {
            file << f << " ";
        }
        for (auto oldF: p.getOldF()) {
            file << oldF << " ";
        }
        file << p.getM() << " " << p.getType() << " " << p.getSigma() << " " << p.getEpsilon() << std::endl;
    }
}

CheckpointFileWriter::CheckpointFileWriter() = default;
