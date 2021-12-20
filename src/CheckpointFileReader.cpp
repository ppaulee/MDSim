//
// Created by Jonas on 20.12.21.
//

#include "CheckpointFileReader.h"

#include "FileReader.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

CheckpointFileReader::CheckpointFileReader() = default;

CheckpointFileReader::~CheckpointFileReader() = default;

void CheckpointFileReader::readFile(SimulationContainer &particles, char *filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    std::array<double, 3> f;
    std::array<double, 3> oldF;
    double m;
    int type;
    double sigma;
    double epsilon;

    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

        getline(input_file, tmp_string);
        std::cout << "Read line: " << tmp_string << std::endl;

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_particles;
        std::cout << "Reading " << num_particles << "." << std::endl;
        getline(input_file, tmp_string);
        std::cout << "Read line: " << tmp_string << std::endl;

        for (int i = 0; i < num_particles; i++) {
            std::istringstream datastream(tmp_string);

            for (auto &xj : x) {
                datastream >> xj;
            }
            for (auto &vj : v) {
                datastream >> vj;
            }
            for (auto &fj : f) {
                datastream >> fj;
            }
            for (auto &oldfj : oldF) {
                datastream >> oldfj;
            }
            if (datastream.eof()) {
                std::cout
                        << "Error reading file: eof reached unexpectedly reading from line "
                        << i << std::endl;
                exit(-1);
            }
            datastream >> m;
            datastream >> type;
            datastream >> sigma;
            datastream >> epsilon;
            auto new_p = Particle(x, v, m, type, sigma, epsilon);
            new_p.setF(f);
            new_p.setOldF(oldF);
            particles.forceInsert(new_p);

            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;
        }
    } else {
        std::cout << "Error: could not open file " << filename << std::endl;
        exit(-1);
    }
}

