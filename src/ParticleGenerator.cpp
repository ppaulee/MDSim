//
// Created by paulu on 10.11.2021.
//
#include "ParticleGenerator.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "vector"

void generateCube(std::array<double, 3> dimension, std::array<double, 3> startPoint, double h, double m, std::array<double, 3> v, double meanV, LinkedCells &container) {
    for (double x = 0; x < dimension[0]; x++) {
        for (double y = 0; y < dimension[1]; y++) {
            for (double z = 0; z < dimension[2]; z++) {
                // Add particle to ParticleContainer
                Particle p = Particle({x*h+startPoint[0],y*h+startPoint[1],z*h+startPoint[2]},v,m);
                container.insert(p);
            }
        }
    }
}

void generateFromFile(LinkedCells &particles, char *filename) {
    int num = 0;

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
        numstream >> num;
        std::cout << "Reading " << num << "." << std::endl;


        for (int i = 0; i < num*2; i++) {
            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;

            if (tmp_string == "Cube:\r" || tmp_string == "Cube:\n") {
                getline(input_file, tmp_string);
                parseCube(tmp_string, particles);
            }
        }
    } else {
        std::cout << "Error: could not open file " << filename << std::endl;
        exit(-1);
    }
}

void parseCube(std::string str, LinkedCells& particleContainer) {
    std::vector<std::string> strings = splitToString(str, ';');
    //Read input to arrays/doubles
    std::array<double, 3> startPoint = convertToFixedArray(splitToDouble(strings[0], ','));
    std::array<double, 3> v = convertToFixedArray(splitToDouble(strings[1], ','));
    std::array<double, 3> dimension = convertToFixedArray(splitToDouble(strings[2], ','));
    double h = std::stod(strings[3]);
    double meanV = std::stod(strings[4]);
    double mass = std::stod(strings[5]);

    generateCube(dimension, startPoint, h, mass, v, meanV, particleContainer);
}

std::vector<std::string> splitToString(std::string str, char delimiter) {
    std::vector<std::string> strings;
    std::istringstream f(str);
    std::string s;
    while(getline(f, s, delimiter)) {
        std::cout << s << std::endl;
        strings.push_back(s);
    }
    return strings;
}

std::vector<double> splitToDouble(std::string str, char delimiter) {
    std::vector<double> strings;
    std::istringstream f(str);
    std::string s;
    while(getline(f, s, delimiter)) {
        std::cout << s << std::endl;
        strings.push_back(std::stod(s));
    }
    return strings;
}

std::array<double,3> convertToFixedArray(std::vector<double> v) {
    return {v[0], v[1], v[2]};
}

