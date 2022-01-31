//
// Created by paulu on 10.11.2021.
//
#include "ParticleGenerator.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <list>
#include "vector"
#include "utils/ArrayUtils.h"


void generateCube(std::array<int, 3> dimension, std::array<double, 3> startPoint, double h, double m,
                  std::array<double, 3> v, double meanV, SimulationContainer &container) {
    for (double x = 0; x < dimension[0]; x++) {
        for (double y = 0; y < dimension[1]; y++) {
            for (double z = 0; z < dimension[2]; z++) {
                // Add particle to ParticleContainer
                Particle p = Particle({x * h + startPoint[0], y * h + startPoint[1], z * h + startPoint[2]}, v, m);
                container.insert(p);
            }
        }
    }
}

void generateCube(std::array<int, 3> dimension, std::array<double, 3> startPoint, double h, double m,
                  std::array<double, 3> v, double meanV, int type, double sigma, double epsilon,
                  SimulationContainer &container) {
    for (double x = 0; x < dimension[0]; x++) {
        for (double y = 0; y < dimension[1]; y++) {
            for (double z = 0; z < dimension[2]; z++) {
                // Add particle to ParticleContainer
                Particle p = Particle({x * h + startPoint[0], y * h + startPoint[1], z * h + startPoint[2]}, v, m, type,
                                      sigma, epsilon);
                container.forceInsert(p);
            }
        }
    }
}

void generateSphere3D(std::array<double, 3> center, std::array<double, 3> v, int r, double h, double m, double meanV,
                      SimulationContainer &container) {
    //first generate cube of size 2*r, then only insert particles that fit into the wanted sphere, "cutting" the sphere out of the cube
    std::array<double, 3> startPoint = {center[0] - (r * h), center[1] - (r * h), center[2] - (r * h)};
    for (double x = 0; x < (2 * r); x++) {
        for (double y = 0; y < (2 * r); y++) {
            for (double z = 0; z < (2 * r); z++) {
                Particle p = Particle({x * h + startPoint[0], y * h + startPoint[1], z * h + startPoint[2]}, v, m);
                if (ArrayUtils::L2Norm(p.getX() - center) <= (r * h)) {
                    container.insert(p);
                }

            }
        }
    }
}

void generateSphere2D(std::array<double, 3> center, std::array<double, 3> v, int r, double h, double m, double meanV,
                      SimulationContainer &container) {
    //first generate cube of size 2*r, then only insert particles that fit into the wanted sphere, "cutting" the sphere out of the cube
    std::array<double, 3> startPoint = {center[0] - (r * h), center[1] - (r * h), center[2]};
    for (double x = 0; x < (2 * r); x++) {
        for (double y = 0; y < (2 * r); y++) {
            for (double z = 0; z < 1; z++) {
                Particle p = Particle({x * h + startPoint[0], y * h + startPoint[1], z * h + startPoint[2]}, v, m);
                if (ArrayUtils::L2Norm(p.getX() - center) <= (r * h)) {
                    container.insert(p);
                }

            }
        }
    }
}

void generateSphere2D(std::array<double, 3> center, std::array<double, 3> v, int r, double h, double m, double meanV,
                      int type, SimulationContainer &container) {
//first generate cube of size 2*r, then only insert particles that fit into the wanted sphere, "cutting" the sphere out of the cube
    std::array<double, 3> startPoint = {center[0] - (r * h), center[1] - (r * h), center[2]};
    for (double x = 0; x < (2 * r); x++) {
        for (double y = 0; y < (2 * r); y++) {
            for (double z = 0; z < 1; z++) {
                Particle p = Particle({x * h + startPoint[0], y * h + startPoint[1], z * h + startPoint[2]}, v, m);
                if (ArrayUtils::L2Norm(p.getX() - center) <= (r * h)) {
                    p.setType(type);
                    container.insert(p);
                }

            }
        }
    }
}

void generateMembrane(std::array<double, 3> center, std::array<double, 3> v,std::array<int, 3> dimension, double h, const std::vector<std::array<int, 3>>& pull, double m, SimulationContainer &container) {
    int id = 0;
    Particle ps[dimension[0]][dimension[1]];
    //std::shared_ptr<Particle> ps_ptr[dimension[0]][dimension[1]];
    for (int x = 0; x < dimension[0]; x++) {
        for (int y = 0; y < dimension[1]; y++) {
            // Add particle to the grid
            ps[x][y] = Particle({x * h + center[0], y * h + center[1], center[2]}, v, m, 0, 1, 1);
            ps[x][y].setID(id);
            //ps_ptr[x][y] = std::make_shared<Particle>(ps[x][y]);
            id++;
        }
    }

    // Indicate that the particles with indices in vector pull should be pulled upwards during the simulation
    for (auto position : pull) {
        ps[position[0]][position[1]].setMembranePull();
    }

    for (int x = 0; x < dimension[0]; x++) {
        for (int y = 0; y < dimension[1]; y++) {
            if (ps[x][y].getX()[1] == 17.2) {
                std::cout << "TEST" << std::endl;
            }
            // set neighbours
            // note that only neighbours in the xy plane matter
            if (x < dimension[0]-1) {
                ps[x][y].addNeighbour(ps[x + 1][y].getID(), false);
            }
            if (x > 0) {
                ps[x][y].addNeighbour(ps[x-1][y].getID(), false);
            }
            if (y < dimension[1]-1) {
                ps[x][y].addNeighbour(ps[x][y+1].getID(), false);
            }
            if (y > 0) {
                ps[x][y].addNeighbour(ps[x][y-1].getID(), false);
            }

            if (x > 0 && y > 0) {
                ps[x][y].addNeighbour(ps[x-1][y-1].getID(), true);
            }
            if (x < dimension[0]-1 && y < dimension[1]-1) {
                ps[x][y].addNeighbour(ps[x+1][y+1].getID(), true);
            }
            if (x > 0 && y < dimension[1]-1) {
                ps[x][y].addNeighbour(ps[x-1][y+1].getID(), true);
            }
            if (x < dimension[0]-1 && y > 0) {
                ps[x][y].addNeighbour(ps[x+1][y-1].getID(), true);
            }
            // add to simulation contaienr
            container.forceInsert(ps[x][y]);
            std::cout << "X: " << x << " Y: " << y << std::endl;
        }
    }

}

void generateFromFileTest(SimulationContainer &particles, std::string f) {
    int num = 0;

    std::string tmp_string;

    auto vec = splitToString(f, '\n');
    for (auto s: vec) {
        if (num == 0) {
            std::cout << "Reading " << s << "." << std::endl;
            num = std::stoi(s);
        } else {
            if (s == "Cube:")
                continue;
            parseCube(s, particles);
        }
    }
}

void generateFromFile(SimulationContainer &particles, char *filename) {
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


        for (int i = 0; i < num * 2; i++) {
            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;

            if (tmp_string == "Cube:\r" || tmp_string == "Cube:\n") {
                getline(input_file, tmp_string);
                parseCube(tmp_string, particles);
            }
            if (tmp_string == "Sphere:\r" || tmp_string == "Sphere:\n") {
                getline(input_file, tmp_string);
                parseSphere2D(tmp_string, particles);
            }
        }
    } else {
        std::cout << "Error: could not open file " << filename << std::endl;
        exit(-1);
    }
}

void parseSphere2D(std::string str, SimulationContainer &particleContainer) {
    std::vector<std::string> strings = splitToString(str, ';');
    //Read input to arrays/doubles
    std::array<double, 3> startPoint = convertToFixedArray(splitToDouble(strings[0], ','));
    std::array<double, 3> v = convertToFixedArray(splitToDouble(strings[1], ','));
    double r = std::stoi(strings[2]);
    double h = std::stod(strings[3]);
    double meanV = std::stod(strings[4]);
    double mass = std::stod(strings[5]);

    generateSphere2D(startPoint, v, r, h, mass, meanV, particleContainer);
}

void parseCube(std::string str, SimulationContainer &particleContainer) {
    std::vector<std::string> strings = splitToString(str, ';');
    //Read input to arrays/doubles
    std::array<double, 3> startPoint = convertToFixedArray(splitToDouble(strings[0], ','));
    std::array<double, 3> v = convertToFixedArray(splitToDouble(strings[1], ','));
    std::array<double, 3> dimension = convertToFixedArray(splitToDouble(strings[2], ','));
    std::array<int, 3> dim = {(int) dimension[0], (int) dimension[1], (int) dimension[2]};
    double h = std::stod(strings[3]);
    double meanV = std::stod(strings[4]);
    double mass = std::stod(strings[5]);

    generateCube(dim, startPoint, h, mass, v, meanV, 0, 1, 1., particleContainer);
}

std::vector<std::string> splitToString(std::string str, char delimiter) {
    std::vector<std::string> strings;
    std::istringstream f(str);
    std::string s;
    while (getline(f, s, delimiter)) {
        std::cout << s << std::endl;
        strings.push_back(s);
    }
    return strings;
}

std::vector<double> splitToDouble(std::string str, char delimiter) {
    std::vector<double> strings;
    std::istringstream f(str);
    std::string s;
    while (getline(f, s, delimiter)) {
        std::cout << s << std::endl;
        strings.push_back(std::stod(s));
    }
    return strings;
}

std::array<double, 3> convertToFixedArray(std::vector<double> v) {
    return {v[0], v[1], v[2]};
}

