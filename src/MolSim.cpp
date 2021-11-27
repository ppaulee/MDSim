#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "ForceCalculation.h"
#include "StoermerVerlet.h"
#include "LennardJones.h"
#include "Log.h"




#include <iostream>
#include <cmath>
#include <getopt.h>
#include <chrono>


/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
//void plotParticles(int iteration);

//double add(double a, double b);

//double multiply(double a, double b);

//double divide(double a, double b);

//double sub(double a, double b);

void printHelp();

constexpr double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;
int outputStep = 100;
// Parameters for Lennard Jones Potential
double epsilon = 5;
double sigma = 1;
// Brownian Motion average velocity
double averageV = 0.1;

std::array<int, 3> dim = {300,300,0};
double mesh = 1;
double cutOff = 1.5;
LinkedCells* particles = new LinkedCells(dim,mesh,cutOff);
// Stores the algorithm used for force calculation between 2 particles
ForceCalculation *algorithm = nullptr;

// Variables for the benchmark
bool benchmark_active = false;
std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point beginAfterIO;

int main(int argc, char *argsv[]) {
    MolSim::Log::Init();
    LOGC_INFO("Hello from MolSim for PSE!");
    if (argc <= 1) {
        printHelp();
        return 1;
    }

    bool cuboids = false;
    char *file = nullptr;
    int c;
    while ((c = getopt(argc, argsv, "hf:s:e:w:a:b")) != -1) {
        if (c == 'h') {
            printHelp();
            return 0;
        }
        if (c == 'f') {
            file = optarg;
        }
        if (c == 's') {
            delta_t = atof(optarg);
        }
        if (c == 'e') {
            end_time = atof(optarg);
        }
        if (c == 'w') {
            outputStep = std::stoi(optarg);
        }
        if (c == 'a') {
            if (std::string("sv") == optarg) {
                algorithm = new StoermerVerlet();
            } else if (std::string("lj") == optarg) {
                algorithm = new LennardJones(epsilon, sigma);
                cuboids = true;
            }
        }
        if (c == 'b') {
            benchmark_active = true;
        }
    }
    if (benchmark_active) {
        begin = std::chrono::steady_clock::now();
    }

    if (file == nullptr) {
        LOGF_ERROR("Error: Path to file is missing, use -h for help");
        return 1;
    }
    if (algorithm == nullptr) {
        LOGC_ERROR("Error: Algorithm missing or erroneous algorithm argument, use -h for help");
        return 1;
    }
    if (!cuboids) {
        /**
            FileReader fileReader;
            fileReader.readFile(particles, file);
         */
    } else {
        generateFromFile(*particles, file);
        //generateCube({40, 8, 1}, {0, 0, 0}, 1.1225, 1, {0, 0, 0}, averageV, particles);
        //generateCube({8, 8, 1}, {15, 15, 0}, 1.1225, 1, {0, -10, 0}, averageV, particles);

        particles->addBrownianMotion(averageV, 2);
    }

    if (benchmark_active) {
        beginAfterIO = std::chrono::steady_clock::now();
    }

    //std::cout << "NUmber of particles: " << particles.getVec().size() << std::endl;
    double current_time = start_time;

    int iteration = 1;
    // for this loop, we assume: current x, current f and current v are known

    if (!benchmark_active) {
        particles->plotParticles(0);
    }

    while (current_time < end_time) {
        particles->simulate(delta_t, algorithm);

        iteration++;
        if (iteration % outputStep == 0 && !benchmark_active) {
            particles->plotParticles(iteration);
        }

        if (!benchmark_active) {
            LOGC_TRACE("Iteration {} finished.", iteration);
        } else if (iteration % outputStep == 0) {
            std::cout << "Iteration " << iteration << "finished." << std::endl;
            std::cout << "Current time " << current_time  << std::endl;
        }


        current_time += delta_t;
    }

    ///std::cout << "NUmber of particles: " << particles.getVec().size() << std::endl;

    if (benchmark_active) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time (including reading the input file and setting everything up) = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << "[s]" << std::endl;
        std::cout << "Time (only calculations) = " << std::chrono::duration_cast<std::chrono::seconds> (end - beginAfterIO).count() << "[s]" << std::endl;
    }

    LOGC_INFO("output written. Terminating...");
    return 0;
}

void printHelp() {
    LOGC_INFO("Usage: -f filename -a algorithm -s step size -e end time -w Output step size -b activate benchmark");
    LOGC_WARN("Filename: path to the input file (required)");
    LOGC_WARN("Step size: size of a timestep in the simulation (optional)");
    LOGC_WARN("Output step size: Every this often steps an output file will be generated (optional)");
    LOGC_WARN("Algorithm: Algorithm used for force calculations (required)");
    LOGC_INFO("Possible Algorithms: sv (Stoermer Verlet), lj (Lennard Jones, generates cuboids)");
    LOGC_INFO("Benchmark: Disables writing files and benchmarks the program");
}

/**
void calculateF() {
    for (auto &p: particles) {
        p.setOldF(p.getF());
        p.setF({0, 0, 0});
    }

    for (long unsigned int i = 0; i < particles.getVec().size(); i++) {
        auto &p1 = particles.getParticle(i);

        for (long unsigned int j = i + 1; j < particles.getVec().size(); j++) {
            auto &p2 = particles.getParticle(j);

            std::array<double, 3> vec = algorithm->calculateF(p1, p2);

            //F_p1 = -F_p2
            //p1.setF(ArrayUtils::elementWisePairOp(p1.getF(), vec, add));
            p1.setF(p1.getF() + vec);
            //p2.setF(ArrayUtils::elementWisePairOp(p2.getF(), vec, sub));
            p2.setF(p2.getF() - vec);
        }
    }
}

void calculateX() {
    for (auto &p: particles) {
        std::array<double, 3> res = (1 / (2 * p.getM())) * p.getOldF();
        res = (delta_t * delta_t) * res;

        std::array<double, 3> res2 = delta_t * p.getV();
        res = res + res2;
        p.setX(res + p.getX());
    }
}

void calculateV() {
    for (auto &p: particles) {
        std::array<double, 3> res = p.getF() + p.getOldF();
        res = (1 / (2 * p.getM())) * res;
        res = delta_t * res;
        res = res + p.getV();
        p.setV(res);
    }
}



void plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    //writer.initializeOutput(particles.size());
    for (auto &p : particles) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}
**/
