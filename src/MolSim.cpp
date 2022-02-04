#include "FileReader.h"
#include "ParticleGenerator.h"
#include "Thermostats.h"
#include "container/LinkedCells.h"
#include "container/ParallelLinkedCells.h"
#include "forceCalculation/Gravitation.h"
#include "forceCalculation/LennardJones.h"
#include "forceCalculation/HarmonicPotential.h"
#include "XMLReader/driver.h"
#include "XMLReader/library.h"
#include "Log.h"
#include "CheckpointFileWriter.h"
#include <iostream>
#include <chrono>
#include <getopt.h>
#include <sstream>
#include <string>
#include "densityAndVelocityProfile/DensityAndVelocityProfile.h"
#include "densityAndVelocityProfile/DensityAndVelocityToCSVWriter.h"

/**** forward declaration of the calculation functions ****/
void printHelp();

constexpr double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;
int outputStep = 100;
// Parameters for Lennard Jones Potential
double epsilon = 1;
double sigma = 1.2;
// Brownian Motion average velocity
double averageV = 2 * sqrt(10);
std::array<int, 3> dim = {60, 60, 60};
//std::array<int, 3> dim = {64, 36, 0};
double mesh = 2.5;
double cutOff = 2.5;
std::array<int, 3> bound = {2, 1, 1};
SimulationContainer *particles = new ParallelLinkedCells(dim, mesh, cutOff, -12.44, bound, 2);
//SimulationContiner
//std::array<int, 3> dim = {148, 148, 148};

// Stores the algorithm used for force calculation between 2 particles
ForceCalculation *algorithm = nullptr;

//thermostats variables

double initial_temp = 0.5;
int stepSize = 1000;
double target_temp = 0.5;
double max_delta_temp = -1;

//simulationType
bool simuMembrane = true;
bool simuFlow = true;

//DensityAndVelocityProfiler
DensityAndVelocityProfile *profile = nullptr;
int profileBinSize;
double profileH;
double profileW = 3;
std::array<double, 3> profileStartPoint;

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
    std::string file;
    int c;
    while ((c = getopt(argc, argsv, "hx:f:s:e:w:a:b")) != -1) {
        if (c == 'h') {
            printHelp();
            return 0;
        }
        if (c == 'x') {
            library::molsim m = generateFromXML(optarg);

            /**
             * Preparing values from xml for simulation
             */
            file = m.input();
            delta_t = m.delta_t();
            end_time = m.endtime();
            outputStep = m.outputStep();
            epsilon = m.epsilon();
            sigma = m.sigma();
            averageV = m.averageV();
            //thermostats
            initial_temp = m.thermostats().initialTemperature();
            stepSize = m.thermostats().stepSize();
            target_temp = m.thermostats().targetTemperature();
            max_delta_temp = m.thermostats().maxDelta();
            //simulationContainer
            std::cout << "container: " << m.simu().containerAlgorithm() << "\n";
            if (m.simu().containerAlgorithm() == std::string("linkedCells")) {
                dim = {m.simu().dimension().x(), m.simu().dimension().y(), m.simu().dimension().z()};
                bound = {m.simu().boundaryConditions().x(), m.simu().boundaryConditions().y(),
                         m.simu().boundaryConditions().z()};
                //parallelization strategy
                int parallelstrategy = m.parStrat();
                if (parallelstrategy == 0)
                    particles = new LinkedCells(dim, m.simu().mesh(), m.simu().cutOff(), m.gravity(), bound);
                if (parallelstrategy == 1) {
                    particles = new ParallelLinkedCells(dim, m.simu().mesh(), m.simu().cutOff(), m.gravity(), bound, 1);
                }
                if (parallelstrategy == 1) {
                    particles = new ParallelLinkedCells(dim, m.simu().mesh(), m.simu().cutOff(), m.gravity(), bound, 2);
                }
            } else if (m.simu().containerAlgorithm() == std::string("naiv")) {
                algorithm = new Gravitation();
            } else {
                LOGC_ERROR("Error: container algorithm doesnt fit, -naiv or -linkedCells");
            }
            //algorithm
            if (m.algorithm() == std::string("sv")) {
                algorithm = new Gravitation();
            } else if (m.algorithm() == std::string("lj")) {
                algorithm = new LennardJones(epsilon, sigma, 0);
            } else {
                LOGC_ERROR("Error:algorithm doesnt fit, -sv or -lj");
            }
            //benchmark
            if (m.benchmark() == std::string("yes")) {
                benchmark_active = true;
            }

            if (m.simulationType() == std::string("simuMembrane")) {
                simuFlow = false;
            } else if (m.simulationType() == std::string("simuFlow")) {
                simuMembrane = false;
            } else {
                simuMembrane = false;
                simuFlow = false;
            }


            //  generate particle from xml
            int currentType = 0;
            std::cout << "nr cubes " << m.particles().cube().size() << "\n";
            for (long unsigned int i = 0; i < m.particles().cube().size(); i++) {
                library::Cube cube = m.particles().cube()[i];
                std::array<int, 3> dim = {cube.dimension().x(), cube.dimension().y(), cube.dimension().z()};
                std::array<double, 3> startPoint = {cube.startPoint().x(), cube.startPoint().y(),
                                                    cube.startPoint().z()};
                std::array<double, 3> velocity = {cube.velocity().x(), cube.velocity().y(), cube.velocity().z()};
                generateCube(dim, startPoint, cube.h(), cube.m(), velocity, averageV, currentType, cube.sigma(),
                             cube.epsilon(), *particles);

                //profiler for flow simulation
                if (simuFlow && !cube.fixed()) {
                    profileBinSize = dim[0];
                    profileH = cube.h();
                    profileStartPoint = startPoint;
                    profile = new DensityAndVelocityProfile(profileBinSize, profileW, profileH, profileStartPoint);
                }

                currentType++;
            }
            for (long unsigned int i = 0; i < m.particles().sphere().size(); i++) {
                library::Sphere sphere = m.particles().sphere()[i];
                std::array<double, 3> center = {sphere.center().x(), sphere.center().y(), sphere.center().z()};
                std::array<double, 3> velocity = {sphere.velocity().x(), sphere.velocity().y(),
                                                  sphere.velocity().z()};
                generateSphere2D(center, velocity, sphere.radius(), sphere.h(), sphere.m(), averageV, currentType,
                                 *particles);
                currentType++;
                /**
                * end setting values from xml for simulation
                */
            }
            if (dim[2] == 0) {
                particles->addBrownianMotion(averageV, 2);
            } else {
                particles->addBrownianMotion(averageV, 3);
            }
            break;
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
                algorithm = new Gravitation();
            } else if (std::string("lj") == optarg) {
                algorithm = new LennardJones(epsilon, sigma, 0);
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
    if (file == "") {
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
        //char *file_char = &file[0];
        //generateFromFile(*particles, file_char);

        //generateCube({250, 50, 1}, {7.5, 8, 0}, 1.2, 1, {0, 0, 0}, 0, 0, 1.2, 1, *particles);

        //generateCube({50, 14, 1}, {5.6, 7, 0}, 1.2, 1, {0, 0, 0}, 2 * sqrt(10), 0, 1, 1, *particles);
        //generateCube({50, 14, 1}, {5.6, 24.5, 0}, 1.2, 2, {0, 0, 0}, 2 * sqrt(5), 1, 0.9412, 1, *particles);


        //3D rayleigh taylor
        generateCube({50, 20, 50}, {6.2, 4.2, 6.2}, 1.2, 1, {0, 0, 0}, 2 * sqrt(10), 0, 1.2, 1, *particles);
        generateCube({50, 20, 50}, {6.2, 32, 6.2}, 1.2, 2, {0, 0, 0}, 2 * sqrt(5), 1, 1.1, 1, *particles);

        //generateCube({4, 4, 1}, {5.6, 7, 0}, 1.2, 1, {0, 0, 0}, 1.2, 0, 1, 1, *particles);
        //generateCube({4, 4, 1}, {5.6, 9, 0}, 1.2, 2, {0, 0, 0}, 1.2, 1, 0.9412, 1, *particles);

        //generateCube({2, 2, 1}, {15, 15, 0}, 1.1225, 1, {-10, 0, 0}, averageV, 0, 1, 5, *particles);

        particles->addBrownianMotion(averageV, 3);
    }

    if (simuMembrane) {
        particles->setMembraneSimulation();
        std::array<double, 3> center = {15, 15, 1.5};
        std::vector<std::array<int, 3>> pull;

        pull.push_back({17, 24, 0});
        pull.push_back({17, 25, 0});
        pull.push_back({18, 24, 0});
        pull.push_back({18, 25, 0});

        //pull.push_back({4,4,0});
        std::array<double, 3> v = {0, 0, 0};
        std::array<int, 3> dim_ = {50, 50, 1};
        //std::array<int, 3> dim_ = {10,10,1};
        generateMembrane(center, v, dim_, 2.2, pull, 1, *particles);
        particles->initMembrane();
    }

    if (simuFlow) {
        particles->setNanoScaleFlowSimulation();
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

    int dimensions = 3;
    if (dim[0] == 0) {
        dimensions = 2;
    }

    Thermostats *thermostat;
    if (!simuMembrane) {
        thermostat = new Thermostats(*particles, delta_t, end_time, initial_temp, stepSize, dimensions, target_temp,
                                     max_delta_temp);

        thermostat->calcCurrentTemperature(*particles);
        thermostat->adjustTemperature(*particles, current_time);
    }

    DensityAndVelocityToCSVWriter csvWriter("profile.csv");

    bool pullState = true;
    while (current_time < end_time) {
        std::cout << "Number particles: " << particles->numberParticles() << "\n";
        if (iteration >= 15000) {
            pullState = false;
        }
        if (simuMembrane) {
            particles->simulateMembrane(delta_t, pullState);
        } else {
            particles->simulate(algorithm, delta_t);
        }

        // Control temperature
        if (!simuMembrane) {
            thermostat->calcCurrentTemperature(*particles);
            thermostat->adjustTemperature(*particles, current_time);
            thermostat->calcCurrentTemperature(*particles);
        }

        iteration++;
        if (iteration % outputStep == 0 && !benchmark_active) {
            particles->plotParticles(iteration);
        }

        if (!benchmark_active) {
            LOGC_TRACE("Iteration {} finished.", iteration);
        } else if (iteration % outputStep == 0) {
            std::cout << "Iteration " << iteration << "finished." << std::endl;
            std::cout << "Current time " << current_time << std::endl;
        }

        if (simuFlow && iteration % 1000 == 0) {
            profile->calculateDnV(*particles);
            csvWriter.writeToCSV(*profile, iteration);
        }

        current_time += delta_t;
    }

    //std::cout << "TEMP: " << thermostat->getCurrentTemperature() << std::endl;
    if (benchmark_active) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time (including reading the input file and setting everything up) = "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
        std::cout << "Time (only calculations) = "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - beginAfterIO).count() << "[s]" << std::endl;
    }
    CheckpointFileWriter w = CheckpointFileWriter();
    w.writeFile(*particles, "stableLiquid");
    LOGC_INFO("output written. Terminating...");
    return 0;
}

void printHelp() {
    LOGC_INFO(
            "Usage: -x xml file -f filename -a algorithm -s step size -e end time -w Output step size -b activate benchmark");
    LOGC_WARN("xmlfile: if you choose to use xml input, you are not required to fill others");
    LOGC_WARN("Filename: path to the input file (required)");
    LOGC_WARN("Step size: size of a timestep in the simulation (optional)");
    LOGC_WARN("Output step size: Every this often steps an output file will be generated (optional)");
    LOGC_WARN("Algorithm: Algorithm used for force calculations (required)");
    LOGC_INFO("Possible Algorithms: sv (Stoermer Verlet), lj (Lennard Jones, generates cuboids)");
    LOGC_INFO("Benchmark: Disables writing files and benchmarks the program");
}

