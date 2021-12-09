#include "FileReader.h"
#include "ParticleGenerator.h"
#include "container/LinkedCells.h"
#include "forceCalculation/Gravitation.h"
#include "forceCalculation/LennardJones.h"
#include "XMLReader/MolSim-pimpl.h"
#include "Log.h"
#include <iostream>
#include <chrono>
#include <getopt.h>


/**** forward declaration of the calculation functions ****/
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

std::array<int, 3> dim = {300, 300, 0};
double mesh = 3;
double cutOff = 3;
SimulationContainer *particles = new LinkedCells(dim, mesh, cutOff, sigma);
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
    std::string file;
    int c;
    while ((c = getopt(argc, argsv, "hx:f:s:e:w:a:b")) != -1) {
        if (c == 'h') {
            printHelp();
            return 0;
        }
        if (c == 'x') {
            try
            {
                // Instantiate individual parsers.
                //
                ::molsim_pimpl molsim_p;
                ::input_file_pimpl input_file_p;
                ::delta_t_pimpl delta_t_p;
                ::end_time_pimpl end_time_p;
                ::output_step_pimpl output_step_p;
                ::algorithm_pimpl algorithm_p;
                ::benchmark_pimpl benchmark_p;

                // Connect the parsers together.
                //
                molsim_p.parsers (input_file_p,
                                  delta_t_p,
                                  end_time_p,
                                  output_step_p,
                                  algorithm_p,
                                  benchmark_p);

                // Parse the XML document.
                //
                ::xml_schema::document doc_p (molsim_p, "molsim");

                molsim_p.pre ();

                doc_p.parse (optarg);
                molsim_p.post_molsim ();

                file = input_file_p.get_input_file();
                delta_t = delta_t_p.get_delta_t();
                end_time = end_time_p.get_end_time();
                outputStep = output_step_p.get_output_step();
                if (std::string("sv") == algorithm_p.get_algorithm()) {
                    algorithm = new Gravitation();
                } else if (std::string("lj") == algorithm_p.get_algorithm()) {
                    algorithm = new LennardJones(epsilon, sigma);
                    cuboids = true;
                }
                if (benchmark_p.get_benchmark()==std::string("yes")) {
                    benchmark_active = true;
                }
            }
            catch (const ::xml_schema::exception& e)
            {
                std::cerr << e << std::endl;
                return 1;
            }
            catch (const std::ios_base::failure&)
            {
                std::cerr << optarg << ": error: io failure" << std::endl;
                return 1;
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
        char* file_char = &file[0];
        generateFromFile(*particles, file_char);
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
        particles->simulate(algorithm, delta_t);

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


        current_time += delta_t;
    }

    if (benchmark_active) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time (including reading the input file and setting everything up) = "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
        std::cout << "Time (only calculations) = "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - beginAfterIO).count() << "[s]" << std::endl;
    }

    LOGC_INFO("output written. Terminating...");
    return 0;
}

void printHelp() {
    LOGC_INFO("Usage: -x xml file -f filename -a algorithm -s step size -e end time -w Output step size -b activate benchmark");
    LOGC_WARN("xmlfile: if you choose to use xml input, you are not required to fill others");
    LOGC_WARN("Filename: path to the input file (required)");
    LOGC_WARN("Step size: size of a timestep in the simulation (optional)");
    LOGC_WARN("Output step size: Every this often steps an output file will be generated (optional)");
    LOGC_WARN("Algorithm: Algorithm used for force calculations (required)");
    LOGC_INFO("Possible Algorithms: sv (Stoermer Verlet), lj (Lennard Jones, generates cuboids)");
    LOGC_INFO("Benchmark: Disables writing files and benchmarks the program");
}
