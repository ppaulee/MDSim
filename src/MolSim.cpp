#include "FileReader.h"
#include "ParticleGenerator.h"
#include "Thermostats.h"
#include "container/LinkedCells.h"
#include "forceCalculation/Gravitation.h"
#include "forceCalculation/LennardJones.h"
#include "XMLReader/MolSimImpl.h"
#include "XMLReader/library.h"
#include "Log.h"
#include "CheckpointFileWriter.h"
#include <iostream>
#include <chrono>
#include <getopt.h>
#include <sstream>
#include <string>


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
double averageV = 0.7;

std::array<int, 3> dim = {302, 180, 0};
double mesh = 3;
double cutOff = 3;
std::array<int, 3> bound = {2, 1, 0};
SimulationContainer *particles = new LinkedCells(dim, mesh, cutOff, -12.44, bound);
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
            try {
                // Instantiate individual parsers.
                //
                ::molsim_pimpl molsim_p;
                ::xml_schema::string_pimpl string_p;
                ::time_pimpl time_p;
                ::xml_schema::int_pimpl int_p;
                ::epsilon_pimpl epsilon_p;
                ::sigma_pimpl sigma_p;
                ::xml_schema::double_pimpl double_p;
                ::algorithm_pimpl algorithm_p;
                ::simulationContainer_pimpl simulationContainer_p;
                ::boundaryConditions_pimpl boundaryConditions_p;
                ::dimension_pimpl dimension_p;
                ::containerAlgorithm_pimpl containerAlgorithm_p;
                ::particles_pimpl particles_p;
                ::Cube_pimpl Cube_p;
                ::point_pimpl point_p;
                ::velocity_pimpl velocity_p;
                ::Sphere_pimpl Sphere_p;
                ::thermostats_pimpl thermostats_p;
                ::benchmark_pimpl benchmark_p;

                // Connect the parsers together.
                //
                molsim_p.parsers(string_p,
                                 time_p,
                                 time_p,
                                 int_p,
                                 epsilon_p,
                                 sigma_p,
                                 double_p,
                                 double_p,
                                 algorithm_p,
                                 simulationContainer_p,
                                 particles_p,
                                 thermostats_p,
                                 benchmark_p);

                simulationContainer_p.parsers(boundaryConditions_p,
                                              dimension_p,
                                              double_p,
                                              double_p,
                                              containerAlgorithm_p);

                dimension_p.parsers(int_p,
                                    int_p,
                                    int_p);

                particles_p.parsers(Cube_p,
                                    Sphere_p);

                Cube_p.parsers(dimension_p,
                               point_p,
                               double_p,
                               double_p,
                               velocity_p,
                               epsilon_p,
                               sigma_p);

                point_p.parsers(double_p,
                                double_p,
                                double_p);

                velocity_p.parsers(double_p,
                                   double_p,
                                   double_p);

                Sphere_p.parsers(point_p,
                                 double_p,
                                 double_p,
                                 double_p,
                                 velocity_p,
                                 epsilon_p,
                                 sigma_p);

                thermostats_p.parsers(double_p,
                                      double_p,
                                      double_p,
                                      int_p);

                // Parse the XML document.
                //
                ::xml_schema::document doc_p(molsim_p, "molsim");

                molsim_p.pre();
                doc_p.parse(optarg);
                library::molsim m = molsim_p.post_molsim();

                //all parameters can be accessed from m from here
                //
                //take a look at the xml schema how to access the parameters
                //

                file = m.input();
                delta_t = m.delta_t();
                end_time = m.endtime();
                outputStep = m.outputStep();
                epsilon = m.epsilon();
                sigma = m.sigma();
                averageV = m.averageV();
                //simulationContainer
                std::cout << "container: " << m.simu().containerAlgorithm() << "\n";
                if (m.simu().containerAlgorithm() == std::string("linkedCells")) {
                    dim = {m.simu().dimension().x(), m.simu().dimension().y(), m.simu().dimension().z()};
                    std::array<int, 3> bounds = {1, 1, 1};
                    particles = new LinkedCells(dim, m.simu().mesh(), m.simu().cutOff(), m.gravity(), bounds);
                } else if (m.simu().containerAlgorithm() == std::string("naiv")) {
                    //TODO make naive container here
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
                }

                /*
                 *  Thermostats
                 *
                int dimensions = 3;
                if (dim[0] == 0) {
                    dimensions = 2;
                auto thermostat = new Thermostats(*particles, delta_t, end_time, m.thermostats().initialTemperature(), m.thermostats().stepSize(), dimensions, m.thermostats().targetTemperature(), m.thermostats().maxDelta());
                **/

                /* sigma and particles can be accessed pro particle

                    m.particles().cube()[0].sigma()
                    m.particles().cube()[0].epsilon()

                 **/

                /*  added gravity parameter
                 *
                 * m.gravity()
                 *
                 * */
            }
            catch (const ::xml_schema::exception &e) {
                std::cerr << e << std::endl;
                return 1;
            }
            catch (const std::ios_base::failure &) {
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

        generateCube({250, 50, 1}, {7.5, 8, 0}, 1.2, 1, {0, 0, 0}, 0, 0, 1.2, 1, *particles);

        //generateCube({50, 14, 1}, {5.6, 7, 0}, 1.2, 1, {0, 0, 0}, 1.2, 0, 1, 1, *particles);
       // generateCube({50, 14, 1}, {5.6, 24, 0}, 1.2, 2, {0, 0, 0}, 1.2, 1, 0.9412, 1, *particles);

        //generateCube({2, 2, 1}, {15, 15, 0}, 1.1225, 1, {-10, 0, 0}, averageV, 0, 1, 5, *particles);

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

    int dimensions = 3;
    if (dim[0] == 0) {
        dimensions = 2;
    }

    // TODO set to XML values
    /*
    <xsd:complexType name="thermostats">
        <xsd:sequence>
            <xsd:element name="initialTemperature" type="xsd:double"/>
            <xsd:element name="targetTemperature" type="xsd:double"/>
            <xsd:element name="maxDelta" type="xsd:double"/>
            <xsd:element name="stepSize" type="xsd:int"/>
        </xsd:sequence>
    </xsd:complexType>
     **/
    double initial_temp = 0.5;
    int stepSize = 1000;
    double target_temp = 0.5;
    // TODO Dieser Parameter ist optional, wenn nicht angegeben wird -1 übergeben
    double max_delta_temp = -1;

    auto thermostat = new Thermostats(*particles, delta_t, end_time, initial_temp, stepSize, dimensions, target_temp,
                                      max_delta_temp);

    thermostat->calcCurrentTemperature(*particles);
    thermostat->adjustTemperature(*particles, current_time);

    while (current_time < end_time) {
        std::cout << "Number particles: " << particles->numberParticles() << "\n";
        particles->simulate(algorithm, delta_t);
        // Control temperature
        thermostat->calcCurrentTemperature(*particles);
        thermostat->adjustTemperature(*particles, current_time);
        thermostat->calcCurrentTemperature(*particles);

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
