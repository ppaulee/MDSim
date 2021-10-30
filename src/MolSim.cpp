
#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <list>
#include <cmath>


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
void plotParticles(int iteration);

/**
 * @param p1 Postition 1
 * @param p2  Position 2
 * @return difference between both vectors
 */
std::array<double, 3> posDifference(std::array<double, 3> p1, std::array<double, 3> p2);

/**
 *
 * @param p Vector
 * @return euclidean norm
 */
double eudlideanNorm(std::array<double, 3> p);

double add(double a, double b);
double multiply(double a, double b);
double divide(double a, double b);
double sub(double a, double b);

constexpr double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;


std::vector<Particle> particles;

int main(int argc, char *argsv[]) {

    std::cout << "Hello from MolSim for PSE!" << std::endl;
    if (argc == 4) {
        end_time = atof(argsv[2]);
        delta_t = atof(argsv[3]);
    } else if (argc != 2) {
        std::cout << "Erroneous programme call! " << std::endl;
        std::cout << "./molsym filename" << std::endl;
    }

    FileReader fileReader;
    fileReader.readFile(particles, argsv[1]);

    double current_time = start_time;

    int iteration = 0;

    // for this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
        // calculate new x
        calculateX();
        // calculate new f
        calculateF();
        // calculate new v
        calculateV();

        iteration++;
        if (iteration % 10 == 0) {
            plotParticles(iteration);
        }
        std::cout << "Iteration " << iteration << " finished." << std::endl;

        current_time += delta_t;
    }

    std::cout << "output written. Terminating..." << std::endl;
    return 0;
}


void calculateF() {
    for(auto &p : particles) {
        p.setOldF(p.getF());
        p.setF({0,0,0});
    }

    for (long unsigned int i = 0; i < particles.size(); i++) {
        auto &p1 = particles[i];

        for (long unsigned int j = i+1; j < particles.size(); j++) {
            auto &p2 = particles[j];

            double f = (p1.getM() * p2.getM()) / pow(ArrayUtils::L2Norm(ArrayUtils::elementWisePairOp(p1.getX(), p2.getX(), sub)), 3);
            std::array<double, 3> vec = ArrayUtils::elementWisePairOp(p2.getX(),p1.getX(), sub);
            vec = ArrayUtils::elementWiseScalarOp(f, vec, multiply);

            //F_p1 = -F_p2
            p1.setF(ArrayUtils::elementWisePairOp(p1.getF(), vec, add));
            p2.setF(ArrayUtils::elementWisePairOp(p2.getF(), vec, sub));
        }
    }
}

void calculateX() {
    for (auto &p: particles) {
        std::array<double, 3> res = ArrayUtils::elementWiseScalarOp((2*p.getM()), p.getOldF(), divide);
        res = ArrayUtils::elementWiseScalarOp((delta_t*delta_t), res, multiply);

        std::array<double, 3> res2 = ArrayUtils::elementWiseScalarOp(delta_t, p.getV(), multiply);
        res = ArrayUtils::elementWisePairOp(res, res2, add);
        p.setX(ArrayUtils::elementWisePairOp(res, p.getX(), add));
    }
}

void calculateV() {
    for (auto &p: particles) {
        std::array<double, 3> res = ArrayUtils::elementWisePairOp(p.getF(), p.getOldF(), add);
        res = ArrayUtils::elementWiseScalarOp((2*p.getM()), res, divide);
        res = ArrayUtils::elementWiseScalarOp(delta_t, res , multiply);
        res = ArrayUtils::elementWisePairOp(res, p.getV(), add);
        p.setV(res);
    }
}

double add(double a, double b){return a+b;};
double multiply(double a, double b){return a*b;};
double divide(double a, double b){return b/a;};
double sub(double a, double b){return a-b;};


void plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(particles.size());
    for(auto &p : particles) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}
