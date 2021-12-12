//
// Created by paulu on 03.11.2021.
//
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "outputWriter/VTKWriter.h"
#include <iostream>


ParticleContainer::ParticleContainer() = default;

Particle& ParticleContainer::getParticle(int i) {
    return particles[i];
}

std::vector<Particle>&  ParticleContainer::getVec() { return particles; }

std::vector< std::pair<Particle, Particle> > ParticleContainer::pairs() {
    std::vector< std::pair<Particle, Particle> > pairs;
    for (long unsigned int i = 0; i < particles.size(); i++) {
        auto &p1 = particles[i];

        for (long unsigned int j = i+1; j < particles.size(); j++) {
            auto &p2 = particles[j];
            pairs.emplace_back(p1, p2);
        }
    }

    return pairs;
}


void ParticleContainer::emplace_back(std::array<double, 3> &x, std::array<double, 3> &v, double &m) {
    particles.emplace_back(x, v, m);
}

int ParticleContainer::size() {
    return particles.size();
}

void ParticleContainer::calculateF(ForceCalculation *algorithm) {
    for (auto &p: particles) {
        p.setOldF(p.getF());
        p.setF({0, 0, 0});
    }

    for (long unsigned int i = 0; i < particles.size(); i++) {
        auto &p1 = particles[i];

        for (long unsigned int j = i + 1; j < particles.size(); j++) {
            auto &p2 = particles[j];

            std::array<double, 3> vec = algorithm->calculateF(p1, p2);

            //F_p1 = -F_p2
            //p1.setF(ArrayUtils::elementWisePairOp(p1.getF(), vec, add));
            p1.setF(p1.getF() + vec);
            //p2.setF(ArrayUtils::elementWisePairOp(p2.getF(), vec, sub));
            p2.setF(p2.getF() - vec);
        }
    }
}

void ParticleContainer::calculateX(double delta_t) {
    for (auto &p: particles) {
        std::array<double, 3> res = (1 / (2 * p.getM())) * p.getOldF();
        res = (delta_t * delta_t) * res;

        std::array<double, 3> res2 = delta_t * p.getV();
        res = res + res2;
        p.setX(res + p.getX());
    }
}

void ParticleContainer::calculateV(double delta_t) {
    for (auto &p: particles) {
        std::array<double, 3> res = p.getF() + p.getOldF();
        res = (1 / (2 * p.getM())) * res;
        res = delta_t * res;
        res = res + p.getV();
        p.setV(res);
    }
}

void ParticleContainer::plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(particles.size());
    //writer.initializeOutput(particles.size());
    for (auto &p : particles) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}

void ParticleContainer::simulate(ForceCalculation *algorithm, double delta_t) {

    std::cout << "Number of particles: " << particles.size() << std::endl;
    // calculate new x
    calculateX(delta_t);

    // calculate new f
    calculateF(algorithm);

    // calculate new v
    calculateV(delta_t);
}

void ParticleContainer::insert(Particle& p) {
    particles.push_back(p);
}

void ParticleContainer::addBrownianMotion(double averageV, int dimension) {
    for (auto &p : particles) {
            p.setV(p.getV() + maxwellBoltzmannDistributedVelocity(averageV, dimension));
    }
}

double ParticleContainer::calcKineticEnergy() {
    double energy = 0;
    for (auto &p : particles) {
        double scalar_product = p.getV()[0] * p.getV()[0] + p.getV()[1] * p.getV()[1] + p.getV()[2] * p.getV()[2];
        energy += p.getM() * scalar_product / 2;
    }
    return energy;
}

int ParticleContainer::numberParticles() {
    return particles.size();
}

void ParticleContainer::scaleVelocity(double scale) {
    for (auto &p : particles) {
        p.setV(scale * p.getV());
    }
}
