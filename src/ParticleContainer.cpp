//
// Created by paulu on 03.11.2021.
//

#include "ParticleContainer.h"



ParticleContainer::ParticleContainer(std::vector<Particle> &p) {
    particles = p;
}

ParticleContainer::ParticleContainer() {}

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

void ParticleContainer::applyPairs(pairFunction func) {
    std::vector< std::pair<Particle, Particle> > pairs;
    for (long unsigned int i = 0; i < particles.size(); i++) {
        auto &p1 = particles[i];

        for (long unsigned int j = i+1; j < particles.size(); j++) {
            auto &p2 = particles[j];
            func(p1, p2);
        }
    }
}

void ParticleContainer::apply(singleFunction func) {
    for (auto &p : particles) {
        func(p);
    }
}

void ParticleContainer::emplace_back(std::array<double, 3> &x, std::array<double, 3> &v, double &m) {
    particles.emplace_back(x, v, m);
}

int ParticleContainer::size() {
    return particles.size();
}


