//
// Created by paulu on 03.11.2021.
//

#pragma once

#include "Particle.h"
#include "vector"
#include "memory"

typedef void (* pairFunction)(Particle a, Particle b);
typedef void (* singleFunction)(Particle p);

class ParticleContainer {
private:
    std::vector<Particle> particles;

public:
    Particle& getParticle(int i);
    std::vector<Particle>& getVec();
    std::vector< std::pair<Particle,Particle> > pairs();
    void applyPairs(pairFunction func);
    void apply(singleFunction func);
    void emplace_back(std::array<double, 3> &x, std::array<double, 3> &v, double &m);
    int size();
    ParticleContainer(std::vector<Particle> &p);
    ParticleContainer();
};

