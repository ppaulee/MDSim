//
// Created by paulu on 03.11.2021.
//

#pragma once

#include "Particle.h"
#include "vector"


class ParticleContainer {
private:
    std::vector<Particle> particles;

public:
    Particle& getParticle(int i);

    /**
     * @return Reference to Vector
     */
    std::vector<Particle>& getVec();

    /**
     * The order of the pair does not matter e.g. (a,b) = (b,a)
     *
     * @return Vector of all pairs
     */
    std::vector< std::pair<Particle,Particle> > pairs();

    void emplace_back(std::array<double, 3> &x, std::array<double, 3> &v, double &m);

    int size();
    ParticleContainer();
};

