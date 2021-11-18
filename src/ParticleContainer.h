//
// Created by paulu on 03.11.2021.
//

#pragma once

#include "Particle.h"
#include "vector"
#include <functional>


class ParticleContainer {
private:
    std::vector<Particle> particles;

public:
    Particle& getParticle(int i);

    using iterator = std::vector<Particle>::iterator;
    using const_iterator = std::vector<Particle>::const_iterator;

    iterator begin()
    {
        return particles.begin();
    }

    iterator end()
    {
        return particles.end();
    }

    const_iterator cbegin() const
    {
        return particles.cbegin();
    }

    const_iterator cend() const
    {
        return particles.cend();
    }

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

