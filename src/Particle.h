/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>
#include <math.h>
#include <vector>
#include <memory>

class Particle {

private:
    /**
     * Position of the particle
     */
    std::array<double, 3> x;

    /**
     * Velocity of the particle
     */
    std::array<double, 3> v;

    /**
     * Force effective on this particle
     */
    std::array<double, 3> f;

    std::array<double, 3> forceBuffer;

    /**
     * Stores all neighbours in the xy plane (only for membrane)
     */
    std::vector<std::reference_wrapper<Particle>> neighbours;
    std::vector<std::reference_wrapper<Particle>> neighboursDiag;

    /**
     * Indicates whether this particle should be pulled up in the membrane simulation
     */
    bool membranePull;

private:

    /**
     * Force which was effective on this particle
     */
    std::array<double, 3> old_f;

    /**
     * Mass of this particle
     */
    double m;

    /**
     * Tells if the particle is marked
     */
    bool marked;

    /**
     * Type of the particle. Use it for whatever you want (e.g. to separate
     * molecules belonging to different bodies, matters, and so on)
     */
    int type;

    /**
     * Parameters for Lennard-Jones force calculation
     */
    double sigma, epsilon;


public:
    explicit Particle(int type = 0);

    //Particle(const Particle &other);

    Particle(
            // for visualization, we need always 3 coordinates
            // -> in case of 2d, we use only the first and the second
            std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
            int type = 0, double sigma_arg = 1, double epsilon_arg = 5);

    double getSigma() const;

    double getEpsilon() const;

    virtual ~Particle();

    const std::array<double, 3> &getX() const;

    void setX(std::array<double, 3> xNew);

    const std::array<double, 3> &getV() const;

    void setV(std::array<double, 3> vNew);

    const std::array<double, 3> &getF() const;

    void setF(std::array<double, 3> fNew);

    const std::array<double, 3> &getOldF() const;

    void setOldF(std::array<double, 3> oldFNew);

    double getM() const;

    void setType(int typeNew);

    int getType() const;

    void mark();

    void unmark();

    bool isMarked();

    bool operator==(Particle &other);

    bool operator==(const Particle &rhs) const;

    bool operator==(const Particle *rhs) const;

    std::string toString() const;

    void setForceBuffer(const std::array<double, 3> &forceBuffer);

    const std::array<double, 3> &getForceBuffer() const;

    void addNeighbour(Particle &p, bool diag);
    std::vector<std::reference_wrapper<Particle>> getNeighbours();
    std::vector<std::reference_wrapper<Particle>> getNeighboursDiag();

    bool isMembranePull();
    void setMembranePull();
};

std::ostream &operator<<(std::ostream &stream, Particle &p);


