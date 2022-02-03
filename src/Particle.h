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
#include <omp.h>
#include <vector>
#include <memory>
#include <functional>

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

    /**
     * Buffer for multiple threads to store their force calculations without blocking
     */
    std::vector<std::array<double, 3>> parallelForceBuffer;

#ifdef _OPENMP
    /**
     * Lock for access to force array of the particle
     */
    omp_lock_t lck;

    /**
    * Lock for access to mark of the particle
    */
    omp_lock_t markLck;
#endif
    /**
     * Buffer to circumvent forces being moved to oldF when calculating forces outside of calculateF
     */
    std::array<double, 3> forceBuffer;

    /**
     * Stores all direct neighbours in the xy plane (only for membrane)
     */
    //std::vector<std::reference_wrapper<Particle>> neighbours;
    std::vector<int> neighbours;

    /**
     * Stores all diagonal neighbours in the xy plane (only for membrane)
     */
    //std::vector<std::reference_wrapper<Particle>> neighboursDiag;
    std::vector<int> neighboursDiag;

    /*
     * ID for a particle
     */
    int id;


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
     * Id of the particle
     */
    int id;

    /**
     * Parameters for Lennard-Jones force calculation
     */
    double sigma, epsilon;


public:

    //std::vector<int> &getCalculated();

public:
    explicit Particle(int type = 0);

    //Particle(const Particle &other);

    Particle(
            // for visualization, we need always 3 coordinates
            // -> in case of 2d, we use only the first and the second
            std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
            int type = 0, double sigma_arg = 1, double epsilon_arg = 5, int id_arg = -1);

    //bool calculatedContains(int i);

    double getSigma() const;

    double getEpsilon() const;

    virtual ~Particle();

    const std::array<double, 3> &getX() const;

    void setX(std::array<double, 3> xNew);

    const std::array<double, 3> &getV() const;

    void setV(std::array<double, 3> vNew);

     std::array<double, 3> &getF() ;

    //void setCalculated(const std::vector<int> &calculated);

    void setF(std::array<double, 3> fNew);

     std::array<double, 3> &getOldF() ;

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

#ifdef _OPENMP
    void setLock();

    void unSetLock();

    int testLock();

    void setMarkLock();

    void unSetMarkLock();

    int testMarkLock();
#endif
    void initParallelBuffer(int numThreads);

    std::array<double, 3> &getParallelForce(int index);

    void setParallelForce(int index, std::array<double, 3> force);

    /**
     * Adds a particle to the neighbour list
     *
     * @param p Particle to add
     * @param diag Should be true iff the particle is a diagonal neighbour, false if not
     */
    void addNeighbour(int id, bool diag);
    //void addNeighbour(std::reference_wrapper<Particle> p, bool diag);

    /**
     *
     * @return All direct neighboured particles
     */
    //std::vector<std::reference_wrapper<Particle>>& getNeighbours();
    //Particle& getNeighbours();
    std::vector<int> getNeighbours();
    /**
     *
     * @return All diagonal neighboured particles
     */
    //std::vector<std::reference_wrapper<Particle>>& getNeighboursDiag();
    //Particle& getNeighboursDiag();
    std::vector<int> getNeighboursDiag();

    /**
     *
     * @return True if the particle should be pulled up in the membrane simulation
     */
    bool isMembranePull();
    /**
     * Sets the membranePull to true
     */
    void setMembranePull();

    void setID(int id_arg);
    int getID();

};

std::ostream &operator<<(std::ostream &stream, Particle &p);


