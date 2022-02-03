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
    int getId();

    //std::vector<int> &getCalculated();

public:
    explicit Particle(int type = 0);

    //Particle(const Particle &other);

    Particle(
            // for visualization, we need always 3 coordinates
            // -> in case of 2d, we use only the first and the second
            std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
            int type = 0, double sigma_arg = 1, double epsilon_arg = 5, int id_arg = -1);

    void setId(int id);

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
};

std::ostream &operator<<(std::ostream &stream, Particle &p);
