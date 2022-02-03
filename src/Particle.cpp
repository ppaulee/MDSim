/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <iostream>
#include "utils/ArrayUtils.h"

Particle::Particle(int type_arg) {
    type = type_arg;
    //std::cout << "Particle generated!" << std::endl;
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    membranePull = false;

    //neighbours = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
}

/*Particle::Particle(const Particle &other) {
    x = other.x;
    v = other.v;
    f = other.f;
    marked = false;
    old_f = other.old_f;
    m = other.m;
    type = other.type;
    sigma = other.sigma;
    epsilon = other.epsilon;
    //std::cout << "Particle generated by copy!" << std::endl;
}*/

// Todo: maybe use initializater list instead of copy?
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, int type_arg, double sigma_arg, double epsilon_arg) {
    x = x_arg;
    v = v_arg;
    m = m_arg;
    marked = false;
    type = type_arg;
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    sigma = sigma_arg;
    epsilon = epsilon_arg;
    membranePull = false;
}

Particle::~Particle() {
    //std::cout << "Particle destructed!" << std::endl;
}

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }

void Particle::mark() { marked = true; }

void Particle::unmark() { marked = false; }

bool Particle::isMarked() { return marked; }


double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

std::string Particle::toString() const {
    std::stringstream stream;
    stream << "Particle: X:" << x << " v: " << v << " f: " << f << " m: " << m
           << " old_f: " << old_f << " type: " << type << " sigma: " << sigma << " epsilon: " << epsilon;
    return stream.str();
}

bool Particle::operator==(Particle &other) {
    return (x == other.x) and (v == other.v) and (f == other.f or isnan(f[0])) and
           (type == other.type) and (m == other.m) and (old_f == other.old_f) and (sigma == other.sigma) and
           (epsilon == other.epsilon);
}

bool Particle::operator==(const Particle &rhs) const {
    return (x == rhs.x);
}

bool Particle::operator==(const Particle *rhs) const {
    return (x == rhs->x);
}

void Particle::setX(std::array<double, 3> xNew) {
    x = xNew;
}

void Particle::setV(std::array<double, 3> vNew) {
    v = vNew;
}

void Particle::setF(std::array<double, 3> fNew) {
    f = fNew;
}

void Particle::setOldF(std::array<double, 3> oldFNew) {
    old_f = oldFNew;
}

void Particle::setType(int typeNew) {
    type = typeNew;
}

double Particle::getSigma() const {
    return sigma;
}

double Particle::getEpsilon() const {
    return epsilon;
}

void Particle::setForceBuffer(const std::array<double, 3> &forceBuffer) {
    Particle::forceBuffer = forceBuffer;
}

const std::array<double, 3> &Particle::getForceBuffer() const {
    return forceBuffer;
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
    stream << p.toString();
    return stream;
}

/*
void Particle::addNeighbour(std::shared_ptr<Particle> p, bool diag) {
    if (diag) {
        neighboursDiag.push_back(std::move(p));
    } else {
        neighbours.push_back(std::move(p));
    }
}
*/


void Particle::addNeighbour(int id, bool diag) {
    if (diag) {
        neighboursDiag.push_back(id);
    } else {
        neighbours.push_back(id);
    }
}


/*void Particle::addNeighbour(std::reference_wrapper<Particle> p, bool diag) {
    if (diag) {
        neighboursDiag.push_back(p);
    } else {
        neighbours.push_back(p);
    }
}*/


/*
std::vector<std::shared_ptr<Particle>>& Particle::getNeighbours() {
    return neighbours;
}

std::vector<std::shared_ptr<Particle>>& Particle::getNeighboursDiag() {
    return neighboursDiag;
}
 */




std::vector<int> Particle::getNeighbours() {
    return neighbours;
}

std::vector<int> Particle::getNeighboursDiag() {
    return neighboursDiag;
}


void Particle::setID(int id_arg) {id = id_arg;}
int Particle::getID() {return id;}

bool Particle::isMembranePull() {return membranePull;}
void Particle::setMembranePull() {membranePull = true;}
