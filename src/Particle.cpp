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
  std::cout << "Particle generated!" << std::endl;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  f = other.f;
  old_f = other.old_f;
  m = other.m;
  type = other.type;
  std::cout << "Particle generated by copy!" << std::endl;
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, int type_arg) {
  x = x_arg;
  v = v_arg;
  m = m_arg;
  type = type_arg;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
  std::cout << "Particle generated!" << std::endl;
}

Particle::~Particle() { std::cout << "Particle destructed!" << std::endl; }

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }



double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f
         << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

bool Particle::operator==(Particle &other) {
  return (x == other.x) and (v == other.v) and (f == other.f) and
         (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

void Particle::setX(std::array<double, 3> xNew) {
    x = xNew;
}

void Particle::setV(std::array<double, 3>  vNew) {
    v = vNew;
}

void Particle::setF(std::array<double, 3>  fNew) {
    f = fNew;
}

void Particle::setOldF(std::array<double, 3> oldFNew) {
    old_f = oldFNew;
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
  stream << p.toString();
  return stream;
}
