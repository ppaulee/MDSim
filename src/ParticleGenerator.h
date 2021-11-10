//
// Created by paulu on 10.11.2021.
//

#pragma once
#include "array"
#include "ParticleContainer.h"


/**
 * Generates a cube of particles with given parameters. Those particles will be added to the given container.
 *
 * @param dimension Number of particles per dimension: N1 × N2 × N3
 * @param startPoint Coordinates of the lower left front-side corner: X × Y × Z
 * @param h Distance h of the particles (mesh width of the grid)
 * @param m Mass m of one particle
 * @param v Initial velocity v of the particles
 * @param meanV The mean-value of the velocity of the Brownian Motion
 * @param container
 */
void generateCube(std::array<double, 3> dimension, std::array<double, 3> startPoint, double h, double m, std::array<double, 3> v, double meanV, ParticleContainer &container);

/**
 * Generates cubes from input file
 * Format:
 *  numberOfConstructs
 *  Cube:
 *  x,y,z;v1,v2,v3;N1,N2,N3;h;meanVelocity;mass;
 *
 * This can be repeated with multiple cubes e.g.
 *  2
 *  Cube:
 *  0,0,0;0,0,0;40,8,1;1.1225;0.1;1;
 *  Cube:
 *  15,15,0;0,-10,0;8,8,1;1.1225;0.1;1;
 *
 *
 *
 * @param particles Particle container
 * @param filename Path to input file
 */
void generateFromFile(ParticleContainer &particles, char *filename);

/**
 * Splits a string into a vector of string
 *
 * @param str
 * @param delimiter
 * @return vector of string
 */
std::vector<std::string> splitToString(std::string str, char delimiter);

/**
 * Splits a string into a vector of doubles
 *
 * @param str
 * @param delimiter
 * @return vector of doubles
 */
std::vector<double> splitToDouble(std::string str, char delimiter);

/**
 * Converts a vector of size 3 to an array of size 3
 *
 * @param v vector of size 3
 * @return array of size 3
 */
std::array<double,3> convertToFixedArray(std::vector<double> v);

/**
 * Parses and generates a cube with given input (see generateFromFile())
 *
 * @param str Input string (format see: generateFromFile())
 * @param particleContainer
 */
void parseCube(std::string str, ParticleContainer& particleContainer);

