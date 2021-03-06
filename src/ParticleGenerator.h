//
// Created by paulu on 10.11.2021.
//

#pragma once

#include "array"
#include "container/SimulationContainer.h"


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
void generateCube(std::array<int, 3> dimension, std::array<double, 3> startPoint, double h, double m,
                  std::array<double, 3> v, double meanV, SimulationContainer &container);

/**
 * Generates a cube of particles with given parameters. Those particles will be added to the given container. Additionally sets type of each particle
 *
 * @param dimension Number of particles per dimension: N1 × N2 × N3
 * @param startPoint Coordinates of the lower left front-side corner: X × Y × Z
 * @param h Distance h of the particles (mesh width of the grid)
 * @param m Mass m of one particle
 * @param v Initial velocity v of the particles
 * @param meanV The mean-value of the velocity of the Brownian Motion
 * @param type type of the particles
 * @param container
 */
void generateCube(std::array<int, 3> dimension, std::array<double, 3> startPoint, double h, double m,
                  std::array<double, 3> v, double meanV, int type, double sigma, double epsilon,
                  SimulationContainer &container);

/**
 * Generates a sphere of particles with given parameters. Those particles will be added to the given container.
 * @param center center point of the sphere
 * @param v Initial velocity v of the particles
 * @param r Radius of the sphere (given in the amount of particles)
 * @param h Distance h of the particles (mesh width of the grid)
 * @param m Mass m of one particle
 * @param meanV The mean-value of the velocity of the Brownian Motion
 * @param container
 */
void generateSphere3D(std::array<double, 3> center, std::array<double, 3> v, int r, double h, double m, double meanV,
                      SimulationContainer &container);

/**
 * Generates a 2D sphere (so a circle) of particles with given parameters. Those particles will be added to the given container.
 *
 * @param center center point of the sphere
 * @param v Initial velocity v of the particles
 * @param r Radius of the sphere (given in the amount of particles)
 * @param h Distance h of the particles (mesh width of the grid)
 * @param m Mass m of one particle
 * @param meanV The mean-value of the velocity of the Brownian Motion
 * @param container
 */
void generateSphere2D(std::array<double, 3> center, std::array<double, 3> v, int r, double h, double m, double meanV,
                      SimulationContainer &container);

/**
 * Generates a 2D sphere (so a circle) of particles with given parameters. Those particles will be added to the given container.
 *
 * @param center center point of the sphere
 * @param v Initial velocity v of the particles
 * @param r Radius of the sphere (given in the amount of particles)
 * @param h Distance h of the particles (mesh width of the grid)
 * @param m Mass m of one particle
 * @param meanV The mean-value of the velocity of the Brownian Motion
 * @param type type of the particles
 * @param container
 */
void generateSphere2D(std::array<double, 3> center, std::array<double, 3> v, int r, double h, double m, double meanV,
                      int type, SimulationContainer &container);
/**
 * Generates a membrane
 *
 * @param center center point of the membrane
 * @param v Initial velocity of all particles
 * @param dimension Dimensions in x,y,z
 * @param h Distance h of the particles
 * @param pull Coordinates of particles which will be pulled
 * @param m Mass
 * @param container SiimulationContainer where the particles will be stored in (LinkedCells!)
 */
void generateMembrane(std::array<double, 3> center, std::array<double, 3> v,std::array<int, 3> dimension, double h, const std::vector<std::array<int, 3>>& pull, double m, SimulationContainer &container);

/**
 * Generates cubes from input file
 * Format:
 *  numberOfConstructs
 *  Cube:
 *  x,y,z;v1,v2,v3;N1,N2,N3;h;meanVelocity;mass;
 *
 * This can be repeated with multiple cubes e.g. (without leading spaces)
    2
    Cube:
    0,0,0;0,0,0;40,8,1;1.1225;0.1;1;
    Cube:
    15,15,0;0,-10,0;8,8,1;1.1225;0.1;1;

 *
 *Instead of Cube you can use Sphere for 2D spheres like this:
 *
 * Sphere:
 * x,y,z;v1,v2,v3;r;h;meanVelocity;mass;
 * note that r (radius in particle count) is an integer
 *
 * @param particles Particle container
 * @param filename Path to input file
 */
void generateFromFile(SimulationContainer &particles, char *filename);

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
std::array<double, 3> convertToFixedArray(std::vector<double> v);

/**
 * Parses and generates a sphere with given input (see generateFromFile())
 *
 * @param str Input string (format see: generateFromFile())
 * @param particleContainer
 */
void parseSphere2D(std::string str, SimulationContainer &particleContainer);

/**
 * Parses and generates a cube with given input (see generateFromFile())
 *
 * @param str Input string (format see: generateFromFile())
 * @param particleContainer
 */
void parseCube(std::string str, SimulationContainer &particleContainer);

/**
 * see generateFromFile(SimulationContainer &particles, char *filename)
 *
 * @param particles
 * @param f file content
 */
void generateFromFileTest(SimulationContainer &particles, std::string f);

