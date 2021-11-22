//
// Created by paulu on 21.11.2021.
//

#include "LinkedCells.h"
#include <iostream>


// TODO 2D Grid
LinkedCells::LinkedCells(std::array<int, 3> dimension, double mesh) {
    size = coordToIndex({dimensions[0], dimensions[1], dimensions[2]});
    // Add boundary and halo cells
    dimensions = {dimension[0] + 2, dimension[1] + 2, dimension[2] + 2};
    meshSize = mesh;
    // Initialise all vectors
    for(int i = 0; i < coordToIndex({dimensions[0], dimensions[1], dimensions[2]}); i++) {
        particles.push_back(std::vector<Particle>());
    }
}

int LinkedCells::coordToIndex(std::array<int, 3> coords) {
    return coords[0] * dimensions[1] * dimensions[2] + coords[1] * dimensions[2] + coords[2];
}

std::array<int, 3> LinkedCells::indexToCoords(int index) {
    int width_index = index / (dimensions[1] * dimensions[2]);
    int height_index = (index - width_index * dimensions[1] * dimensions[2]) / dimensions[2];
    int depth_index = index - width_index * dimensions[1] * dimensions[2] - height_index * dimensions[2];
    return {width_index, height_index, depth_index};
}

void LinkedCells::test() {
    Particle* p = new Particle({2.2,1.2,1.2},{1,0,0},1,0);
    insert(*p);
    Particle* p1 = new Particle({2.3,1.2,1.2},{1,0,0},1,0);
    insert(*p1);
    std::cout << "### TEST: " << coordToIndex(getCellCoords(*p));
    std::cout << "### TEST: " << getCellCoords(*p)[0] << getCellCoords(*p)[1] << getCellCoords(*p)[2] << std::endl;
    remove(*p);
    std::cout << "TEST";
}

void LinkedCells::insert(Particle &p) {
    particles[coordToIndex(getCellCoords(p))].push_back(p);
}

std::array<int, 3> LinkedCells::getCellCoords(Particle &p) {
    int x = ((int) (p.getX()[0] / meshSize));
    int y = ((int) (p.getX()[1] / meshSize) );
    int z = ((int) (p.getX()[2] / meshSize) );
    return {x,y,z};
}

void LinkedCells::remove(Particle& p) {
    int currentIndex = 0;
    // Search for right element
    for (auto p1 : particles[coordToIndex(getCellCoords(p))]) {
        if (p1 == p) {
            break;
        }
        currentIndex++;
    }
    // Remove particle by index
    particles[coordToIndex(getCellCoords(p))].erase(particles[coordToIndex(getCellCoords(p))].begin() + currentIndex);
}

bool LinkedCells::isHaloCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0]-1 || coords[1] == dimensions[1]-1 || coords[2] == dimensions[2]-1;
}

bool LinkedCells::isBoundaryCell(std::array<int, 3> coords) {
    return coords[0] == dimensions[0]-2 || coords[1] == dimensions[1]-2 || coords[2] == dimensions[2]-2;
}

std::vector<Particle>& LinkedCells::get(std::array<double, 3> coords) {
    int x = ((int) (coords[0] / meshSize));
    int y = ((int) (coords[1] / meshSize) );
    int z = ((int) (coords[2] / meshSize) );
    return particles[coordToIndex({x,y,z})];
}

void LinkedCells::calculateF(ForceCalculation *algorithm) {

}


