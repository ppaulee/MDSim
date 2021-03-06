/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "Particle.h"
#include "container/SimulationContainer.h"
#include <vector>

class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  void readFile(SimulationContainer &particles, char *filename);
};
