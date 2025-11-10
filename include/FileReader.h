/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <string_view>

#include "ParticleContainer.h"

namespace FileReader {
void readFile(ParticleContainer& particles, std::string_view filename);
}
