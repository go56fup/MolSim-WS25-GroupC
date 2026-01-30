/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include <fstream>
#include <span>
#include <string_view>

#include "physics/particle.hpp"

namespace output_writer::xyz {
void plotParticles(
	std::span<const Particle> particles, std::string_view filename, unsigned iteration
);
}
