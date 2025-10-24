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

#include "Particle.h"

namespace outputWriter {

class XYZWriter {
public:
	XYZWriter();

	virtual ~XYZWriter();

	static void plotParticles(std::span<const Particle> particles, std::string_view filename, int iteration);
};

}  // namespace outputWriter
