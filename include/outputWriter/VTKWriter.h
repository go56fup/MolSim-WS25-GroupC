/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once
#ifdef ENABLE_VTK_OUTPUT

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <span>
#include <string_view>

#include "Concepts.h"
#include "Particle.h"
#include "ParticleContainer.h"

namespace outputWriter::VTKWriter {
/**
 * Write VTK output of particles.
 * @param particles Particles to add to the output
 * @param filename Output filename
 * @param iteration Current iteration number
 */
void plotParticles(ParticleContainer& particles, std::string_view filename, unsigned iteration);
}  // namespace outputWriter::VTKWriter
#endif
