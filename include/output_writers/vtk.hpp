/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once
#ifdef ENABLE_VTK_OUTPUT
#include <iomanip>
#include <sstream>
#include <string_view>

#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "grid/particle_container/fwd.hpp"

namespace output_writer::vtk {
/**
 * Write VTK output of particles.
 * @param particles Particles to add to the output
 * @param filename Output filename
 * @param iteration Current iteration number
 */
inline void
plot_particles(particle_container& particles, std::string_view filename, unsigned iteration) {
	// Initialize points
	auto points = vtkSmartPointer<vtkPoints>::New();

	// Create and configure data arrays
	const vtkNew<vtkFloatArray> massArray;
	massArray->SetName("mass");
	massArray->SetNumberOfComponents(1);

	const vtkNew<vtkFloatArray> velocityArray;
	velocityArray->SetName("velocity");
	velocityArray->SetNumberOfComponents(3);

	const vtkNew<vtkFloatArray> forceArray;
	forceArray->SetName("force");
	forceArray->SetNumberOfComponents(3);

	const auto& system = particles.system();
	for (particle_id i = 0; i < system.size(); ++i) {
		points->InsertNextPoint(system.serialize_position(i).data());
		massArray->InsertNextValue(static_cast<float>(system.mass[i]));
		velocityArray->InsertNextTuple(system.serialize_velocity(i).data());
		forceArray->InsertNextTuple(system.serialize_force(i).data());
	}

	// Set up the grid
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(points);

	// Add arrays to the grid
	grid->GetPointData()->AddArray(massArray);
	grid->GetPointData()->AddArray(velocityArray);
	grid->GetPointData()->AddArray(forceArray);

	// Create filename with iteration number
	std::stringstream strstr;
	strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".vtu";

	// Create writer and set data
	const vtkNew<vtkXMLUnstructuredGridWriter> writer;
	writer->SetFileName(strstr.str().c_str());
	writer->SetInputData(grid);
	writer->SetDataModeToAscii();

	// Write the file
	writer->Write();
}
}  // namespace output_writer::vtk
#endif
