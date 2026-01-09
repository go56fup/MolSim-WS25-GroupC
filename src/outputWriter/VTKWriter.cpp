/*
 * VTKWriter.cpp
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */
#ifdef ENABLE_VTK_OUTPUT

#include "outputWriter/VTKWriter.h"
#include "Particle.h"
#include "ParticleContainer.h"

#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <iomanip>
#include <sstream>
#include <string_view>

namespace outputWriter {

void VTKWriter::plotParticles(ParticleContainer& particles, std::string_view filename, unsigned iteration) {
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

	const vtkNew<vtkIntArray> typeArray;
	typeArray->SetName("type");
	typeArray->SetNumberOfComponents(1);

	for (const auto& cell : particles) {
		for (const auto& p : cell) {
			points->InsertNextPoint(p.x.data());
			massArray->InsertNextValue(static_cast<float>(p.m));
			velocityArray->InsertNextTuple(p.v.data());
			forceArray->InsertNextTuple(p.f.data());
			typeArray->InsertNextValue(p.type);
		}
	}

	// Set up the grid
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(points);

	// Add arrays to the grid
	grid->GetPointData()->AddArray(massArray);
	grid->GetPointData()->AddArray(velocityArray);
	grid->GetPointData()->AddArray(forceArray);
	grid->GetPointData()->AddArray(typeArray);

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
}  // namespace outputWriter
#endif
