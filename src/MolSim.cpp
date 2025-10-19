#include <cmath>
#include <iostream>
#include <vector>

#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"
#include "Vector.h"

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

std::vector<Particle> particles;

constexpr void update_values(std::span<Particle> ps) noexcept {
	for (auto& p : ps) {
		p.setOldF(p.getF());
	}
}

#include <iostream>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

int main(int argc, char* argsv[]) {
	std::cout << "Hello from MolSim for PSE!" << std::endl;
	if (argc != 2) {
		std::cout << "Erroneous programme call! " << std::endl;
		std::cout << "./molsym filename" << std::endl;
	}

	FileReader fileReader;
	fileReader.readFile(particles, argsv[1]);

	double current_time = start_time;

	int iteration = 0;

	// for this loop, we assume: current x, current f and current v are known
	while (current_time < end_time) {
		// calculate new x
		calculateX();
		// calculate new f
		calculateF();
		// calculate new v
		calculateV();

		iteration++;
		if (iteration % 10 == 0) {
			plotParticles(iteration);
		}
		std::cout << "Iteration " << iteration << " finished." << std::endl;
		update_values(particles);
		current_time += delta_t;
	}

	std::cout << "output written. Terminating..." << std::endl;
	return 0;
}

double cube(double x) {
	return x * x * x;
}

vec calculate_component(const Particle& p1, const Particle& p2) {
	const auto xi = p1.getX();
	const auto xj = p2.getX();
	const double reciprocal = cube((xi - xj).euclidian_norm());
	const double scaling_factor = p1.getM() * p2.getM() / reciprocal;
	return scaling_factor * (xj - xi);
}

void calculateF() {
	for (auto& p1 : particles) {
		vec new_f;
		for (auto& p2 : particles) {
			if (p1 == p2) continue;
			const auto f_ij = calculate_component(p1, p2);
			new_f += f_ij;
		}
		p1.setF(new_f);
	}
}

void calculateX() {
	for (auto& p : particles) {
		const auto force_scalar = std::pow(delta_t, 2) / (2 * p.getM());
		const auto new_x = p.getX() + delta_t * p.getV() + p.getF() * force_scalar;
		p.setX(new_x);
	}
}

void calculateV() {
	for (auto& p : particles) {
		const auto velocity_scalar = delta_t / (2 * p.getM());
		const auto new_v = p.getV() + velocity_scalar * (p.getOldF() + p.getF());
		p.setV(new_v);
	}
}

void plotParticles(int iteration) {
	std::string_view out_name = "MD_vtk";

	outputWriter::VTKWriter writer;
	writer.plotParticles(particles, out_name, iteration);
}
