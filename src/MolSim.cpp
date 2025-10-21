#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iostream>
#include <vector>
#include <ranges>
#include <charconv>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>

#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"
#include "Vector.h"


/*
 *	@param p1 represents the particle with coordinates xi
 	@param p2 represents the particle with coordinates xj
	This function calculates the force acting on i from j
 */
constexpr vec calculate_component(const Particle& p1, const Particle& p2) noexcept {
	const auto xi = p1.getX();
	const auto xj = p2.getX();
	const double reciprocal = std::pow((xi - xj).euclidian_norm(), 3);
	const double scaling_factor = p1.getM() * p2.getM() / reciprocal;
	return scaling_factor * (xj - xi);
}


/// This function sums over all the component forces acting on each particle.
constexpr void calculateF(std::span<Particle> particles) noexcept {
	for (std::size_t i = 0; i < particles.size(); ++i) {
		auto& p1 = particles[i];
		vec new_f{};
		for (std::size_t j = 0; j < particles.size(); ++j) {
			const auto& p2 = particles[j];
			if (i == j) continue;
			const auto f_ij = calculate_component(p1, p2);
			new_f += f_ij;
		}
		p1.setF(new_f);
	}
}

/// This function updates the coordinate of each particle
constexpr void calculateX(std::span<Particle> particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto force_scalar = std::pow(delta_t, 2) / (2 * p.getM());
		const auto new_x = p.getX() + delta_t * p.getV() + p.getOldF() * force_scalar;
		p.setX(new_x);
	}
}

/// This function updates the velocity of each particle
constexpr void calculateV(std::span<Particle> particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto velocity_scalar = delta_t / (2 * p.getM());
		const auto new_v = p.getV() + velocity_scalar * (p.getOldF() + p.getF());
		p.setV(new_v);
	}
}

void plotParticles(std::span<Particle> particles, int iteration) {
	std::string_view out_name = "MD_vtk";

	outputWriter::VTKWriter writer;
	writer.plotParticles(particles, out_name, iteration);
}

/// This function updates all values for the start of the next tick (currently only the old Force)
constexpr void update_values(std::span<Particle> ps) noexcept {
	for (auto& p : ps) {
		p.setOldF(p.getF());
	}
}

[[noreturn]] void die() {
	std::cout << "Erroneous programme call! " << std::endl;
	std::cout << "./molsym filename" << std::endl;
	std::exit(1);
}

constexpr double get_double_from_argv(std::string_view str) {
	auto invalid_input = []{
		std::cout << "Given decimal number is invalid.\n";
		die();
	};
	double result = 0;
#if __cpp_lib_to_chars == 201611L
	const auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.length(), result);
	if (*ptr != '\0' || ec == std::errc::invalid_argument || ec == std::errc::result_out_of_range) {
		invalid_input();
	}
#else
	std::size_t pos = 0;
	try {
		result = std::stod(str.data(), &pos);
	} catch (std::exception&) {
		invalid_input();
	}
	if (pos != str.length()) invalid_input();
#endif
	return result;
}

struct get_parameters_args {
	bool use_default;
	const char* input;
	double default_;
};

constexpr double get_parameters(const get_parameters_args& args) {
	if (args.use_default) return args.default_;
	return get_double_from_argv(args.input);
}

int main(int argc, char* argv[]) {
	std::cout << "Hello from MolSim for PSE!" << std::endl;
	if (argc < 2 || argc > 5) die();
	const bool defaulted = argc == 2;
	const double delta_t = get_parameters({.use_default = defaulted, .input = argv[2], .default_ = 0.014});
	const double end_time = get_parameters({.use_default = defaulted, .input = argv[3], .default_ = 1000});

	FileReader fileReader;
	std::vector<Particle> particles;
	fileReader.readFile(particles, argv[1]);

	double current_time = 0;
	int iteration = 0;

	// for this loop, we assume: current x, current f and current v are known
	while (current_time < end_time) {
		calculateF(particles);
		calculateX(particles, delta_t);
		calculateV(particles, delta_t);

		iteration++;
		if (iteration % 10 == 0) {
			plotParticles(particles, iteration);
		}
		std::cout << "Iteration " << iteration << " finished." << std::endl;
		update_values(particles);
		current_time += delta_t;
	}

	std::cout << "output written. Terminating..." << std::endl;
	return 0;
}
