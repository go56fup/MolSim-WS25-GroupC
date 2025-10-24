#if __cpp_lib_to_chars == 201611L
#include <charconv>
#endif
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <span>
#include <string>
#include <string_view>

#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include "Vector.h"

constexpr vec calculate_component(const Particle& p1, const Particle& p2) noexcept {
	const auto xi = p1.getX();
	const auto xj = p2.getX();
	const double reciprocal = std::pow((xi - xj).euclidian_norm(), 3);
	const double scaling_factor = p1.getM() * p2.getM() / reciprocal;
	return scaling_factor * (xj - xi);
}

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

constexpr void calculateX(std::span<Particle> particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto force_scalar = std::pow(delta_t, 2) / (2 * p.getM());
		const auto new_x = p.getX() + delta_t * p.getV() + p.getOldF() * force_scalar;
		p.setX(new_x);
	}
}

constexpr void calculateV(std::span<Particle> particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto velocity_scalar = delta_t / (2 * p.getM());
		const auto new_v = p.getV() + velocity_scalar * (p.getOldF() + p.getF());
		p.setV(new_v);
	}
}

void plotParticles(std::span<const Particle> particles, int iteration) {
	static constexpr std::string_view out_name = "MD_vtk";
	outputWriter::VTKWriter::plotParticles(particles, out_name, iteration);
}

constexpr void update_values(std::span<Particle> ps) noexcept {
	for (auto& p : ps) {
		p.setOldF(p.getF());
	}
}

[[noreturn]] void die() {
	std::cout << "Erroneous programme call!\n./molsym filename\n";
	std::exit(1);  // NOLINT(*mt-unsafe)
}

constexpr double get_double_from_argv(std::string_view str) {
	auto invalid_input = [] {
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
		// Because the string_view points to argv entries, which are guaranteed to be null-terminated,
		// this is fine.
		// NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
		result = std::stod(str.data(), &pos);
	} catch (std::exception&) {
		invalid_input();
	}
	if (pos != str.length()) invalid_input();
#endif
	return result;
}

constexpr double get_parameters(bool use_default, const char* input, double default_) {
	if (use_default) return default_;
	return get_double_from_argv(input);
}

int main(int argc, char* argv[]) {
	static constexpr int argc_with_everything_defaulted = 2;
	static constexpr int argc_with_everything_explicit = 5;
	std::cout << "Hello from MolSim for PSE!\n";
	if (argc < argc_with_everything_defaulted || argc > argc_with_everything_explicit) die();
	const bool defaulted = argc == 2;
	// NOLINTBEGIN(*pointer-arithmetic)
	const double delta_t = get_parameters(defaulted, argv[2], 0.014);
	const double end_time = get_parameters(defaulted, argv[3], 1000);

	ParticleContainer particles;
	FileReader::readFile(particles, argv[1]);
	// NOLINTEND(*pointer-arithmetic)

	double current_time = 0;
	int iteration = 0;

	// for this loop, we assume: current x, current f and current v are known
	while (current_time < end_time) {
		calculateF(particles);
		calculateX(particles, delta_t);
		calculateV(particles, delta_t);

		iteration++;
		// NOLINTNEXTLINE(*magic-numbers)
		if (iteration % 10 == 0) {
			plotParticles(particles, iteration);
		}
		std::cout << "Iteration " << iteration << " finished.\n";
		update_values(particles);
		current_time += delta_t;
	}

	std::cout << "output written. Terminating...\n";
	return 0;
}
