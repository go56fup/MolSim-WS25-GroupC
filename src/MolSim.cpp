#if __cpp_lib_to_chars == 201611L
#include <charconv>
#endif
#include <cmath>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <span>
#include <string>
#include <string_view>

#include "FileReader.h"
#include "ForceCalculators.h"
#include "IOProviders.h"
#include "outputWriter/VTKWriter.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include "Vector.h"

/**
 * @brief Calculates forces on a collection of particles using the provided force calculator.
 *
 * This function delegates the force calculation to the provided @p calculator callable, allowing
 * different calculation methods.
 *
 * @param calculator A callable object that computes forces for a span of particles.
 * Its type must satisfy force_calculator.
 * @param particles The span over particles on which forces will be calculated.
 *
 */
constexpr void calculateF(force_calculator auto&& calculator, std::span<Particle> particles) noexcept(
	noexcept(std::invoke(calculator, particles))
) {
	std::invoke(calculator, particles);
}

/**
 * @brief Updates the position of each particle.
 *
 * Implements the first step of the velocity Störmer-Verlet algorithm:
 * \f[
 * \vec{x}(t + \Delta t) = \vec{x}(t) + \Delta t \cdot \vec{v}(t) + \frac{\Delta t^2}{2m} \vec{F}(t)
 * \f]
 *
 * @param particles Mutable span over particles to update.
 * @param delta_t The time step for integration.
 */
constexpr void calculateX(std::span<Particle> particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto force_scalar = std::pow(delta_t, 2) / (2 * p.getM());
		const auto new_x = p.getX() + delta_t * p.getV() + p.getOldF() * force_scalar;
		p.setX(new_x);
	}
}

/**
 * @brief Updates the velocity of each particle.
 *
 * Implements the second step of the velocity Störmer-Verlet algorithm:
 * \f[
 * \vec{v}(t + \Delta t) = \vec{v}(t) + \frac{\Delta t}{2m} \left[\vec{F}(t) + \vec{F}(t + \Delta t)\right]
 * \f]
 *
 * @param particles Mutable span over particles to update.
 * @param delta_t The time step for integration.
 */
constexpr void calculateV(std::span<Particle> particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto velocity_scalar = delta_t / (2 * p.getM());
		const auto new_v = p.getV() + velocity_scalar * (p.getOldF() + p.getF());
		p.setV(new_v);
	}
}

/**
 * @brief Exports particle data using the given I/O provider.
 *
 * This function delegates to the given @p io_provider callable, allowing different
 * data exporting strategies.
 *
 * @param io_provider Callable object responsible for exporting particle data. Its type must satisfy
 * particle_io_provider.
 * @param particles Constant span of particles to export data from.
 * @param out_name  Base name for the output file.
 * @param iteration Current simulation iteration (used for output naming).
 */
void plotParticles(
	particle_io_provider auto&& io_provider, std::span<const Particle> particles, std::string_view out_name,
	int iteration
) {
	io_provider(particles, out_name, iteration);
}

/**
 * @brief Prepares particles for the next iteration.
 *
 * Currently only sets the old force of each particle to the force calculated within the
 * most recent iteration at the end of each iteraton.
 *
 * @param particles Mutable span over particles to update.
 */
constexpr void update_values(std::span<Particle> particles) noexcept {
	for (auto& p : particles) {
		p.setOldF(p.getF());
	}
}

/**
 * @brief Prints an error message and terminates the program.
 *
 * Used when command-line arguments are invalid or missing.
 *
 * @note This function never returns.
 */
[[noreturn]] void die() {
	std::cout << "Erroneous programme call!\n./molsym filename\n";
	std::exit(1);  // NOLINT(*mt-unsafe)
}

/**
 * @brief Parses an argument from argv into a double value.
 *
 * Falls back to `std::stod` if `std::from_chars` is unavailable.
 * Exits the program if parsing fails.
 *
 * @param str View of the passed command-line argument.
 * @return The parsed double value.
 */
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

/**
 * @brief Returns a simulation parameter from either the default or a command-line argument.
 *
 * @param use_default Whether to use the default value.
 * @param input Pointer to a command-line argument to parse (ignored if @p use_default is `true`).
 * @param default_ The default value to use.
 * @return The parsed or default parameter.
 */
constexpr double get_sim_parameter(bool use_default, const char* input, double default_) {
	if (use_default) return default_;
	return get_double_from_argv(input);
}

/**
 * @brief Entry point of our simulation program.
 *
 * The program performs a basic molecular dynamics simulation:
 *  - Reads initial particle data from a file.
 *  - Calculates particle states over time.
 *  - Periodically writes simulation snapshots to VTK files.
 *
 * Command-line usage:
 * ```
 * ./MolSim input_file [delta_t end_time]
 * ```
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit status code (0 on success).
 */
int main(int argc, char* argv[]) {
	static constexpr int argc_with_everything_defaulted = 2;
	static constexpr int argc_with_everything_explicit = 5;
	std::cout << "Hello from MolSim for PSE!\n";
	if (argc < argc_with_everything_defaulted || argc > argc_with_everything_explicit) die();
	const bool defaulted = argc == 2;
	// NOLINTBEGIN(*pointer-arithmetic)
	const double delta_t = get_sim_parameter(defaulted, argv[2], 0.014);
	const double end_time = get_sim_parameter(defaulted, argv[3], 1000);

	ParticleContainer particles;
	FileReader::readFile(particles, argv[1]);
	// NOLINTEND(*pointer-arithmetic)

	double current_time = 0;
	int iteration = 0;

	// for this loop, we assume: current x, current f and current v are known
	while (current_time < end_time) {
		calculateF(gravitational_force, particles);
		calculateX(particles, delta_t);
		calculateV(particles, delta_t);

		iteration++;
		// NOLINTNEXTLINE(*magic-numbers)
		if (iteration % 10 == 0) {
			plotParticles(outputWriter::VTKWriter::plotParticles, particles, "MD_vtk", iteration);
		}
		std::cout << "Iteration " << iteration << " finished.\n";
		update_values(particles);
		current_time += delta_t;
	}

	std::cout << "output written. Terminating...\n";
	return 0;
}
