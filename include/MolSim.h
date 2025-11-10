#pragma once

#include <cmath>
#include <functional>
#include <span>
#include <string_view>

#include "ForceCalculators.h"
#include "IOProviders.h"
#include "outputWriter/outputWriters.h"
#include "ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"

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
		const auto force_scalar = std::pow(delta_t, 2) / (2 * p.m);
		const auto new_x = p.x + delta_t * p.v + p.old_f * force_scalar;
		p.x = new_x;
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
		const auto velocity_scalar = delta_t / (2 * p.m);
		const auto new_v = p.v + velocity_scalar * (p.old_f + p.f);
		p.v = new_v;
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
inline void plotParticles(
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
		p.old_f = p.f;
	}
}

// TODO(anyone): document me
constexpr void run_simulation(std::span<Particle> particles, double delta_t, double end_time) noexcept {
	double current_time = 0;
	int iteration = 0;

	while (current_time < end_time) {
		calculateX(particles, delta_t);
		calculateF(gravitational_force, particles);
		calculateV(particles, delta_t);

		iteration++;
		if not consteval {
			// NOLINTNEXTLINE(*magic-numbers)
			if (iteration % 10 == 0) {
				// TODO(tuna): Stop spewing VTK files into cwd for worksheets after 1.
				plotParticles(OUTPUT_WRITER::plotParticles, particles, "MD_vtk", iteration);
			}
		}
		update_values(particles);
		current_time += delta_t;
	}
}

// TODO(anyone): document me
constexpr ParticleContainer cuboid_generator(const vec& origin, const vec& scale, int distance, int mass, double mean) {
	ParticleContainer particles;

	for (int i = 0; i < scale.x; ++i) {
		for (int j = 0; j < scale.y; ++j) {
			for (int k = 0; k < scale.z; ++k) {
				particles.emplace_back(
					vec{origin.x + (i * distance), origin.y + (j * distance), origin.z + (k * distance)},
					maxwellBoltzmannDistributedVelocity<3>(mean), mass
				);
			}
		}
	}

	return particles;
}
