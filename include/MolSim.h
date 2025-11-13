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
	noexcept(std::invoke(calculator, Particle{}, Particle{}))
) {
	// TODO(tuna): see if this actually compiles down to a mask over simd like it should
	for (auto& p : particles) {
		p.f = {};
	}

	for (auto&& [p1, p2] : unique_pairs(particles)) {
		const auto f_ij = std::invoke(calculator, p1, p2);
		p1.f += f_ij;
		p2.f -= f_ij;
	}
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

struct sim_args {
	double delta_t;
	double end_time;
};

/**
 * @brief Run a tick of the simulation.
 *
 * @param force_calc A callable object that computes forces for a span of particles.
 * Its type must satisfy force_calculator.
 * @param particles Mutable span over particles to run simulation tick for
 * @param delta_t Tick length
 */
template <force_calculator ForceCalculator>
constexpr void run_sim_iteration(ForceCalculator&& force_calc, std::span<Particle> particles, double delta_t) noexcept {
	calculateX(particles, delta_t);
	calculateF(std::forward<ForceCalculator>(force_calc), particles);
	calculateV(particles, delta_t);
	update_values(particles);
}

/** @brief Compile-time options for the simulation. */
struct sim_traits {
	/// Whether to call `plotParticles` and generate output.
	bool create_output = true;
};

/**
 * @brief Start the simulation.
 *
 * Calculates new positions, velocities and forces at each time tick
 * and plots the updated particles.
 *
 * @param particles Mutable span over particles to run simulation for
 * @param args Simulation parameters (@p delta_t and @p end_time)
 * @param output_path Path to put simulation result files into
 **/
template <sim_traits Traits = {}>
constexpr void
run_simulation(std::span<Particle> particles, const sim_args& args, std::string_view output_path) noexcept {
	double current_time = 0;
	int iteration = 0;
	const std::string output_prefix = std::string(output_path) + "/MD_vtk";

	while (current_time < args.end_time) {
		if constexpr (Traits.create_output) {
			// NOLINTNEXTLINE(*magic-numbers)
			if (iteration % 10 == 0) {
				plotParticles(OUTPUT_WRITER::plotParticles, particles, output_prefix, iteration);
			}
		}

		run_sim_iteration(lennard_jones_force, particles, args.delta_t);
		current_time += args.delta_t;
		++iteration;
	}
}

/**
 * @brief Create a 3D grid of particles, representing one body.
 *
 * @param particles Particle container to add generated particles onto
 * @param origin Position of the lower left front-side corner
 * @param initial_velocity Initial additional velocity of each particle
 * @param scale Number of particles in each direction
 * @param distance Relative distance between two particles
 * @param mass Mass of one particle
 * @param mean_brownian Average velocity of the Brownian Motion
 **/
constexpr void cuboid_generator(
	ParticleContainer& particles, const vec& origin, const vec_3d<int>& scale, double distance,
	const vec& initial_velocity, double mass, double brownian_mean
) {
	for (int i = 0; i < scale.x; ++i) {
		for (int j = 0; j < scale.y; ++j) {
			for (int k = 0; k < scale.z; ++k) {
				particles.emplace_back(
					vec{origin.x + (i * distance), origin.y + (j * distance), origin.z + (k * distance)},
					maxwellBoltzmannDistributedVelocity<3>(brownian_mean) + initial_velocity, mass
				);
			}
		}
	}
}
