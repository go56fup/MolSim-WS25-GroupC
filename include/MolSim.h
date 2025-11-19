#pragma once

#include <cmath>
#include <functional>
#include <ranges>
#include <span>
#include <string_view>

#include "ForceCalculators.h"
#include "IOProviders.h"
#include "outputWriter/outputWriters.h"
#include "ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"

template <typename Range, typename Element>
concept range_of = std::ranges::range<Range> && (std::same_as<std::ranges::range_value_t<Range>, Element> ||
                                                 std::same_as<std::ranges::range_reference_t<Range>, Element>);


constexpr void calculate_force(force_calculator auto calculator, Particle& p1, Particle& p2) noexcept
	// TODO(tuna): see if this actually compiles down to a mask over simd like it should
	for (auto& p : particles) {
		p.f = {};
	}

	for (auto&& [p1, p2] : unique_pairs(particles)) {
		const auto f_ij = std::invoke(calculator, p1, p2);
		p1.f += f_ij;
		p2.f -= f_ij;
	}

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
constexpr void calculate_forces(force_calculator auto calculator, ParticleContainer& particles) noexcept
) {
	for (auto&& [current_cell, target_cell] : particles.directional_interactions()) {
		for (auto&& [p1, p2] : unique_pairs(current_cell)) {
			calculate_force(calculator, p1, p2);
		}
		for (auto&& [p1, p2] : std::views::cartesian_product(current_cell, dest_cell)) {
			calculate_force(calculator, p1, p2);
		}
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
constexpr void calculateX(force_calculator auto calculator, ParticleContainer& particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto force_scalar = std::pow(delta_t, 2) / (2 * p.m);
		// TODO(tuna): check if expression templates are needed here
		p.x = p.x + delta_t * p.v + p.old_f * force_scalar;
		while (p.on_or_beyond_boundary_cell(particles.domain())) {
			Particle ghost_particle(/* mirror */);
			calculate_force(calculator, p, ghost_particle);
			MUSTTAIL return calculateX(calculator, particles, delta_t);
		}
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
constexpr void calculateV(range_of<Particle> auto&& particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto velocity_scalar = delta_t / (2 * p.m);
		p.v = p.v + velocity_scalar * (p.old_f + p.f);
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
constexpr void plotParticles(
	particle_io_provider auto&& io_provider, range_of<Particle> auto&& particles, std::string_view out_name,
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

// TODO(anyone): do tparam docs
/** @brief Compile-time options for the simulation. */
template <particle_io_provider IOProvider, force_calculator ForceCalc>
struct sim_traits {
	ForceCalc force = lennard_jones_force;
	IOProvider io = OUTPUT_WRITER::plotParticles;
};

// TODO(anyone): update span references
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
template <sim_traits Traits = {.force = lennard_jones_force, .io = OUTPUT_WRITER::plotParticles}>
constexpr void
run_simulation(ParticleContainer& particles, const sim_args& args, std::string_view output_path) noexcept {
	double current_time = 0;
	int iteration = 0;
	constexpr auto plot_every_nth_iter = 10;
	const std::string output_prefix = std::string(output_path) + "/MD_vtk";

	while (current_time < args.end_time) {
		if (iteration % plot_every_nth_iter == 0) {
			plotParticles(Traits.io, particles.view(), output_prefix, iteration);
		}
		run_sim_iteration(Traits.force, particles, args.delta_t);
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
				particles.emplace(
					vec{origin.x + (i * distance), origin.y + (j * distance), origin.z + (k * distance)},
					maxwellBoltzmannDistributedVelocity<2>(brownian_mean) + initial_velocity, mass
				);
			}
		}
	}
}
