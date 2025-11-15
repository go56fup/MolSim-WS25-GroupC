#pragma once

#include <cmath>
#include <functional>
#include <ranges>
#include <span>
#include <string_view>

#include <spdlog/spdlog.h>

#include "Concepts.h"
#include "ForceCalculators.h"
#include "IOProviders.h"
#include "outputWriter/outputWriters.h"
#include "outputWriter/VTKWriter.h"
#include "ParticleContainer.h"

constexpr void update_particle_forces(force_calculator auto calculator, Particle& p1, Particle& p2) noexcept {
	const auto f_ij = std::invoke(calculator, p1, p2);
	p1.f += f_ij;
	p2.f -= f_ij;
}

enum class bounds_type : std::uint8_t { none = 0, x_min, y_min, z_min, x_max, y_max, z_max };

constexpr bounds_type
bounds_check(const ParticleContainer::index& current, const ParticleContainer::index& boundaries) noexcept {
	if (current.x == 0) {
		return bounds_type::x_min;
	}
	if (current.y == 0) {
		return bounds_type::y_min;
	}
	if (current.z == 0) {
		return bounds_type::z_min;
	}
	if (current.x == boundaries.x) {
		return bounds_type::x_max;
	}
	if (current.y == boundaries.y) {
		return bounds_type::y_max;
	}
	if (current.z == boundaries.z) {
		return bounds_type::z_max;
	}
	return bounds_type::none;
};

constexpr void if_cell_at_border_reflect_via_ghost_particle(
	ParticleContainer::cell& current_cell, const ParticleContainer::index& current_cell_idx,
	const ParticleContainer::index& bounds, force_calculator auto calculator
) noexcept {
	static constexpr double ghost_particle_threshold = 1;

#define REFLECT_IF(check, x_mod, y_mod, z_mod)                                                                         \
	do {                                                                                                               \
		for (auto& p : current_cell) {                                                                                 \
			if (p.x.check) {                                                                                           \
				p.f += std::invoke(calculator, p, Particle{{p.x.x_mod, p.x.y_mod, p.x.z_mod}, {}, p.m});               \
			}                                                                                                          \
		}                                                                                                              \
	} while (0)

	switch (bounds_check(current_cell_idx, bounds)) {
	case bounds_type::none:
		return;
	case bounds_type::x_min:
		REFLECT_IF(x <= ghost_particle_threshold, x * -1, y, z);
		break;
	case bounds_type::y_min:
		REFLECT_IF(y <= ghost_particle_threshold, x, y * -1, z);
		break;
	case bounds_type::z_min:
		REFLECT_IF(z <= ghost_particle_threshold, x, y, z * -1);
		break;
	case bounds_type::x_max:
		REFLECT_IF(x >= bounds.x - ghost_particle_threshold, x + bounds.x, y, z);
		break;
	case bounds_type::y_max:
		REFLECT_IF(y >= bounds.y - ghost_particle_threshold, x, y + bounds.y, z);
		break;
	case bounds_type::z_max:
		REFLECT_IF(z >= bounds.z - ghost_particle_threshold, x, y, z + bounds.z);
		break;
	}
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
constexpr void calculate_forces(force_calculator auto calculator, ParticleContainer& particles) noexcept {
	for (const auto& [current_cell_idx, target_cell_idx] : particles.directional_interactions()) {
		// TODO(tuna): see if after the implementation of the border iterator whether we still need the indices
		SPDLOG_TRACE("Directionally interacting cells: {} -> {}", current_cell_idx, target_cell_idx);
		auto& current_cell = particles[current_cell_idx];
		auto& target_cell = particles[target_cell_idx];

		for (auto&& [p1, p2] : unique_pairs(current_cell)) {
			update_particle_forces(calculator, p1, p2);
		}
		for (auto& p1 : current_cell) {
			for (auto& p2 : target_cell) {
				update_particle_forces(calculator, p1, p2);
			}
		}
		// TODO(tuna): there must be a more efficient way to check this rather than checking for every cell
		if_cell_at_border_reflect_via_ghost_particle(current_cell, current_cell_idx, particles.grid_size(), calculator);
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
constexpr void calculateX(range_of<Particle> auto&& particles, double delta_t) noexcept {
	for (auto& p : particles) {
		const auto force_scalar = std::pow(delta_t, 2) / (2 * p.m);
		// TODO(tuna): check if expression templates are needed here
		p.x = p.x + delta_t * p.v + p.old_f * force_scalar;
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

// TODO(tuna): fix docs
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
	particle_io_provider auto&& io_provider, ParticleContainer& particles, std::string_view out_name, int iteration
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
constexpr void update_values(range_of<Particle> auto&& particles) noexcept {
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
constexpr void
run_sim_iteration(force_calculator auto force_calc, ParticleContainer& particles, double delta_t) noexcept {
	calculateX(particles.view(), delta_t);
	calculate_forces(force_calc, particles);
	calculateV(particles.view(), delta_t);
	update_values(particles.view());
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
			plotParticles(Traits.io, particles, output_prefix, iteration);
		}
		run_sim_iteration(Traits.force, particles, args.delta_t);
		current_time += args.delta_t;
		++iteration;
	}
}
