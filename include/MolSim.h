#pragma once

#include <cmath>
#include <functional>
#include <ranges>
#include <span>
#include <string_view>

#include <spdlog/spdlog.h>

#include "Concepts.h"
#include "FileReader.h"
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

// The wrapping in do-while increases this metric dramatically. Otherwise, the rest of the logic is impossible
// to express in more granular functions without having to recompute information. The template for facility
// from C++26 would have come in useful here.
// NOLINTBEGIN(*cognitive-complexity)
// Useful to create a block inside the macro.
// NOLINTBEGIN(*avoid-do-while)
constexpr void reflect_via_ghost_particle(
	ParticleContainer::cell& current_cell, const ParticleContainer::index& bounds, boundary_type border_type,
	force_calculator auto calculator, double ghost_particle_threshold
) noexcept {
#define REFLECT_IF(check, x_mod, y_mod, z_mod)                                                                         \
	do {                                                                                                               \
		for (auto& p : current_cell) {                                                                                 \
			if (p.x.check) {                                                                                           \
				SPDLOG_TRACE("{} reflecting, bc: {} -- ghost: {} {} {}", p, #check, #x_mod, #y_mod, #z_mod);           \
				p.f += std::invoke(calculator, p, Particle{{p.x.x_mod, p.x.y_mod, p.x.z_mod}, {}, p.m});               \
			}                                                                                                          \
		}                                                                                                              \
	} while (0)

	using enum boundary_type;
	switch (border_type) {

	case x_min:
		REFLECT_IF(x <= ghost_particle_threshold, x * -1, y, z);
		break;
	case y_min:
		REFLECT_IF(y <= ghost_particle_threshold, x, y * -1, z);
		break;
	case z_min:
		REFLECT_IF(z <= ghost_particle_threshold, x, y, z * -1);
		break;
	case x_max:
		REFLECT_IF(x >= bounds.x - ghost_particle_threshold, x + bounds.x, y, z);
		break;
	case y_max:
		REFLECT_IF(y >= bounds.y - ghost_particle_threshold, x, y + bounds.y, z);
		break;
	case z_max:
		REFLECT_IF(z >= bounds.z - ghost_particle_threshold, x, y, z + bounds.z);
		break;
	case none:
		// TODO(tuna): make the iterator smarter so it doesn't have to check for anything, and just knows
		assert(false && "This should never be called with a non-border, the iterator is wrong.");
	}
}

// NOLINTEND(*avoid-do-while)
// NOLINTEND(*cognitive-complexity)

template <axis Axis>
constexpr bool out_of_bounds(const Particle& p, const ParticleContainer::index& bounds) noexcept {
	static constexpr auto proj_double = axis_ptr<double>(Axis);
	static constexpr auto proj_size_t = axis_ptr<ParticleContainer::size_type>(Axis);
	const auto& pos = p.x.*proj_double;
	return pos <= 0 || pos >= bounds.*proj_size_t;
}

constexpr void delete_ouflowing_particles(ParticleContainer::cell& cell, const ParticleContainer::index& bounds) {
	auto end = cell.end();
	for (auto iter = cell.begin(); iter != end; ++iter) {
		if (out_of_bounds<axis::x>(*iter, bounds) || out_of_bounds<axis::y>(*iter, bounds) ||
		    out_of_bounds<axis::z>(*iter, bounds)) {
			cell.erase(iter);
		}
	}
}

constexpr void handle_boundary_condition(
	ParticleContainer::cell& cell, boundary_type border_type, const ParticleContainer::index& bounds,
	force_calculator auto force_calc, const sim_configuration& config
) noexcept {
	switch (config.boundary_behavior[border_type]) {
	case boundary_condition::reflecting:
		reflect_via_ghost_particle(cell, bounds, border_type, force_calc, config.meshwidth);
		break;
	case boundary_condition::outflow:
		delete_ouflowing_particles(cell, bounds);
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
constexpr void calculate_forces(
	force_calculator auto calculator, ParticleContainer& particles, const sim_configuration& config
) noexcept {
	for (auto& cell : particles) {
		for (auto&& [p1, p2] : unique_pairs(cell)) {
			SPDLOG_TRACE("Got pair of particles: {}, {}", p1, p2);
			update_particle_forces(calculator, p1, p2);
		}
	}
	for (const auto& [current_cell_idx, target_cell_idx] : particles.directional_interactions()) {
		// TODO(tuna): see if after the implementation of the border iterator whether we still need the indices
		SPDLOG_TRACE("Directionally interacting cells: {} -> {}", current_cell_idx, target_cell_idx);
		auto& current_cell = particles[current_cell_idx];
		auto& target_cell = particles[target_cell_idx];

		for (auto& p1 : current_cell) {
			for (auto& p2 : target_cell) {
				SPDLOG_TRACE("Found cross-cell particles to interact: {} -> {}", p1, p2);
				update_particle_forces(calculator, p1, p2);
			}
		}
	}
	for (auto&& [cell_idx, boundary_type] : particles.border_cells()) {
		auto& cell = particles[cell_idx];
		SPDLOG_TRACE("Doing border cells: {}", std::to_underlying(boundary_type));
		const auto& bounds = particles.grid_size();
		handle_boundary_condition(cell, boundary_type, bounds, calculator, config);
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
constexpr void calculateX(ParticleContainer& particles, double delta_t) noexcept {
	// TODO(tuna): rewrite in terms of a common Bounds.h that has this displacement logic,
	// switch over a boundary_type and set new i, j, k -> at the end of the loop do the
	// index, push_back, erase.

	const double cell_width = particles.cutoff_radius();
	auto move_to_cell = [&](const ParticleContainer::index& new_cell_idx, ParticleContainer::cell& current_cell,
	                        ParticleContainer::cell::iterator particle_iter) {
		SPDLOG_TRACE("Moving {} to {}", *particle_iter, new_cell_idx);
		ParticleContainer::cell& new_cell = particles[new_cell_idx];
		new_cell.push_back(*particle_iter);
		current_cell.erase(particle_iter);
	};

	for (auto&& [idx, cell] : particles.enumerate_cells()) {
		const vec cell_origin{idx.x * cell_width, idx.y * cell_width, idx.z * cell_width};
		// out of cell check
		auto end = cell.end();
		for (auto iter = cell.begin(); iter != end; ++iter) {
			auto& particle = *iter;

			const auto force_scalar = std::pow(delta_t, 2) / (2 * particle.m);
			// TODO(tuna): check if expression templates are needed here
			particle.x = particle.x + delta_t * particle.v + particle.old_f * force_scalar;

			if (particle.x.x < cell_origin.x) {  // for x direction
				SPDLOG_TRACE("performing below x move");
				move_to_cell({idx.x - 1, idx.y, idx.z}, cell, iter);
			} else if (particle.x.x >= cell_origin.x + cell_width) {
				SPDLOG_TRACE("performing above x move");
				move_to_cell({idx.x + 1, idx.y, idx.z}, cell, iter);
			} else if (particle.x.y < cell_origin.y) {  // for y direction
				SPDLOG_TRACE("performing below y move");
				move_to_cell({idx.x, idx.y - 1, idx.z}, cell, iter);
			} else if (particle.x.y >= cell_origin.y + cell_width) {
				SPDLOG_TRACE("performing above y move");
				move_to_cell({idx.x, idx.y + 1, idx.z}, cell, iter);
			} else if (particle.x.z < cell_origin.z) {  // for z direction
				SPDLOG_TRACE("performing below z move");
				move_to_cell({idx.x, idx.y, idx.z - 1}, cell, iter);
			} else if (particle.x.z >= cell_origin.z + cell_width) {
				SPDLOG_TRACE("performing above z move");
				move_to_cell({idx.x, idx.y, idx.z + 1}, cell, iter);
			}
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
	particle_io_provider auto&& io_provider, ParticleContainer& particles, std::string_view out_name, unsigned iteration
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

/**
 * @brief Run a tick of the simulation.
 *
 * @param force_calc A callable object that computes forces for a span of particles.
 * Its type must satisfy force_calculator.
 * @param particles Mutable span over particles to run simulation tick for
 * @param delta_t Tick length
 */
constexpr void run_sim_iteration(
	force_calculator auto force_calc, ParticleContainer& particles, const sim_configuration& config
) noexcept {
	calculateX(particles, config.delta_t);
	calculate_forces(force_calc, particles, config);
	calculateV(particles.view(), config.delta_t);
	update_values(particles.view());
}

// TODO(anyone): do tparam docs
/** @brief Compile-time options for the simulation. */
template <particle_io_provider IOProvider>
struct sim_traits {
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
template <sim_traits Traits = {.io = OUTPUT_WRITER::plotParticles}>
constexpr void run_simulation(
	ParticleContainer& particles, const sim_configuration& config, force_calculator auto force_calc,
	std::string_view output_path
) noexcept {
	double current_time = 0;
	decltype(sim_configuration::write_frequency) iteration = 0;
	const std::string output_prefix = fmt::format("{}/{}", output_path, config.base_name);

	while (current_time < config.end_time) {
		if (iteration % config.write_frequency == 0) {
			plotParticles(Traits.io, particles, output_prefix, iteration);
		}
		run_sim_iteration(force_calc, particles, config);
		current_time += config.delta_t;
		++iteration;
	}
}
