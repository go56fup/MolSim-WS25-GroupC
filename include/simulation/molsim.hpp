#pragma once

#include <cmath>
#include <functional>
#include <omp.h>
#include <ranges>
#include <span>
#include <string_view>

#include <fmt/compile.h>
#include <spdlog/spdlog.h>

#include "config/parse.hpp"

#include "grid/bounds/conditions.hpp"
#include "grid/bounds/operations.hpp"
#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
#include "grid/particle_container/particle_container.hpp"
#include "grid/particle_container/system.hpp"
#include "iterators/pairwise.hpp"
#include "iterators/periodic.hpp"
#include "output_writers/io.hpp"
#include "physics/forces.hpp"
#include "physics/particle.hpp"
#include "simulation/config/entities.hpp"
#include "simulation/config/parse.hpp"
#include "simulation/thermostat.hpp"
#include "utility/concepts.hpp"
#include "utility/macros.hpp"
#include "utility/tracing/macros.hpp"

/**
 * @brief Updates forces of two particles where one particle acts on the other.
 * @param calculator A callable object that computes forces for the particles.
 * Its type must satisfy `force_calculator`.
 * @param p1 The particle being acted upon.
 * @param p2 The particle exerting the force.
 */
/**
 * @brief Calculates forces on a collection of particles using the provided force calculator.
 *
 * This function delegates the force calculation to the provided @p calculator callable,
 * allowing different calculation methods.
 *
 * @param calculator A callable object that computes forces for a span of particles.
 * Its type must satisfy force_calculator.
 * @param particles The span over particles on which forces will be calculated.
 *
 */

constexpr void
calculate_forces_batched(particle_container& container, const sim_configuration& config) noexcept {
	auto use_batch_piecewise = [&](particle_batch batch_p1, particle_batch batch_p2,
	                               std::size_t up_to = batch_size) {
		for (std::size_t i = 0; i < up_to; ++i) {
			// TODO(gabriel): Also we must fully remove Newtons Axiom from our code.
			lennard_jones_force_soa(container, batch_p1[i], batch_p2[i]);
		}
	};

	auto use_batch = [&](particle_batch batch_p1, particle_batch batch_p2) {
		lennard_jones_force_soa_batchwise(container, batch_p1, batch_p2);
	};

#ifndef SINGLETHREADED
#pragma omp parallel for schedule(dynamic)
#endif
	for (auto& cell : container.cells()) {
		std::size_t count = 0;
		particle_batch batch_p1;
		particle_batch batch_p2;

		for (auto [p1_idx, p2_idx] : unique_pairs(cell)) {
			TRACE_FORCES(
				"Putting {} {} into batch, count={} on thread {}", p1_idx, p2_idx, count,
				omp_get_thread_num()
			);
			batch_p1[count] = p1_idx;
			batch_p2[count] = p2_idx;
			++count;

			if (count == batch_size) {
				use_batch(batch_p1, batch_p2);
				count = 0;
			}
		}
		use_batch_piecewise(batch_p1, batch_p2, count);
	}

	for (const auto& [current_cell_idx, target_cell_idx] : container.directional_interactions()) {
		// TODO(tuna): see if after the implementation of the border iterator whether we still
		// need the indices
		const auto& current_cell = container[current_cell_idx];
		const auto& target_cell = container[target_cell_idx];

		std::size_t count = 0;
		particle_batch batch_p1;
		particle_batch batch_p2;

		for (auto p1_idx : current_cell) {
			for (auto p2_idx : target_cell) {
				batch_p1[count] = p1_idx;
				batch_p2[count] = p2_idx;
				++count;

				if (count == batch_size) {
					use_batch(batch_p1, batch_p2);
					count = 0;
				}
			}
		}
		use_batch_piecewise(batch_p1, batch_p2, count);
	}

	handle_boundaries(container, config);
}

/**
 * @brief Updates the position of each particle.
 *
 * Implements the first step of the velocity Störmer-Verlet algorithm:
 * \f[
 * \vec{x}(t + \Delta t) = \vec{x}(t) + \Delta t \cdot \vec{v}(t) + \frac{\Delta t^2}{2m}
 * \vec{F}(t)
 * \f]
 *
 * @param particles Mutable span over particles to update.
 * @param delta_t The time step for integration.
 */
constexpr void
calculate_x(particle_container& container, const sim_configuration& config) noexcept(false) {
	const double cell_width = container.cutoff_radius();
	auto move_to_cell = [&](const particle_container::index& new_cell_idx,
	                        particle_container::cell& current_cell,
	                        particle_container::cell::size_type i) noexcept {
		const auto& particle = current_cell[i];
		TRACE_GRID("Moving {} to {}", particle, new_cell_idx);
		particle_container::cell& new_cell = container[new_cell_idx];
		new_cell.push_back(particle);
		current_cell.erase(std::next(
			current_cell.begin(), static_cast<particle_container::cell::difference_type>(i)
		));
	};
	const auto& grid = container.grid_size();
	auto& system = container.system();

#ifndef SINGLETHREADED
	// #pragma omp parallel for schedule(dynamic)
#endif
	for (auto&& [cell_idx, cell] : container.enumerate_cells()) {
		const vec cell_origin{
			cell_idx.x * cell_width, cell_idx.y * cell_width, cell_idx.z * cell_width
		};

		for (std::size_t particle_idx = 0; particle_idx < cell.size(); ++particle_idx) {
			particle_container::signed_index displacement{};
			auto p_idx = cell[particle_idx];

			const auto force_scalar = (config.delta_t * config.delta_t) / (2 * system.mass[p_idx]);
			system.x[p_idx] +=
				config.delta_t * system.vx[p_idx] + system.old_fx[p_idx] * force_scalar;
			system.y[p_idx] +=
				config.delta_t * system.vy[p_idx] + system.old_fy[p_idx] * force_scalar;
			system.z[p_idx] +=
				config.delta_t * system.vz[p_idx] + system.old_fz[p_idx] * force_scalar;
			auto new_cell_idx = cell_idx;

			// TODO(tuna): see if reformulating interactions' do_displacement via something common
			// to this is possible
			auto update_displacement = [&]<boundary_type b> {
				static constexpr auto axis_ = boundary_type_to_axis(b);
				static constexpr auto is_min = boundary_type_to_extremum(b) == extremum::min;
				static constexpr auto increment = is_min ? -1 : 1;
				displacement[axis_] = increment;
				const bool becomes_illegal = cell_idx[axis_] == (is_min ? 0 : grid[axis_] - 1);
				TRACE_GRID(
					"{}, illegality={}, current: {}, displacement: {}, origin: {}", p_idx,
					becomes_illegal, cell_idx, displacement, cell_origin
				);
				if (becomes_illegal) {
					if (config.boundary_behavior[b] == boundary_condition::periodic) {
						// Zero out the current position first, teleport to 0 if overflowing, to end
						// of grid if underflowing.
						TRACE_PERIODIC(
							"Performing periodic teleport on {}, currently stored in {}", p_idx,
							cell_idx
						);
						displacement[axis_] =
							-static_cast<particle_container::difference_type>(cell_idx[axis_]);
						displacement[axis_] += is_min ? grid[axis_] - 1 : 0;
						const double update =
							increment * static_cast<double>(grid[axis_]) * cell_width;
						system.position_component(axis_)[p_idx] -= update;
						TRACE_PERIODIC(
							"Periodic displacement: {}, pos update: -{}, new position {}",
							displacement, update, system.serialize_position(p_idx)
						);
					} else {
						// TODO(tuna): see if putting a std::unreachable here and only throwing on
						// debug makes any difference
						throw std::out_of_range(fmt::format(
							"Trying to move {} to non-existant cell with index {} + {} on "
							"grid {}",
							p_idx, cell_idx, displacement, grid
						));
					}
				}
				new_cell_idx[axis_] = apply_difference(cell_idx[axis_], displacement[axis_]);
			};

			using enum boundary_type;

			if (system.x[p_idx] < cell_origin.x) {
				update_displacement.template operator()<x_min>();
			} else if (system.x[p_idx] >= cell_origin.x + cell_width) {
				update_displacement.template operator()<x_max>();
			} else if (system.y[p_idx] < cell_origin.y) {
				update_displacement.template operator()<y_min>();
			} else if (system.y[p_idx] >= cell_origin.y + cell_width) {
				update_displacement.template operator()<y_max>();
			} else if (system.z[p_idx] < cell_origin.z) {
				update_displacement.template operator()<z_min>();
			} else if (system.z[p_idx] >= cell_origin.z + cell_width) {
				update_displacement.template operator()<z_max>();
			} else {
				// Skip moving if all checks succeed
				continue;
			}
			move_to_cell(new_cell_idx, cell, particle_idx);
		}
	}
}

/**
 * @brief Updates the velocity of each particle.
 *
 * Implements the second step of the velocity Störmer-Verlet algorithm:
 * \f[
 * \vec{v}(t + \Delta t) = \vec{v}(t) + \frac{\Delta t}{2m} \left[\vec{F}(t) + \vec{F}(t +
 * \Delta t)\right]
 * \f]
 *
 * @param particles Mutable span over particles to update.
 * @param delta_t The time step for integration.
 */
constexpr void calculate_v(particle_container& container, double delta_t) noexcept {
	auto& system = container.system();

#ifndef SINGLETHREADED
#pragma omp parallel for simd schedule(static)
#endif
	for (particle_system::particle_id i = 0; i < container.system().size(); ++i) {
		// TODO(gabriel): Maybe add a invM array and * 0.5 instead of / 2 to remove all these costly
		// divisions
		const auto velocity_scalar = delta_t / (2 * system.mass[i]);
		system.vx[i] += velocity_scalar * (system.old_fx[i] + system.fx[i]);
		system.vy[i] += velocity_scalar * (system.old_fy[i] + system.fy[i]);
		system.vz[i] += velocity_scalar * (system.old_fz[i] + system.fz[i]);
	}
}

// TODO(tuna): fix docs
/**
 * @brief Exports particle data using the given I/O provider.
 *
 * This function delegates to the given @p io_provider callable, allowing different
 * data exporting strategies.
 *
 * @param io_provider Callable object responsible for exporting particle data. Its type must
 * satisfy particle_io_provider.
 * @param particles Constant span of particles to export data from.
 * @param out_name  Base name for the output file.
 * @param iteration Current simulation iteration (used for output naming).
 */
constexpr void plot_particles(
	particle_io_provider auto&& io_provider, particle_container& particles,
	std::string_view out_name, unsigned iteration
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
constexpr void update_values(particle_system& system) noexcept {
#ifndef SINGLETHREADED
#pragma omp parallel for simd schedule(static)
#endif
	for (particle_system::particle_id p = 0; p < system.size(); ++p) {
		system.old_fx[p] = system.fx[p];
		system.old_fy[p] = system.fy[p];
		system.old_fz[p] = system.fz[p];
		system.fx[p] = 0;
		system.fy[p] = 0;
		system.fz[p] = 0;
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
// TODO(tuna): now that iteration is being passed in, isFirstIteration should probably be removed
template <bool IsFirstIteration = false>
constexpr void run_sim_iteration(
	particle_container& container, const sim_configuration& config, sim_iteration_t iteration
) noexcept {
	if constexpr (!IsFirstIteration) {
		update_values(container.system());
	}

	// TODO(tuna): mark if_not_testing
	static const bool has_thermostat = config.thermostat.has_value();
	static const bool has_gravity = config.gravitational_constant != 0;
	calculate_x(container, config);
	calculate_forces_batched(container, config);

	if (has_thermostat && iteration % config.thermostat->application_frequency == 0) {
		run_thermostat(container, *config.thermostat, config.dimensions);
	}

	if (has_gravity) {
		apply_gravity(container, config.gravitational_constant);
	}

	calculate_v(container, config.delta_t);
#if LOG_SIM
	const auto& system = container.system();
	for (particle_system::particle_id i = 0; i < system.size(); ++i) {
		TRACE_SIM(
			"End of iteration, particle {}: pos={}, v={}, f={}, old_f={}", i,
			system.serialize_position(i), system.serialize_velocity(i), system.serialize_force(i),
			system.serialize_old_force(i)
		);
	}
#endif
}

// TODO(anyone): update span references
/**
 * @brief Start the simulation.
 *
 * Calculates new positions, velocities and forces at each time tick
 * and plots the updated particles.
 *
 * @param container Particles to simulate.
 * @param config Simulation parameters.
 * @param force_calc Force calculation method to use for simulation.
 * @param output_path Path to put simulation result files into.
 **/
constexpr sim_iteration_t run_simulation(
	particle_container& container, const sim_configuration& config, std::string_view output_prefix
) noexcept {
	double current_time = 0;
	sim_iteration_t iteration = 0;

	struct sim_traits {
		bool is_first_iteration = false;
	};

	auto run_sim_pass = [&]<sim_traits Traits = {}> {
		TRACE_SIM("beginning iteration, current_time={}", current_time);
		if (iteration % config.write_frequency == 0) {
			WRITE_VTK_OUTPUT(plot_particles, container, output_prefix, iteration);
		}

		run_sim_iteration<Traits.is_first_iteration>(container, config, iteration);
		current_time += config.delta_t;
		++iteration;
	};

	run_sim_pass.template operator()<{.is_first_iteration = true}>();
	while (current_time < config.end_time) {
		run_sim_pass.template operator()();
	}
	return iteration;
}
