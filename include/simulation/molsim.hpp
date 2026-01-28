#pragma once

#include <cmath>
#include <functional>
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
#include "iterators/pairwise.hpp"
#include "iterators/periodic.hpp"
#include "output_writers/io.hpp"
#include "physics/forces.hpp"
#include "physics/particle.hpp"
#include "simulation/config/entities.hpp"
#include "simulation/config/parse.hpp"
#include "simulation/thermostat.hpp"
#include "utility/concepts.hpp"
#include "utility/tracing/macros.hpp"

/**
 * @brief Updates forces of two particles where one particle acts on the other.
 * @param calculator A callable object that computes forces for the particles.
 * Its type must satisfy `force_calculator`.
 * @param p1 The particle being acted upon.
 * @param p2 The particle exerting the force.
 */
constexpr void
apply_force_interaction(force_calculator auto calculator, particle& p1, particle& p2) noexcept {
	const auto f_ij = std::invoke(calculator, p1, p2);
	p1.f += f_ij;
	p2.f -= f_ij;
}

constexpr bool hasPeriodic(const sim_configuration& config) {
	using enum boundary_type;
	return config.boundary_behavior[x_min] == boundary_condition::periodic ||
	       config.boundary_behavior[y_min] == boundary_condition::periodic ||
	       config.boundary_behavior[z_min] == boundary_condition::periodic;
}

constexpr bool
isPeriodic(const sim_configuration& config, const std::vector<boundary_type>& bt_array) {
	bool result = true;
	for (auto boundary : bt_array) {
		result = result && config.boundary_behavior[boundary] == boundary_condition::periodic;
	}

	return result;
}

// Define the signature: void return, taking (cell, boundary, config)
using BoundaryFunc = std::function<
	void(particle_container::cell&, boundary_type, const particle_container::index& idx)>;

// TODO(tuna): move away from std::function
// Capture 'particles' and 'config' by reference (they live outside the function)
// Capture 'calculator' by VALUE (since it's a template/auto param that might be local)
BoundaryFunc funcDefiner(
	particle_container& particles, const sim_configuration& config,
	force_calculator auto calculator, boundary_type bType
) {

	auto behavior = config.boundary_behavior[bType];

	if (behavior == boundary_condition::reflecting) {
		return
			[&particles, calculator,
		     &config](particle_container::cell& cell, boundary_type b, const particle_container::index&) {
				reflect_via_ghost_particle(cell, particles.domain(), b, calculator);
			};
	} else if (behavior == boundary_condition::outflow) {
		return
			[&particles](particle_container::cell& cell, boundary_type b, const particle_container::index&) {
				delete_ouflowing_particles(cell, particles.domain(), b);
			};
	} else {
		assert(behavior == boundary_condition::periodic);
		return [&particles, calculator, &config](
				   particle_container::cell& cell, boundary_type b,
				   const particle_container::index& idx
			   ) { periodic(cell, b, calculator, idx, particles); };
	}
}

constexpr void loop(
	particle_container& particles, const sim_configuration& config, force_calculator auto calculator
) {
	const auto& grid = particles.grid_size();
	using enum boundary_type;
	// TODO: Map to function

	// x faces
	BoundaryFunc func = funcDefiner(particles, config, calculator, x_min);
	for (unsigned i = 0; i < grid.y; ++i) {
		for (unsigned j = 0; j < grid.z; ++j) {
			const particle_container::index min_coords{0, i, j};
			const particle_container::index max_coords{grid.x - 1, i, j};
			TRACE_BORDER_CELL_ITER("min: {}, max: {}", min_coords, max_coords);
			std::invoke(func, particles[min_coords], x_min, max_coords);
			std::invoke(func, particles[max_coords], x_max, max_coords);
		}
	}

	// y faces
	func = funcDefiner(particles, config, calculator, y_min);
	for (unsigned i = 0; i < grid.x; ++i) {
		for (unsigned j = 0; j < grid.z; ++j) {
			const particle_container::index min_coords{i, 0, j};
			const particle_container::index max_coords{i, grid.y - 1, j};
			TRACE_BORDER_CELL_ITER("min: {}, max: {}", min_coords, max_coords);
			std::invoke(func, particles[min_coords], y_min, min_coords);
			std::invoke(func, particles[max_coords], y_max, max_coords);
		}
	}

	// z faces
	func = funcDefiner(particles, config, calculator, z_min);
	for (unsigned i = 0; i < grid.x; ++i) {
		for (unsigned j = 0; j < grid.y; ++j) {
			const particle_container::index min_coords{i, j, 0};
			const particle_container::index max_coords{i, j, grid.z - 1};
			TRACE_BORDER_CELL_ITER("min: {}, max: {}", min_coords, max_coords);
			std::invoke(func, particles[min_coords], z_min, max_coords);
			std::invoke(func, particles[max_coords], z_max, max_coords);
		}
	}

	if (!hasPeriodic(config)) {
		return;
	}
	func = [&](particle_container::cell& cell, boundary_type b, const particle_container::index& idx
	       ) { periodic(cell, b, calculator, idx, particles); };

	// x-y edges
	bool cond = isPeriodic(config, {x_min, y_min});
	if (cond) {
		for (unsigned i = 0; i < grid.z; ++i) {
			std::invoke(
				func, particles[{0, 0, i}], x_min | y_min, particle_container::index{0, 0, i}
			);
			std::invoke(
				func, particles[{0, grid.y, i}], x_min | y_max,
				particle_container::index{0, grid.y, i}
			);
			std::invoke(
				func, particles[{grid.x, 0, i}], x_max | y_min,
				particle_container::index{grid.x, 0, i}
			);
			std::invoke(
				func, particles[{grid.x, grid.y, i}], x_max | y_max,
				particle_container::index{grid.x, grid.y, i}
			);
		}
	}

	// x-z edges
	cond = isPeriodic(config, {x_min, z_min});
	if (cond) {
		for (unsigned i = 0; i < grid.y; ++i) {
			std::invoke(
				func, particles[{0, i, 0}], x_min | z_min, particle_container::index{0, i, 0}
			);
			std::invoke(
				func, particles[{0, i, grid.z}], x_min | z_max,
				particle_container::index{0, i, grid.z}
			);
			std::invoke(
				func, particles[{grid.x, i, 0}], x_max | z_min,
				particle_container::index{grid.x, i, 0}
			);
			std::invoke(
				func, particles[{grid.x, i, grid.z}], x_max | z_max,
				particle_container::index{grid.x, i, grid.z}
			);
		}
	}

	// y-z edges
	cond = isPeriodic(config, {y_min, z_min});
	if (cond) {
		for (unsigned i = 0; i < grid.y; ++i) {
			std::invoke(
				func, particles[{i, 0, 0}], y_min | z_min, particle_container::index{i, 0, 0}
			);
			std::invoke(
				func, particles[{i, 0, grid.z}], y_min | z_max,
				particle_container::index{i, 0, grid.z}
			);
			std::invoke(
				func, particles[{i, grid.y, 0}], y_max | z_min,
				particle_container::index{i, grid.y, 0}
			);
			std::invoke(
				func, particles[{i, grid.y, grid.z}], y_max | z_max,
				particle_container::index{i, grid.y, grid.z}
			);
		}
	}

	// x-y-z corners
	cond = isPeriodic(config, {x_min, y_min, z_min});
	if (cond) {
		std::invoke(
			func, particles[{0, 0, 0}], x_min | y_min | z_min, particle_container::index{0, 0, 0}
		);
		std::invoke(
			func, particles[{0, 0, grid.z}], x_min | y_min | z_max,
			particle_container::index{0, 0, grid.z}
		);
		std::invoke(
			func, particles[{0, grid.y, 0}], x_min | y_max | z_min,
			particle_container::index{0, grid.y, 0}
		);
		std::invoke(
			func, particles[{0, grid.y, grid.z}], x_min | y_max | z_max,
			particle_container::index{0, grid.y, grid.z}
		);
		std::invoke(
			func, particles[{grid.x, 0, 0}], x_max | y_min | z_min,
			particle_container::index{grid.x, 0, 0}
		);
		std::invoke(
			func, particles[{grid.x, 0, grid.z}], x_max | y_min | z_max,
			particle_container::index{grid.x, 0, grid.z}
		);
		std::invoke(
			func, particles[{grid.x, grid.y, 0}], x_max | y_max | z_min,
			particle_container::index{grid.x, grid.y, 0}
		);
		std::invoke(
			func, particles[{grid.x, grid.y, grid.z}], x_max | y_max | z_max,
			particle_container::index{grid.x, grid.y, grid.z}
		);
	}
}

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
constexpr void calculate_forces(
	force_calculator auto calculator, particle_container& container, const sim_configuration& config
) noexcept {
	for (auto& cell : container.cells()) {
		for (auto&& [p1, p2] : unique_pairs(cell)) {
			apply_force_interaction(calculator, p1, p2);
		}
	}

	for (const auto& [current_cell_idx, target_cell_idx] : container.directional_interactions()) {
		// TODO(tuna): see if after the implementation of the border iterator whether we still
		// need the indices
		auto& current_cell = container[current_cell_idx];
		auto& target_cell = container[target_cell_idx];

		for (auto& p1 : current_cell) {
			for (auto& p2 : target_cell) {
				apply_force_interaction(calculator, p1, p2);
			}
		}
	}
	loop(container, config, calculator);
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
calculate_x(particle_container& particles, const sim_configuration& config) noexcept(false) {
	// TODO(tuna): rewrite in terms of a common Bounds.h that has this displacement logic,
	// switch over a boundary_type and set new i, j, k -> at the end of the loop do the
	// index, push_back, erase.

	const double cell_width = particles.cutoff_radius();
	auto move_to_cell = [&](const particle_container::index& new_cell_idx,
	                        particle_container::cell& current_cell,
	                        particle_container::cell::size_type i) noexcept {
		const auto& particle = current_cell[i];
		TRACE_GRID("Moving {} to {}", particle, new_cell_idx);
		particle_container::cell& new_cell = particles[new_cell_idx];
		new_cell.push_back(particle);
		current_cell.erase(std::next(
			current_cell.begin(), static_cast<particle_container::cell::difference_type>(i)
		));
	};
	const auto& grid = particles.grid_size();

	for (auto&& [cell_idx, cell] : particles.enumerate_cells()) {
		const vec cell_origin{
			cell_idx.x * cell_width, cell_idx.y * cell_width, cell_idx.z * cell_width
		};

		for (std::size_t particle_idx = 0; particle_idx < cell.size(); ++particle_idx) {
			particle_container::index_diff displacement{};
			auto& p = cell[particle_idx];

			const auto force_scalar = (config.delta_t * config.delta_t) / (2 * p.m);
			// TODO(tuna): check if expression templates are needed here
			p.x = p.x + config.delta_t * p.v + p.old_f * force_scalar;
			auto& pos = p.x;
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
					"{}, illegality={}, current: {}, displacement: {}, origin: {}", p,
					becomes_illegal, cell_idx, displacement, cell_origin
				);
				if (becomes_illegal) {
					if (config.boundary_behavior[b] == boundary_condition::periodic) {
						// Zero out the current position first, teleport to 0 if overflowing, to end
						// of grid if underflowing.
						TRACE_PERIODIC(
							"Performing periodic teleport on {}, currently stored in {}", p,
							cell_idx
						);
						displacement[axis_] = -cell_idx[axis_];
						displacement[axis_] += is_min ? grid[axis_] - 1 : 0;
						const double update =
							increment * static_cast<double>(grid[axis_]) * cell_width;
						pos[axis_] -= update;
						TRACE_PERIODIC(
							"Periodic displacement: {}, pos update: -{}, new position {}",
							displacement, update, pos
						);
					} else {
						// TODO(tuna): see if putting a std::unreachable here and only throwing on
						// debug makes any difference
						throw std::out_of_range(fmt::format(
							"Trying to move {} to non-existant cell with index {} + {} on "
							"grid {}",
							p, cell_idx, displacement, grid
						));
					}
				}
				new_cell_idx[axis_] = apply_difference(cell_idx[axis_], displacement[axis_]);
			};

			using enum boundary_type;

			if (pos.x < cell_origin.x) {
				update_displacement.template operator()<x_min>();
			} else if (pos.x >= cell_origin.x + cell_width) {
				update_displacement.template operator()<x_max>();
			} else if (pos.y < cell_origin.y) {
				update_displacement.template operator()<y_min>();
			} else if (pos.y >= cell_origin.y + cell_width) {
				update_displacement.template operator()<y_max>();
			} else if (pos.z < cell_origin.z) {
				update_displacement.template operator()<z_min>();
			} else if (pos.z >= cell_origin.z + cell_width) {
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
constexpr void calculate_v(range_of<particle> auto&& particles, double delta_t) noexcept {
	for (auto& p : particles) {
		TRACE_SIM("Old V: {}", p);
		const auto velocity_scalar = delta_t / (2 * p.m);
		p.v = p.v + velocity_scalar * (p.old_f + p.f);
		TRACE_SIM("New V: {}", p);
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
constexpr void update_values(range_of<particle> auto&& particles) noexcept {
	for (auto& p : particles) {
		p.old_f = p.f;
		p.f = {};
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
	force_calculator auto force_calc, particle_container& container,
	const sim_configuration& config, sim_iteration_t iteration
) noexcept {
	if constexpr (!IsFirstIteration) {
		update_values(container.particles());
	}

	static const bool has_thermostat = config.thermostat.has_value();
	static const bool has_gravity = config.gravitational_constant != 0;
	calculate_x(container, config);
	calculate_forces(force_calc, container, config);

	if (has_thermostat && iteration % config.thermostat->application_frequency == 0) {
		run_thermostat(container, *config.thermostat, config.dimensions);
	}

	if (has_gravity) {
		apply_gravity(container.particles(), config.gravitational_constant);
	}

	calculate_v(container.particles(), config.delta_t);
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
constexpr void run_simulation(
	particle_container& container, const sim_configuration& config,
	force_calculator auto force_calc, std::string_view output_path
) noexcept {
	double current_time = 0;
	sim_iteration_t iteration = 0;
	const std::string output_prefix = std::string(output_path) + "/" + config.base_name.c_str();

	struct sim_traits {
		bool is_first_iteration = false;
	};

	auto run_sim_pass = [&]<sim_traits Traits = {}> {
		TRACE_SIM("beginning iteration, current_time={}", current_time);
		if (iteration % config.write_frequency == 0) {
			WRITE_VTK_OUTPUT(plot_particles, container, output_prefix, iteration);
		}

		run_sim_iteration<Traits.is_first_iteration>(force_calc, container, config, iteration);
		current_time += config.delta_t;
		++iteration;
	};

	auto start = std::chrono::steady_clock::now();
	run_sim_pass.template operator()<{.is_first_iteration = true}>();
	while (current_time < config.end_time) {
		run_sim_pass.template operator()();
	}
	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;

	SPDLOG_INFO(
		"Runtime: {} ms, iterations: {}, tick length: {} ms, MUPS/s: {}", elapsed.count(),
		iteration, elapsed.count() / iteration, (iteration / elapsed.count()) * 1000
	);

	if (config.create_checkpoint) {
		// TODO(tuna): specify both in terms of plot_particles, where iteration is used in the
		// filename Or just remove plot_particles
		WRITE_CHECKPOINT(container, fmt::format("{}_checkpoint.json", output_prefix));
	}
}
