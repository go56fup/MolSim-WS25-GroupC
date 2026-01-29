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
#include "utility/tracing/macros.hpp"

#include <__xlocale.h>

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
template <std::size_t BatchSize>
constexpr void calculate_forces_batched(
	force_calculator auto calculator, particle_container& container, const sim_configuration& config
) noexcept {
	using batch = std::array<particle_system::particle_id, BatchSize>;

	auto use_batch = [&](batch& batch_p1, batch& batch_p2, std::size_t up_to = BatchSize) {
		TRACE_FORCES(
			"Using &batch_p1={}, &batch_p2={} on thread {}", static_cast<void*>(&batch_p1),
			static_cast<void*>(&batch_p2), omp_get_thread_num()
		);
		for (std::size_t i = 0; i < up_to; ++i) {
			// TODO(gabriel): Also we must fully remove Newtons Axiom from our code.
			// TODO(gabriel): does it rly need to be an invoke do we ever use this not as
			// lennard jones forces?
			std::invoke(
				calculator, container, batch_p1[i], batch_p2[i]
			);
		}
	};

#pragma omp parallel for schedule(dynamic)
	for (auto& cell : container.cells()) {
		TRACE_FORCES("Doing cell {} on thread {}", static_cast<void*>(&cell), omp_get_thread_num());

		std::size_t count = 0;
		batch batch_p1;
		batch batch_p2;

		for (auto [p1_idx, p2_idx] : unique_pairs(cell)) {
			TRACE_FORCES(
				"Putting {} {} into batch, count={} on thread {}", p1_idx, p2_idx, count,
				omp_get_thread_num()
			);
			batch_p1[count] = p1_idx;
			batch_p2[count] = p2_idx;
			++count;

			if (count == BatchSize) {
				use_batch(batch_p1, batch_p2);
				count = 0;
			}
		}
		use_batch(batch_p1, batch_p2, count);
	}

	for (const auto& [current_cell_idx, target_cell_idx] : container.directional_interactions()) {
		// TODO(tuna): see if after the implementation of the border iterator whether we still
		// need the indices
		const auto& current_cell = container[current_cell_idx];
		const auto& target_cell = container[target_cell_idx];

		std::size_t count = 0;
		batch batch_p1;
		batch batch_p2;

		for (auto p1_idx : current_cell) {
			for (auto p2_idx : target_cell) {
				batch_p1[count] = p1_idx;
				batch_p2[count] = p2_idx;
				++count;

				if (count == BatchSize) {
					use_batch(batch_p1, batch_p2);
					count = 0;
				}
			}
		}
		use_batch(batch_p1, batch_p2, count);
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

#pragma omp parallel for schedule(dynamic)
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

#pragma omp parallel for simd schedule(static)
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
#pragma omp parallel for simd schedule(static)
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
	force_calculator auto force_calc, particle_container& container,
	const sim_configuration& config, sim_iteration_t iteration
) noexcept {
	if constexpr (!IsFirstIteration) {
		update_values(container.system());
	}

	static const bool has_thermostat = config.thermostat.has_value();
	static const bool has_gravity = config.gravitational_constant != 0;
	calculate_x(container, config);
	calculate_forces_batched<8>(force_calc, container, config);

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

/**
 * @brief Computes the Radial Distribution Function statistic.
 * In case of periodic boundaries the effective shortest distance accounting for mirror particles is used.
 *
 * @param container Particles simulated.
 * @param config Simulation parameters.
 * @param delta_r Sample width for the distance intervals.
 * @return Calculated local densities.
 */
constexpr std::vector<double> radial_distribution_function(particle_container& container, const sim_configuration& config,
	const double delta_r) {
	auto& system = container.system();
	const auto& grid = container.grid_size();
	const double cell_width = container.cutoff_radius();
	const double x_size = cell_width * grid[axis::x];
	const double y_size = cell_width * grid[axis::y];
	const double z_size = cell_width * grid[axis::z];
	const double r_max = std::sqrt(x_size * x_size + y_size * y_size + z_size * z_size);
	const auto number_of_interval = static_cast<std::size_t>(std::ceil(r_max / delta_r));
	std::vector<double> local_densities (number_of_interval, 0);

	for (particle_system::particle_id p1_idx = 0; p1_idx < system.size(); p1_idx++) {
		for (particle_system::particle_id p2_idx = p1_idx + 1; p2_idx < system.size(); p2_idx++) {
			double diff_x = system.x[p1_idx] - system.x[p2_idx];
			if (config.boundary_behavior.x_max == boundary_condition::periodic ||
				config.boundary_behavior.x_min == boundary_condition::periodic ) {
				if (diff_x >  0.5 * x_size) diff_x -= x_size;
				if (diff_x < -0.5 * x_size) diff_x += x_size;
				}
			double diff_y = system.y[p1_idx] - system.y[p2_idx];
			if (config.boundary_behavior.y_max == boundary_condition::periodic ||
				config.boundary_behavior.y_min == boundary_condition::periodic ) {
				if (diff_y >  0.5 * y_size) diff_y -= y_size;
				if (diff_y < -0.5 * y_size) diff_y += y_size;
				}
			double diff_z = system.z[p1_idx] - system.z[p2_idx];
			if (config.boundary_behavior.z_max == boundary_condition::periodic ||
				config.boundary_behavior.z_min == boundary_condition::periodic ) {
				if (diff_z >  0.5 * z_size) diff_z -= z_size;
				if (diff_z < -0.5 * z_size) diff_z += z_size;
			}
			const double r = std::sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
			const auto interval_idx = static_cast<std::size_t>(r/delta_r);
			local_densities[interval_idx] += 1;
		}
	}

	for (unsigned i = 0; i < number_of_interval; i++) {
		const double r_i = i * delta_r;
		const double r_i_next = r_i + delta_r;
		local_densities[i] *= 3/(4 * std::numbers::pi * (r_i_next * r_i_next * r_i_next - r_i * r_i * r_i));
	}
	return local_densities;
}


/**
 * @brief Computes the Diffusion statistic.
 * The absolute distances with respect to periodic boundaries is used.
 *
 * @param container Particles simulated.
 * @param config Simulation parameters.
 * @param old_x X positions at the reference time t0.
 * @param old_y Y positions at the reference time t0.
 * @param old_z Z positions at the reference time t0.
 * @return Calculated variance.
 */
constexpr double diffusion(particle_container& container, const sim_configuration& config,
	std::vector<double>& old_x, std::vector<double>& old_y, std::vector<double>& old_z) {
	double variance = 0;
	auto& system = container.system();
	const auto& grid = container.grid_size();
	const double cell_width = container.cutoff_radius();
	const double x_size = cell_width * grid[axis::x];
	const double y_size = cell_width * grid[axis::y];
	const double z_size = cell_width * grid[axis::z];


	for (particle_system::particle_id p = 0; p< system.size(); p++) {
		if (config.boundary_behavior.x_max == boundary_condition::periodic ||
		    config.boundary_behavior.x_min == boundary_condition::periodic ) {
			double current_x = system.x[p] + x_size * system.x_boundary_crosses[p];
			variance += (current_x - old_x[p]) * (current_x - old_x[p]);
			system.x_boundary_crosses[p] = 0;
		} else {
			variance += (system.x[p] - old_x[p]) * (system.x[p] - old_x[p]);
		}

		if (config.boundary_behavior.y_max == boundary_condition::periodic ||
			config.boundary_behavior.y_min == boundary_condition::periodic ) {
			double current_y = system.y[p] + y_size * system.y_boundary_crosses[p];
			variance += (current_y - old_y[p]) * (current_y - old_y[p]);
			system.y_boundary_crosses[p] = 0;
		} else {
			variance += (system.y[p] - old_y[p]) * (system.y[p] - old_y[p]);
		}

		if (config.boundary_behavior.z_max == boundary_condition::periodic ||
			config.boundary_behavior.z_min == boundary_condition::periodic ) {
			double current_z = system.z[p] + z_size * system.z_boundary_crosses[p];
			variance += (current_z - old_z[p]) * (current_z - old_z[p]);
			system.z_boundary_crosses[p] = 0;
			} else {
				variance += (system.z[p] - old_z[p]) * (system.z[p] - old_z[p]);
			}
	}

	variance = variance / system.size();
	old_x = system.x;
	old_y = system.y;
	old_z = system.z;
	return variance;
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
		"Runtime: {} ms, iterations: {}, tick length: {} ms", elapsed.count(), iteration,
		elapsed.count() / iteration
	);

	if (config.create_checkpoint) {
		// TODO(tuna): specify both in terms of plot_particles, where iteration is used in the
		// filename Or just remove plot_particles
		WRITE_CHECKPOINT(container, fmt::format("{}_checkpoint.json", output_prefix));
	}
}
