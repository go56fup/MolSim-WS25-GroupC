#pragma once

#include <cmath>
#include <mutex>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string_view>

#include <fmt/compile.h>
#include <spdlog/spdlog.h>

#include "grid/bounds/conditions.hpp"
#include "grid/bounds/operations.hpp"
#include "grid/particle_container/fwd.hpp"
#include "grid/particle_container/particle_container.hpp"
#include "grid/particle_container/system.hpp"
#include "iterators/neighbors.hpp"
#include "iterators/pairwise.hpp"
#include "output_writers/io.hpp"
#include "physics/forces.hpp"
#include "simulation/entities.hpp"
#include "simulation/openmp.hpp"
#include "simulation/statistics.hpp"
#include "simulation/thermostat.hpp"
#include "utility/macros.hpp"
#include "utility/tracing/macros.hpp"

/**
 * @brief Calculates interaction forces on a collection of particles.
 * @param container Particles on which forces will be calculated.
 */
constexpr void
calculate_forces_batched(particle_container& container, const sim_configuration& config) noexcept {
	auto use_batch_piecewise = [&](particle_batch batch_p1, particle_batch batch_p2,
	                               std::size_t up_to = batch_size) {
		for (std::size_t i = 0; i < up_to; ++i) {
			// TODO(gabriel): Also we must fully remove Newtons Axiom from our code.
			force_calculator::soa(container.system(), config, batch_p1[i], batch_p2[i]);
		}
	};

	auto use_batch = [&](particle_batch batch_p1, particle_batch batch_p2) {
		force_calculator::batch(container.system(), config, batch_p1, batch_p2);
	};

#if !SINGLETHREADED || OVERRIDE_CALCULATE_F
#pragma omp parallel for schedule(dynamic)
#endif
	for (auto& cell : container.cells()) {
		std::size_t count = 0;
		particle_batch batch_p1;
		particle_batch batch_p2;

		for (auto [p1_idx, p2_idx] : unique_pairs(cell)) {
			TRACE_BATCHING(
				"Putting {} {} into batch, count={} on thread {}", p1_idx, p2_idx, count,
				get_thread_num()

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
#if (!SINGLETHREADED && !DETERMINISTIC) || OVERRIDE_CALCULATE_F
#pragma omp parallel
#endif
	{
#if (!SINGLETHREADED && !DETERMINISTIC) || OVERRIDE_CALCULATE_F
#pragma omp single
#endif
		for (const auto& [current_cell_idx, target_cell_idx] :
		     container.directional_interactions()) {
// TODO(tuna): see if after the implementation of the border iterator whether we still
// need the indices
#if (!SINGLETHREADED && !DETERMINISTIC) || OVERRIDE_CALCULATE_F
#pragma omp task firstprivate(current_cell_idx, target_cell_idx)
#endif
			{
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
		}
	}
}

constexpr void calculate_forces_no_atomic(
	particle_container& container, [[maybe_unused]] const sim_configuration& config
) noexcept {
	const auto& grid = container.grid_size();

	auto use_batch_piecewise = [&](particle_batch batch_p1, particle_batch batch_p2,
	                               std::size_t up_to = batch_size) {
		for (std::size_t i = 0; i < up_to; ++i) {
			lennard_jones_force_soa_no_atomic(container.system(), batch_p1[i], batch_p2[i], {});
		}
	};

	auto use_batch = [&](particle_batch batch_p1, particle_batch batch_p2) {
		lennard_jones_force_soa_batchwise_no_atomic(container.system(), batch_p1, batch_p2, {});
	};

#if !SINGLETHREADED || OVERRIDE_CALCULATE_F
#pragma omp parallel for schedule(dynamic)
#endif
	for (auto& cell : container.cells()) {
		std::size_t count = 0;
		particle_batch batch_p1;
		particle_batch batch_p2;
		for (auto p1_idx : cell) {
			for (auto p2_idx : cell) {
				if (p1_idx == p2_idx) continue;
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
	}

#if (!SINGLETHREADED && !DETERMINISTIC) || OVERRIDE_CALCULATE_F
#pragma omp parallel
#endif
	{
#if (!SINGLETHREADED && !DETERMINISTIC) || OVERRIDE_CALCULATE_F
#pragma omp single
#endif
		for (particle_container::index current_cell_idx{0, 0, 0}; current_cell_idx.z < grid.z;
		     current_cell_idx = next_3d_index(current_cell_idx, grid)) {
			for (const auto& target_cell_idx : neighbors_range(container, current_cell_idx)) {
#if (!SINGLETHREADED && !DETERMINISTIC) || OVERRIDE_CALCULATE_F
#pragma omp task firstprivate(current_cell_idx, target_cell_idx)
#endif
				{
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
			}
		}
	}
}

/**
 * @brief Updates the position of each particle and move particles inside grid.
 *
 * Implements the first step of the velocity Störmer-Verlet algorithm:
 * \f[
 * \vec{x}(t + \Delta t) = \vec{x}(t) + \Delta t \cdot \vec{v}(t) + \frac{\Delta t^2}{2m}
 * \vec{F}(t)
 * \f]
 *
 * @param container Particles to calculate positional updates for.
 * @param config Simulation parameters.
 */
constexpr void
calculate_x(particle_container& container, const sim_configuration& config) noexcept(false) {
	const double cell_width = container.cutoff_radius();
	STATIC_IF_NOT_TESTING const bool collect_stats = config.statistics.has_value();

	auto move_to_cell = [&](const particle_container::index& new_cell_idx,
	                        particle_container::cell& current_cell,
	                        particle_container::cell::size_type i) {
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

#if !SINGLETHREADED
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
							"Performing periodic teleport on {} at {}, currently stored in {}",
							p_idx, system.serialize_position(p_idx), cell_idx
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
						if (collect_stats) {
							system.crossing_statistics_component(axis_)[p_idx] += increment;
							TRACE_STATS(
								"Logged crossing for {}: {} on {}", p_idx, increment, axis_
							);
						}

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
			bool was_hit = false;

			if (system.x[p_idx] < cell_origin.x) {
				update_displacement.template operator()<x_min>();
				was_hit = true;
			} else if (system.x[p_idx] >= cell_origin.x + cell_width) {
				update_displacement.template operator()<x_max>();
				was_hit = true;
			}

			if (system.y[p_idx] < cell_origin.y) {
				update_displacement.template operator()<y_min>();
				was_hit = true;
			} else if (system.y[p_idx] >= cell_origin.y + cell_width) {
				update_displacement.template operator()<y_max>();
				was_hit = true;
			}

			if (system.z[p_idx] < cell_origin.z) {
				update_displacement.template operator()<z_min>();
				was_hit = true;
			} else if (system.z[p_idx] >= cell_origin.z + cell_width) {
				update_displacement.template operator()<z_max>();
				was_hit = true;
			}
			if (was_hit) {
#if !SINGLETHREADED
				// #pragma omp critical
#endif
				move_to_cell(new_cell_idx, cell, particle_idx);
			}
#if LOG_GRID
			const auto pos = system.serialize_position(p_idx);
			const auto& domain = container.domain();
			if (out_of_bounds<boundary_type::x_min>(pos, domain) ||
			    out_of_bounds<boundary_type::x_max>(pos, domain) ||
			    out_of_bounds<boundary_type::y_min>(pos, domain) ||
			    out_of_bounds<boundary_type::y_max>(pos, domain) ||
			    out_of_bounds<boundary_type::z_min>(pos, domain) ||
			    out_of_bounds<boundary_type::z_max>(pos, domain)) {
				throw std::out_of_range(fmt::format(
					"Trying to move {} with v={}, f={}, old_f={} to position {}, which "
					"is out of "
					"bounds for domain {}.",
					p_idx, system.serialize_velocity(p_idx), system.serialize_force(p_idx),
					system.serialize_old_force(p_idx), pos, domain
				));
			}
#endif
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
 * @param system Particle data to update.
 * @param delta_t The time step for integration.
 */
constexpr void calculate_v(particle_system& system, double delta_t) {
#if !SINGLETHREADED
#pragma omp parallel for simd schedule(static)
#endif
	for (particle_id i = 0; i < system.size(); ++i) {
		// TODO(gabriel): Maybe add a invM array and * 0.5 instead of / 2 to remove all these costly
		// divisions
		const auto velocity_scalar = delta_t / (2 * system.mass[i]);
		system.vx[i] += velocity_scalar * (system.old_fx[i] + system.fx[i]);
		system.vy[i] += velocity_scalar * (system.old_fy[i] + system.fy[i]);
		system.vz[i] += velocity_scalar * (system.old_fz[i] + system.fz[i]);
	}
#if LOG_GRID
	static constexpr double thresh = 100;
	for (particle_id i = 0; i < system.size(); ++i) {
		if (system.vx[i] > thresh || system.vy[i] > thresh || system.vz[i] > thresh) {
			throw std::logic_error(fmt::format("Illegal velocity: {}", system.representation(i)));
		}
	}
#endif
}

/**
 * @brief Prepares particles for the next iteration.
 *
 * Sets the old force of each particle to the force calculated within the
 * most recent iteration at the end of each iteraton.
 *
 * @param system Particle data to update.
 */
constexpr void update_values(particle_system& system) noexcept {
#if !SINGLETHREADED
#pragma omp parallel for simd schedule(static)
#endif
	for (particle_id p = 0; p < system.size(); ++p) {
		system.old_fx[p] = system.fx[p];
		system.old_fy[p] = system.fy[p];
		system.old_fz[p] = system.fz[p];
		system.fx[p] = 0;
		system.fy[p] = 0;
		system.fz[p] = 0;
	}
}

constexpr void calculate_membrane_forces(
	particle_container& container, const sim_configuration& config, double current_time
) {
	static constexpr double upwards_force_until_timepoint = 150;
	// STATIC_IF_NOT_TESTING bool still_upwards = true;
	const auto& membrane_particles = container.membrane_particles();

	const auto& scale = container.membrane_scale();
	const auto& width = scale.y;
	const auto& height = scale.x;

	using enum direction;
	auto perform_interaction = [&]<direction d>(std::size_t current_idx) {
		static constexpr bool is_orthogonal = d == down || d == up || d == left || d == right;
		particle_container::difference_type offset = 0;
		if constexpr ((d & left) == left) {
			offset -= 1;
		}
		if constexpr ((d & right) == right) {
			offset += 1;
		}
		if constexpr ((d & up) == up) {
			offset += static_cast<decltype(offset)>(width);
		}
		if constexpr ((d & down) == down) {
			offset -= static_cast<decltype(offset)>(width);
		}
		const auto other_idx = apply_difference(current_idx, offset);
		TRACE_MEMBRANE_NEIGHBORS(
			"Performing interaction towards {}; offset={}. Got current={}, target={} -> [{}, {}] "
			"for scale={}",
			d, offset, current_idx, other_idx, other_idx / width, other_idx % width, scale
		);

		const particle_id other = membrane_particles[other_idx];
		const particle_id current = membrane_particles[current_idx];
		if constexpr (is_orthogonal) {
			harmonic_force_orthogonal(
				container.system(), current, other, *config.membrane_parameters
			);
		} else {
			harmonic_force_diagonal(
				container.system(), current, other, *config.membrane_parameters
			);
		}
	};
	auto perform_horizontal = [&]<direction d>(std::size_t idx, bool has_up, bool has_down) {
		perform_interaction.template operator()<d>(idx);
		if (has_down) {
			perform_interaction.template operator()<d | down>(idx);
		}
		if (has_up) {
			perform_interaction.template operator()<d | up>(idx);
		}
	};

	for (std::size_t idx = 0; idx < membrane_particles.size(); ++idx) {
		const std::size_t x = idx / width;
		const std::size_t y = idx % width;

		const bool has_left = (y > 0);
		const bool has_right = (y + 1 < width);
		const bool has_down = (x > 0);
		const bool has_up = (x + 1 < height);

		TRACE_MEMBRANE_NEIGHBORS(
			"idx={} -> [{}, {}] with scale={}: left={}, right={}, down={}, up={}", idx, x, y, scale,
			has_left, has_right, has_down, has_up
		);

		if (has_left) {
			perform_horizontal.template operator()<left>(idx, has_up, has_down);
		}
		if (has_right) {
			perform_horizontal.template operator()<right>(idx, has_up, has_down);
		}
		if (has_up) {
			perform_interaction.template operator()<up>(idx);
		}
		if (has_down) {
			perform_interaction.template operator()<down>(idx);
		}
	}
	if (current_time < upwards_force_until_timepoint) {
		for (auto id : container.upwards_moving_membrane_members()) {
			container.system().fz[id] += config.membrane_parameters->upwards_force;
		}
	}
#if 0
	if (still_upwards) {
		if (current_time < upwards_force_until_timepoint) {
			for (auto id : container.upwards_moving_membrane_members()) {
				container.system().fz[id] += config.membrane_parameters->upwards_force;
			}
		} else {
			still_upwards = false;
		}
	}
#endif
}

/**
 * @brief Run a tick of the simulation.
 *
 * @param container Particles to simulate.
 * @param config Simulation parameters.
 * @param iteration The index of iteration to simulate.
 */
// TODO(tuna): now that iteration is being passed in, isFirstIteration should probably be
// removed
template <bool IsFirstIteration = false>
constexpr void run_sim_iteration(
	particle_container& container, const sim_configuration& config, sim_iteration_t iteration,
	double current_time
) {
	if constexpr (!IsFirstIteration) {
		update_values(container.system());
	}

	STATIC_IF_NOT_TESTING const bool has_thermostat = config.thermostat.has_value();
	STATIC_IF_NOT_TESTING const bool has_gravity = config.gravitational_constant != 0;
	STATIC_IF_NOT_TESTING const bool has_membrane = config.membrane_parameters.has_value();

	calculate_x(container, config);
	calculate_forces_batched(container, config);

	if (has_membrane) {
		calculate_membrane_forces(container, config, current_time);
	}

	handle_boundaries(container, config);

	if (has_thermostat && iteration % config.thermostat->application_frequency == 0) {
		run_thermostat(container, *config.thermostat, config.dimensions);
	}

	if (has_gravity) {
		if (has_membrane) {
			apply_gravity<axis::z>(container, config.gravitational_constant);
		} else {
			apply_gravity<axis::y>(container, config.gravitational_constant);
		}
	}

	calculate_v(container.system(), config.delta_t);

#if LOG_SIM_STATE
	const auto& system = container.system();
	for ([[maybe_unused]] particle_id i = 0; i < system.size(); ++i) {
		TRACE_SIM_STATE(
			"End of iteration, particle {}: pos={}, v={}, f={}, old_f={}", i,
			system.serialize_position(i), system.serialize_velocity(i), system.serialize_force(i),
			system.serialize_old_force(i)
		);
	}
#endif
}

/**
 * @brief Start the simulation.
 *
 * Calculates new positions, velocities and forces at each time tick
 * and plots the updated particles.
 *
 * @param container Particles to simulate.
 * @param config Simulation parameters.
 * @param output_prefix Path to put simulation result files into.
 * @return Number of iterations taken to reach specified end time.
 **/
constexpr sim_iteration_t run_simulation(
	particle_container& container, const sim_configuration& config,
	[[maybe_unused]] std::string_view output_prefix
) {
	double current_time = 0;
	sim_iteration_t iteration = 0;

	struct sim_traits {
		bool is_first_iteration = false;
	};

	auto run_sim_pass = [&]<sim_traits Traits = {}> {
		STATIC_IF_NOT_TESTING const bool collect_stats = config.statistics.has_value();
		TRACE_SIM("beginning iteration {}, current_time={}", iteration, current_time);
		if (iteration % config.write_frequency == 0) {
			WRITE_VTK_OUTPUT(container, output_prefix, iteration);
		}

		run_sim_iteration<Traits.is_first_iteration>(container, config, iteration, current_time);

		if (collect_stats && iteration % config.statistics->calculation_frequency == 0) {
			WRITE_STATISTICS_OUTPUT(container, config, output_prefix, iteration);
		}

		current_time += config.delta_t;
		++iteration;
	};

	run_sim_pass.template operator()<{.is_first_iteration = true}>();
	while (current_time < config.end_time) {
		run_sim_pass.template operator()();
	}
	return iteration;
}
