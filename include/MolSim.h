#pragma once

#include <cmath>
#include <functional>
#include <ranges>
#include <span>
#include <string_view>

#include "Concepts.h"
#include "ForceCalculators.h"
#include "IOProviders.h"
#include "outputWriter/outputWriters.h"
#include "ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"

constexpr void update_particle_forces(force_calculator auto calculator, Particle& p1, Particle& p2) noexcept {
	const auto f_ij = std::invoke(calculator, p1, p2);
	p1.f += f_ij;
	p2.f -= f_ij;
}

constexpr void if_cell_at_border_reflect_via_ghost_particle(
	ParticleContainer::cell& current_cell, const index_3d& current_cell_idx, const index_3d& domain,
	force_calculator auto calculator
) noexcept {
	static constexpr double ghost_particle_threshold = 1;
	const auto& dom_x = static_cast<double>(domain.x);
	const auto& dom_y = static_cast<double>(domain.y);
	const auto& dom_z = static_cast<double>(domain.z);

	if (current_cell_idx.x == 0) {
		for (auto& p : current_cell) {
			const auto& pos = p.x;
			if (pos.x <= ghost_particle_threshold) {
				Particle ghost_particle({-pos.x, pos.y, pos.z}, {}, p.m);
				p.f += std::invoke(calculator, p, ghost_particle);
			}
		}

	} else if (current_cell_idx.y == 0) {
		for (auto& p : current_cell) {
			const auto& pos = p.x;
			if (pos.y <= ghost_particle_threshold) {
				Particle ghost_particle({pos.x, -pos.y, pos.z}, {}, p.m);
				p.f += std::invoke(calculator, p, ghost_particle);
			}
		}

	} else if (current_cell_idx.z == 0) {
		for (auto& p : current_cell) {
			const auto& pos = p.x;
			if (pos.z <= ghost_particle_threshold) {
				Particle ghost_particle({pos.x, pos.y, -pos.z}, {}, p.m);
				p.f += std::invoke(calculator, p, ghost_particle);
			}
		}

	} else if (current_cell_idx.x == domain.x) {
		for (auto& p : current_cell) {
			const auto& pos = p.x;
			if (pos.x >= dom_x - ghost_particle_threshold) {
				Particle ghost_particle({dom_x + pos.x, pos.y, pos.z}, {}, p.m);
				p.f += std::invoke(calculator, p, ghost_particle);
			}
		}
	} else if (current_cell_idx.y == domain.y) {
		for (auto& p : current_cell) {
			const auto& pos = p.x;
			if (pos.y >= dom_y - ghost_particle_threshold) {
				Particle ghost_particle({pos.x, dom_y + pos.y, pos.z}, {}, p.m);
				p.f += std::invoke(calculator, p, ghost_particle);
			}
		}
	} else if (current_cell_idx.z == domain.z) {
		for (auto& p : current_cell) {
			const auto& pos = p.x;
			if (pos.z >= dom_z - ghost_particle_threshold) {
				Particle ghost_particle({pos.x, pos.y, dom_z + pos.z}, {}, p.m);
				p.f += std::invoke(calculator, p, ghost_particle);
			}
		}
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
	for (auto&& [current_cell_idx, target_cell_idx] : particles.directional_interactions()) {
		const auto& current_cell = particles[current_cell_idx];
		const auto& target_cell = particles[target_cell_idx];

		for (auto&& [p1, p2] : unique_pairs(current_cell)) {
			update_particle_forces(calculator, p1, p2);
		}
		for (const auto& p1 : current_cell) {
			for (const auto& p2 : target_cell) {
				update_particle_forces(calculator, p1, p2);
			}
		}
		// TODO(tuna): there must be a more efficient way to check this rather than checking for every cell
		if_cell_at_border_reflect_via_ghost_particle(current_cell_idx, particles.domain(), calculator);
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

	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int k = 0;

	const double cell_width = 0;

	for (auto& cell : particles) {

		//for j (iterator goes from top to bottom while coordinate system goes from bottom up)
		const vec cell_origin = {i * cell_width, j * cell_width, k * cell_width};

		//out of cell check
		for (const auto& particle : cell){

			const auto force_scalar = std::pow(delta_t, 2) / (2 * particle.m);
			// TODO(tuna): check if expression templates are needed here
			particle.x = particle.x + delta_t * particle.v + particle.old_f * force_scalar;

			if (particle.x.x < cell_origin.x) {  //for x direction
                ParticleContainer::cell& new_cell = particles[i - 1, j, k];
				new_cell.push_back(particle);
                //TODO(tuna) :ADD delete()

			} else if (particle.x.x >= cell_origin.x + cell_width) {
				ParticleContainer::cell& new_cell = particles[i + 1, j, k];
				new_cell.push_back(particle);
				//TODO(tuna) :ADD delete()

			} else if (particle.x.y < cell_origin.y) { //for y direction
				ParticleContainer::cell& new_cell = particles[i , j - 1, k];
				new_cell.push_back(particle);
				//TODO(tuna) :ADD delete()

			} else if (particle.x.y >= cell_origin.y + cell_width) {
				ParticleContainer::cell& new_cell = particles[i, j + 1, k];
				new_cell.push_back(particle);
				//TODO(tuna) :ADD delete()

			} else if (particle.x.z < cell_origin.z) { //for z direction
				ParticleContainer::cell& new_cell = particles[i , j, k - 1];
				new_cell.push_back(particle);
				//TODO(tuna) :ADD delete()

			} else if (particle.x.z >= cell_origin.z + cell_width) {
				ParticleContainer::cell& new_cell = particles[i, j , k + 1];
				new_cell.push_back(particle);
				//TODO(tuna) :ADD delete()

			}
		}

		//update cell coordinate (normalised)
		i++;
		if (i == particles.domain().x) {
			i = 0;
			j++;
		}
		if (j == particles.domain().y) {
			j = 0;
			k++;
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
	particle_io_provider auto&& io_provider, range_of<const Particle> auto&& particles, std::string_view out_name,
	int iteration
) {
	io_provider(particles, out_name, iteration);
}

/**
 * @brief Prepares particles for the next iteration.
 *
 * Currently, it only sets the old force of each particle to the force calculated within the
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
template <force_calculator ForceCalculator>
constexpr void run_sim_iteration(ForceCalculator&& force_calc, ParticleContainer& particles, double delta_t) noexcept {
	calculateX(particles.view(), delta_t);
	calculateF(std::forward<ForceCalculator>(force_calc), particles);
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
				// TODO(tuna): make the 2 not a template parameter and incorporate into input file format
				particles.emplace(
					vec{origin.x + (i * distance), origin.y + (j * distance), origin.z + (k * distance)},
					maxwellBoltzmannDistributedVelocity<2>(brownian_mean) + initial_velocity, mass
				);
			}
		}
	}
}
