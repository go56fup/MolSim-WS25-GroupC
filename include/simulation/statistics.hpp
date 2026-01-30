#pragma once

#include <numbers>
#include <vector>

#include "grid/particle_container/fwd.hpp"
#include "simulation/entities.hpp"

/**
 * @brief Computes the Radial Distribution Function statistic.
 * In case of periodic boundaries the effective shortest distance accounting for mirror particles is
 * used.
 *
 * @param container Particles simulated.
 * @param config Simulation parameters.
 * @param delta_r Sample width for the distance intervals.
 * @return Calculated local densities.
 */
constexpr std::vector<double>
radial_distribution_function(particle_container& container, const sim_configuration& config) {
	auto& system = container.system();
	const auto& grid = container.grid_size();
	const double cell_width = container.cutoff_radius();
	const double delta_r = config.statistics->delta_r;
	const double x_size = cell_width * grid.x;
	const double y_size = cell_width * grid.y;
	const double z_size = cell_width * grid.z;
	const double r_max = std::sqrt(x_size * x_size + y_size * y_size + z_size * z_size);
	const auto number_of_intervals = static_cast<std::size_t>(std::ceil(r_max / delta_r));
	std::vector<double> local_densities(number_of_intervals, 0);

	for (particle_id p1_idx = 0; p1_idx < system.size(); p1_idx++) {
		for (particle_id p2_idx = p1_idx + 1; p2_idx < system.size(); p2_idx++) {
			double diff_x = system.x[p1_idx] - system.x[p2_idx];
			if (config.boundary_behavior.x_min == boundary_condition::periodic) {
				if (diff_x > 0.5 * x_size) diff_x -= x_size;
				if (diff_x < -0.5 * x_size) diff_x += x_size;
			}
			double diff_y = system.y[p1_idx] - system.y[p2_idx];
			if (config.boundary_behavior.y_min == boundary_condition::periodic) {
				if (diff_y > 0.5 * y_size) diff_y -= y_size;
				if (diff_y < -0.5 * y_size) diff_y += y_size;
			}
			double diff_z = system.z[p1_idx] - system.z[p2_idx];
			if (config.boundary_behavior.z_min == boundary_condition::periodic) {
				if (diff_z > 0.5 * z_size) diff_z -= z_size;
				if (diff_z < -0.5 * z_size) diff_z += z_size;
			}
			const double r = std::sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
			const auto interval_idx = static_cast<std::size_t>(r / delta_r);
			local_densities[interval_idx] += 1;
		}
	}

	for (unsigned i = 0; i < number_of_intervals; i++) {
		const double r_i = i * delta_r;
		const double r_i_next = r_i + delta_r;
		local_densities[i] *=
			3 / (4 * std::numbers::pi * (r_i_next * r_i_next * r_i_next - r_i * r_i * r_i));
	}
	return local_densities;
}

/**
 * @brief Computes the Diffusion statistic.
 * The absolute distances with respect to periodic boundaries is used.
 *
 * @param container Particles simulated.
 * @param config Simulation parameters.
 * @return Calculated variance.
 */
constexpr double diffusion(particle_container& container, const sim_configuration& config) {
	double variance = 0;
	auto& system = container.system();
	const auto& grid = container.grid_size();
	const double cell_width = container.cutoff_radius();
	const double x_size = cell_width * grid.x;
	const double y_size = cell_width * grid.y;
	const double z_size = cell_width * grid.z;

	for (particle_id p = 0; p < system.size(); p++) {
		if (config.boundary_behavior.x_min == boundary_condition::periodic) {
			double current_x =
				system.x[p] + x_size * static_cast<double>(system.x_boundary_crosses[p]);
			variance += (current_x - system.old_x[p]) * (current_x - system.old_x[p]);
			system.x_boundary_crosses[p] = 0;
		} else {
			variance += (system.x[p] - system.old_x[p]) * (system.x[p] - system.old_x[p]);
		}

		if (config.boundary_behavior.y_min == boundary_condition::periodic) {
			double current_y =
				system.y[p] + y_size * static_cast<double>(system.y_boundary_crosses[p]);
			variance += (current_y - system.old_y[p]) * (current_y - system.old_y[p]);
			system.y_boundary_crosses[p] = 0;
		} else {
			variance += (system.y[p] - system.old_y[p]) * (system.y[p] - system.old_y[p]);
		}

		if (config.boundary_behavior.z_min == boundary_condition::periodic) {
			double current_z =
				system.z[p] + z_size * static_cast<double>(system.z_boundary_crosses[p]);
			variance += (current_z - system.old_z[p]) * (current_z - system.old_z[p]);
			system.z_boundary_crosses[p] = 0;
		} else {
			variance += (system.z[p] - system.old_z[p]) * (system.z[p] - system.old_z[p]);
		}
	}

	variance /= static_cast<double>(system.size());
	system.old_x = system.x;
	system.old_y = system.y;
	system.old_z = system.z;
	return variance;
}
