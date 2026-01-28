#pragma once

#include <algorithm>
#include <cmath>

#include "grid/particle_container/particle_container.hpp"
#include "physics/particle.hpp"
#include "simulation/config/entities.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/concepts.hpp"

CONSTEXPR_IF_GCC double
get_temperature(particle_container& container, decltype(sim_configuration::dimensions) dimensions) {
	double result = 0;
	for (particle& p : container.particles()) {
		const double squared_norm = (p.v.x * p.v.x) + (p.v.y * p.v.y) + (p.v.z * p.v.z);
		result += p.m * squared_norm;
	}
	result /= dimensions * static_cast<double>(container.size());
	return result;
}

/**
 * @brief One thermostat iteration. Changes the temperature by scaling the velocity of particles.
 *
 * @param particles Particles of which to alter the temperature of.
 * @param target_temp Target temperature.
 * @param max_change Maximum temperature change for one application. Controls gradual/direct
 * velocity scaling.
 * @param dimensions Degrees of freedom in the simulation.
 */
CONSTEXPR_IF_GCC void run_thermostat(
	particle_container& container, const thermostat_parameters& therm,
	decltype(sim_configuration::dimensions) dimensions
) {
	const double current_temp = get_temperature(container, dimensions);

	const double& target = therm.target_temperature;
	const double& max_change = therm.max_temperature_difference;

	const double new_temp =
		std::clamp(target, current_temp - max_change, current_temp + max_change);
	TRACE_THERMOSTAT("Current temperature: {},  new temperature: {}", current_temp, new_temp);
	if (const double beta = std::sqrt(new_temp / current_temp); beta != 1) {
		for (particle& p : container.particles()) {
			p.v *= beta;
		}
	}
}
