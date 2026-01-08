#pragma once

#include <algorithm>
#include <cmath>

#include "grid/particle_container/fwd.hpp"
#include "physics/particle.hpp"
#include "simulation/config/entities.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/concepts.hpp"

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
	double current_temp = 0;
	for (particle& p : container.particles()) {
		const double squared_norm = (p.v.x * p.v.x) + (p.v.y * p.v.y) + (p.v.z * p.v.z);
		current_temp += p.m * squared_norm;
	}
	current_temp /= dimensions * static_cast<double>(container.size());

	const double& target = therm.target_temperature;
	const double& max_change = therm.max_temperature_difference;

	// I am fairly sure that these two are equivalent, but I'm leaving it on the original
	// implementation (by not defining USE_CLAMP) just in case.

	// #define USE_CLAMP
#if USE_CLAMP
	const double new_temp =
		std::clamp(target, current_temp - max_change, current_temp + max_change);
#else
	double new_temp = target;
	if (target - current_temp > max_change) {
		new_temp = current_temp + max_change;
	} else if (current_temp - target > max_change) {
		new_temp = current_temp - max_change;
	}
#endif
	const double beta = std::sqrt(new_temp / current_temp);
	for (particle& p : container.particles()) {
		p.v *= beta;
	}
}
