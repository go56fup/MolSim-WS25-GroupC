#pragma once

#include <algorithm>
#include <cmath>

#include "grid/particle_container/particle_container.hpp"
#include "physics/particle.hpp"
#include "simulation/config/entities.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/concepts.hpp"

CONSTEXPR_IF_GCC inline double
get_temperature(particle_container& container, decltype(sim_configuration::dimensions) dimensions) {
	double result = 0;
	const auto& system = container.system();
	#pragma omp parallel for simd schedule(static) reduction(+ : result)
	for (particle_system::particle_id p = 0; p < system.size(); ++p) {
		const double squared_norm = (system.vx[p] * system.vx[p]) + (system.vy[p] * system.vy[p]) +
		                            (system.vz[p] * system.vz[p]);
		result += container.material_for_particle(p).mass * squared_norm; //TODO(gabriel): make mass vectorisable
	}
	result /= dimensions * static_cast<double>(system.size());
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
CONSTEXPR_IF_GCC inline void run_thermostat(
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
		auto& system = container.system();

		#pragma omp parallel for simd schedule(static)
		for (particle_system::particle_id p = 0; p < system.size(); ++p) {
			system.vx[p] *= beta;
			system.vy[p] *= beta;
			system.vz[p] *= beta;
		}
	}
}
