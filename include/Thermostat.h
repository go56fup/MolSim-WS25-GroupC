#pragma once
#include "Particle.h"
#include "Concepts.h"
#include "CompilerTraits.h"

/**
 * @brief One thermostat iteration. Changes the temperature by scaling the velocity of particles.
 *
 * @param particles Particles of which to alter the temperature of.
 * @param target_temp Target temperature.
 * @param max_change Maximum temperature change for one application. Controls gradual/direct velocity scaling.
 * @param dimensions Degrees of freedom in the simulation.
 */
CONSTEXPR_IF_GCC void update_thermostat(range_of<Particle> auto&& particles, const double target_temp, const double max_change, const int dimensions) {
	double current_temp = 0;
	for(Particle& particle : particles){
		current_temp += particle.m * particle.v.squared_norm();
	}
	current_temp /= dimensions * std::ranges::size(particles);

	double new_temp = target_temp;
	if (target_temp - current_temp > max_change) {
		new_temp = current_temp + max_change;
	}else if (current_temp - target_temp > max_change) {
		new_temp = current_temp - max_change;
	}

	const double beta = std::sqrt(new_temp/current_temp);
	for (Particle& particle: particles) {
		particle.v *= beta;
	}
}

