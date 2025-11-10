/*
 * MaxwellBoltzmannDistribution.h
 *
 * @Date: 13.12.2019
 * @Author: F. Gratl
 */

#pragma once

#include <random>

#include "Vector.h"

/**
 * Generate a random velocity vector according to the Maxwell-Boltzmann distribution, with a given average velocity.
 *
 * @param averageVelocity The average velocity of the brownian motion for the system.
 * @tparam Dimensions Number of dimensions for which the velocity vector shall be generated, can be 2 or 3.
 * @return The generated velocity vector.
 */
template <std::size_t Dimensions>
	requires(Dimensions == 2 || Dimensions == 3)
constexpr vec maxwellBoltzmannDistributedVelocity(double averageVelocity) {
	// we use a constant seed for repeatability.
	// random engine needs static lifetime otherwise it would be recreated for every call.
	static constexpr auto seed = 42;
	// The seed being constant is on purpose.
	// NOLINTNEXTLINE(cert-msc*)
	static std::default_random_engine randomEngine(seed);

	// when adding independent normally distributed values to all velocity components
	// the velocity change is maxwell boltzmann distributed
	std::normal_distribution<double> normalDistribution{0, 1};
	vec randomVelocity;

	randomVelocity.x = averageVelocity * normalDistribution(randomEngine);
	randomVelocity.y = averageVelocity * normalDistribution(randomEngine);

	if constexpr (Dimensions == 3) {
		randomVelocity.z = averageVelocity * normalDistribution(randomEngine);
	}

	return randomVelocity;
}
