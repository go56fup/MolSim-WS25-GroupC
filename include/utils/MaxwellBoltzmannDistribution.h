/*
 * MaxwellBoltzmannDistribution.h
 *
 * @Date: 13.12.2019
 * @Author: F. Gratl
 */

#pragma once

#include <random>

#include "Vector.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
inline constexpr auto numbers = [] consteval {
	// Currently, the two implementations diverge after 333 particles.
	// TODO(tuna): see if this needs to be increased
	const std::uint8_t data[] {
#embed "random.bin"
	};
	return std::bit_cast<std::array<double, 1000>>(data);
}();
#pragma GCC diagnostic pop

constexpr double get_random(std::size_t no) {
	if not consteval {
		// we use a constant seed for repeatability.
		static constexpr auto seed = 42;
		// The seed being constant is on purpose.
		// NOLINTNEXTLINE(cert-msc*)
		// random engine needs static lifetime otherwise it would be recreated for every call.
		static std::default_random_engine randomEngine(seed);

		// when adding independent normally distributed values to all velocity components
		// the velocity change is maxwell boltzmann distributed
		std::normal_distribution<double> normalDistribution{0, 1};
		return normalDistribution(randomEngine);
	} else {
		return numbers[no];
	}
}

/**
 * Generate a random velocity vector according to the Maxwell-Boltzmann distribution, with a given average velocity.
 *
 * @param averageVelocity The average velocity of the brownian motion for the system.
 * @tparam Dimensions Number of dimensions for which the velocity vector shall be generated, can be 2 or 3.
 * @return The generated velocity vector.
 */
template <std::size_t Dimensions>
	requires(Dimensions == 2 || Dimensions == 3)
constexpr vec maxwellBoltzmannDistributedVelocity(double averageVelocity, std::size_t& seq_no) {
	vec randomVelocity;

	randomVelocity.x = averageVelocity * get_random(seq_no++);
	randomVelocity.y = averageVelocity * get_random(seq_no++);

	if constexpr (Dimensions == 3) {
		randomVelocity.z = averageVelocity * get_random(seq_no++);
	}

	return randomVelocity;
}
