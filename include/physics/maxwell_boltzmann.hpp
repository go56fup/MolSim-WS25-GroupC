/*
 * MaxwellBoltzmannDistribution.h
 *
 * @Date: 13.12.2019
 * @Author: F. Gratl
 */

#pragma once

#include <random>

#include "physics/vec_3d.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/constants.hpp"

namespace detail {
#if IS_GCC && !HAS_EMBED
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
#pragma GCC diagnostic ignored "-Woverflow"
#endif
inline constexpr auto random_numbers = [] consteval {
	// C arrays are how #embed is meant to be used.
	// NOLINTNEXTLINE(*avoid-c-arrays)
	const std::uint8_t data[]{
#if HAS_EMBED
#embed "random.bin"
#else
#include "utility/random_preprocessed.hpp"
#endif
	};
	return std::bit_cast<std::array<double, random_table_size>>(data);
}();
#if IS_GCC
#pragma GCC diagnostic pop
#endif

inline double get_random_runtime() {
	// we use a constant seed for repeatability.
	static constexpr auto seed = 42;
	// The seed being constant is on purpose.
	// random engine needs static lifetime otherwise it would be recreated for every call.
	// NOLINTNEXTLINE(cert-msc*)
	static std::default_random_engine randomEngine(seed);

	// when adding independent normally distributed values to all velocity components
	// the velocity change is maxwell boltzmann distributed
	std::normal_distribution<double> normalDistribution{0, 1};
	return normalDistribution(randomEngine);
}

constexpr double get_random([[maybe_unused]] std::size_t no) {
	if not consteval {
		return get_random_runtime();
	} else {
		return random_numbers[no];
	}
}

}  // namespace detail

/**
 * Generate a random velocity vector according to the Maxwell-Boltzmann distribution, with a given
 * average velocity.
 *
 * @param average_velocity The average velocity of the brownian motion for the system.
 * @param seq_no Mutable reference to a sequence number, used to achieve stateless random number
 * generation.
 * @tparam Dimensions Number of dimensions for which the velocity vector shall be generated, can be
 * 2 or 3.
 * @return The generated velocity vector.
 */
template <std::size_t Dimensions>
	requires(Dimensions == 2 || Dimensions == 3)
constexpr vec maxwell_boltzmann_distributed_velocity(double average_velocity, std::size_t& seq_no) {
	vec randomVelocity;

	randomVelocity.x = average_velocity * detail::get_random(seq_no++);
	randomVelocity.y = average_velocity * detail::get_random(seq_no++);

	if constexpr (Dimensions == 3) {
		randomVelocity.z = average_velocity * detail::get_random(seq_no++);
	}

	return randomVelocity;
}
