#pragma once

#include <concepts>
#include <span>
#include <string_view>

#include "Particle.h"

/**
 * @brief Concept that defines a particle I/O provider.
 *
 * A callable satisfies this concept if it can be invoked as-if it shares this signature:
 * @verbatim
 void operator()(std::span<const Particle> particles, std::string_view base_name, int iteration)
 @endverbatim
 * where @p particles is the particles to output, @p base_name is the base name for the output files being
 * generated and @p iteration is the current simulation iteration.
 *
 * @tparam Candidate Type being tested for the particle I/O provider interface.
 *
 * @return `true` @a iff the callable can be invoked as specified.
 */
template <typename Candidate>
concept particle_io_provider =
	requires(Candidate f, std::span<const Particle> particles, std::string_view out_name, int iteration) {
		{ std::invoke(f, particles, out_name, iteration) } -> std::same_as<void>;
	};
