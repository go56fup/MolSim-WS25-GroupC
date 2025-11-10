#pragma once

#include <cmath>
#include <concepts>
#include <functional>
#include <span>

#include "Particle.h"
#include "CompilerTraits.h"

/**
 * @brief Concept that defines a force calculator.
 *
 * A callable satisfies this concept if it can be invoked as-if it shares this signature:
 * @verbatim void operator()(std::span<Particle> particles) @endverbatim
 * where @p particles is the span of particles for which the forces should be calculated and updated.
 *
 * @tparam Candidate Type being tested for the force calculator interface.
 *
 * @return `true` @a iff the callable can be invoked as specified.
 */
template <typename Candidate>
concept force_calculator = requires(Candidate f, std::span<Particle> ps) {
	{ std::invoke(f, ps) } -> std::same_as<void>;
};

namespace detail {
/**
 * @brief Computes the force component exerted on one particle by another.
 *
 * Calculates the gravitational-like force acting on particle @p p1 by particle @p p2
 * according to the inverse-cube law:
 * \f[
 * \vec{F}_{ij} = G \frac{m_i m_j}{r_{ij}^3} (\vec{x}_j - \vec{x}_i)
 * \f]
 * (where G is taken as 1 for simplicity).
 *
 * @param p1 The particle being acted upon (particle i).
 * @param p2 The particle exerting the force (particle j).
 * @return A vector representing the force acting on @p p1 from @p p2.
 */
CONSTEXPR_IF_GCC inline vec calculate_component(const Particle& p1, const Particle& p2) noexcept {
	const auto& xi = p1.x;
	const auto& xj = p2.x;
	const double reciprocal = std::pow((xi - xj).euclidian_norm(), 3);
	const double scaling_factor = p1.m * p2.m / reciprocal;
	return scaling_factor * (xj - xi);
}
}  // namespace detail

/**
 * @brief Calculates the total gravitational force acting on each particle.
 *
 * Iterates over all particles and computes the net force on each particle
 * as the sum of pairwise interaction forces with all other particles.
 *
 * @param particles Mutable span over particles whose forces will be updated.
 */
constexpr void gravitational_force(std::span<Particle> particles) noexcept {
	for (std::size_t i = 0; i < particles.size(); ++i) {
		auto& p1 = particles[i];
		vec new_f{};
		for (std::size_t j = 0; j < particles.size(); ++j) {
			const auto& p2 = particles[j];
			if (i == j) continue;
			const auto f_ij = detail::calculate_component(p1, p2);
			new_f += f_ij;
		}
		p1.f = new_f;
	}
}
