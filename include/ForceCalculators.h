#pragma once

#include <cmath>
#include <concepts>
#include <functional>

#include "CompilerTraits.h"
#include "Particle.h"

/**
 * @brief Concept that defines a force calculator.
 *
 * A callable satisfies this concept if it can be invoked as-if it shares this signature:
 * @verbatim vec_3d<double> operator()(Particle& p1, Particle& p2) @endverbatim
 * where @p p1 and @p p2 are the particles to calculate a given force for.
 *
 * @tparam Candidate Type being tested for the force calculator interface.
 *
 * @return `true` @a iff the callable can be invoked as specified.
 */
template <typename Candidate>
concept force_calculator = requires(Candidate f, const Particle& p1, const Particle& p2) {
	{ std::invoke(f, p1, p2) } -> std::same_as<vec>;
};

/**
 * @brief Computes the gravitational force exerted on one particle by another.
 *
 * Calculates the gravitational-like force acted upon particle @p p1 by particle @p p2
 * according to the inverse-cube law:
 * \f[
 * \vec{F}_{ij} = G \frac{m_i m_j}{r_{ij}^3} (\vec{x}_j - \vec{x}_i)
 * \f]
 * (where G is taken as 1 for simplicity).
 *
 * @param p1 The particle being acted upon (particle i).
 * @param p2 The particle exerting the force (particle j).
 * @return A vector representing the force acted upon @p p1 by @p p2.
 */
CONSTEXPR_IF_GCC inline vec gravitational_force(const Particle& p1, const Particle& p2) noexcept {
	const auto& xi = p1.x;
	const auto& xj = p2.x;
	const double reciprocal = std::pow((xi - xj).euclidian_norm(), 3);
	const double scaling_factor = p1.m * p2.m / reciprocal;
	return scaling_factor * (xj - xi);
}

/**
 * @brief Calculates the force resulting from the Lennard-Jones potential.
 *
 * The force acted upon particle `p1` by particle `p2` is calculated according to the
 * negative gradient of the Lennard-Jones potential:
 \f[
 *   \vec{F}_{ij}
 *   = -\frac{24 \varepsilon}{r_{ij}^2}
 *     \left[
 *       \left( \frac{\sigma}{r_{ij}} \right)^6
 *       - 2 \left( \frac{\sigma}{r_{ij}} \right)^{12}
 *     \right]
 *     (\vec{x}_i - \vec{x}_j)
 * \f]
 * (with \f$\varepsilon = 5\f$ and \f$\sigma = 12\f$ hardcoded).
 *
 * @param p1 The particle being acted upon (particle i).
 * @param p2 The particle exerting the force (particle j).
 * @return A vector representing the force acted upon @p p1 by @p p2.

 **/
CONSTEXPR_IF_GCC inline vec lennard_jones_force(const Particle& p1, const Particle& p2) noexcept {
	const auto& xi = p1.x;
	const auto& xj = p2.x;
	static constexpr int eps = 5;
	static constexpr int sig = 1;
	const double norm = (xi - xj).euclidian_norm();
	const double scaling_factor =
		24 * eps / std::pow(norm, 2) * (std::pow(sig / norm, 6) - 2 * (std::pow(sig / norm, 12)));
	return scaling_factor * (xj - xi);
}
