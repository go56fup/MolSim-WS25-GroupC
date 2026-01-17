#pragma once

#include <cmath>
#include <concepts>
#include <functional>
#include <stdexcept>
#include <numbers>

#include "physics/particle.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/tracing/macros.hpp"

/**
 * @brief Concept that defines a force calculator.
 *
 * A callable satisfies this concept if it can be invoked as-if it shares this signature:
 * @verbatim vec_3d<double> operator()(particle& p1, particle& p2) @endverbatim
 * where @p p1 and @p p2 are the particles to calculate a given force for.
 *
 * @tparam Candidate Type being tested for the force calculator interface.
 *
 * @return `true` @a iff the callable can be invoked as specified.
 */
template <typename Candidate>
concept force_calculator =
	requires(Candidate f, const particle& p1, const particle& p2) {
		{ std::invoke(f, p1, p2) } -> std::same_as<vec>;
	} && std::is_nothrow_invocable_v<Candidate, const particle&, const particle&> &&
	std::is_trivially_copyable_v<Candidate>;

[[gnu::const]] constexpr double cubed(double x) noexcept {
	return x * x * x;
}

[[gnu::const]] constexpr double sixth_power(double x) {
	return x * x * x * x * x * x;
}


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
CONSTEXPR_IF_GCC inline vec gravitational_force(const particle& p1, const particle& p2) noexcept {
	const auto& xi = p1.x;
	const auto& xj = p2.x;
	const double reciprocal = cubed((xi - xj).euclidian_norm());
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
CONSTEXPR_IF_GCC inline vec lennard_jones_force(const particle& p1, const particle& p2) noexcept {
	TRACE_FORCES("Calculating Lennard-Jones forces for:\n{}\n{}\n", p1, p2);
	const auto xi = p1.x;
	const auto& xj = p2.x;
	const double norm = (xi - xj).euclidian_norm();
	assert(norm != 0 && "Two particles at the same position cannot interact.");
	const double sigma = (p1.sigma + p2.sigma) / 2;
	const double eps = std::sqrt(p1.epsilon * p2.epsilon);
	const double normalized_sigma = sixth_power(sigma / norm);
	const double scaling_factor =
		(24 * eps / (norm * norm)) * (normalized_sigma - 2 * normalized_sigma * normalized_sigma);
	const auto result = scaling_factor * (xj - xi);
	TRACE_FORCES("Resulting force: {}", result);
	return result;
}


/**
 * @brief Calculates the force resulting from the smoothed Lennard-Jones potential.
 *
 * The force acted upon particle `p1` by particle `p2` is calculated according to the
 * negative gradient of the smoothed Lennard-Jones potential:
 *
 *
 * @param p1 The particle being acted upon (particle i).
 * @param p2 The particle exerting the force (particle j).
 * @param cutoff_radius Cutoff radius in the simulation.
 * @param lower_radius Radius resulting in the regular lennard jones force calculation.
 * @return A vector representing the force acted upon @p p1 by @p p2.
 */
CONSTEXPR_IF_GCC inline vec smoothed_lennard_jones_force(const particle& p1, const particle& p2, const double cutoff_radius, const double lower_radius) noexcept {
	const auto& xi = p1.x;
	const auto& xj = p2.x;
	const double norm = (xi - xj).euclidian_norm();
	if (norm >= cutoff_radius) {
		return {0,0,0};
	}
	if (norm <= lower_radius) {
		return lennard_jones_force(p1, p2);
	}
	const double eps = std::sqrt(p1.epsilon * p2.epsilon);
	const double sigma_sixth_power = sixth_power((p1.sigma + p2.sigma) / 2);
	const double norm_sixth_power = sixth_power(norm);
	const double norm_seventh_power = norm * norm_sixth_power;
	const double norm_fourteen_power = norm_seventh_power * norm_seventh_power;
	const double scaling_factor = -24 * sigma_sixth_power * eps/ (norm_fourteen_power * cubed(cutoff_radius - lower_radius)) * (cutoff_radius - norm)
	* ((cutoff_radius * cutoff_radius) * (2 * sigma_sixth_power - norm_sixth_power)
	+ cutoff_radius *  (3 * lower_radius - norm) * (norm_sixth_power - 2 * sigma_sixth_power)
	+ norm * (5 * lower_radius * sigma_sixth_power - 2 * lower_radius * norm_sixth_power - 3* sigma_sixth_power * norm + norm_seventh_power));
	const auto result = scaling_factor * (xj - xi);
	return result;
}

/**
 * @brief Computes the harmonic force exerted on one particle by one of its horizontal or vertical neighbours in a membrane.
 *
\f[
 *   \vec{F}_{ij}
 *   = k
 *       \left(r_{ij} - r_0 \right)
 *     \frac{(\vec{x}_j - \vec{x}_i)}{r_{ij}}
 * \f]
 *
 * @param p1 The particle being acted upon (particle i).
 * @param p2 The particle exerting the force (particle j).
 * @param average_bond_length The average bond length of a particle pair.
 * @param k Stiffness constant.
 * @return A vector representing the force acted upon @p p1 by @p p2.
 */
CONSTEXPR_IF_GCC inline vec harmonic_force_orthogonal(const particle& p1, const particle& p2, const double average_bond_length, const double k) {
	const auto& xi = p1.x;
	const auto& xj = p2.x;
	const double norm = (xi - xj).euclidian_norm();
	const double scaling_factor = k * (norm - average_bond_length) / norm;
	const auto result = scaling_factor * (xj - xi);
	return result;
}

/**
 * @brief Computes the harmonic force exerted on one particle by one of its diagonal neighbours in a membrane.
 *
 \f[
 *   \vec{F}_{ij}
 *   = k
 *       \left(r_{ij} - \sqrt{2} r_0\right)
 *     \frac{(\vec{x}_j - \vec{x}_i)}{r_{ij}}
 * \f]
 *
 * @param p1 The particle being acted upon (particle i).
 * @param p2 The particle exerting the force (particle j).
 * @param average_bond_length The average bond length of a particle pair.
 * @param k Stiffness constant.
 * @return A vector representing the force acted upon @p p1 by @p p2.
 */
CONSTEXPR_IF_GCC inline vec harmonic_force_diagonal(const particle& p1, const particle& p2, const double average_bond_length, const double k) {
	const auto& xi = p1.x;
	const auto& xj = p2.x;
	const double norm = (xi - xj).euclidian_norm();
	const double scaling_factor = k * (norm  - std::numbers::sqrt2 * average_bond_length) / norm;
	const auto result = scaling_factor * (xj - xi);
	return result;
}

/**
 * @brief Computes an external force in the z-direction acting on certain particles in a membrane.
 *
 * @param particles Particles on which the force is applied.
 * @param upward_force The magnitude of the force.
 */
constexpr void apply_upward_force(range_of<particle> auto&& particles, const double upward_force) noexcept {
	for (particle& p : particles) {
		p.f += {0, 0, upward_force};
	}
}


/**
 * @brief Computes an external gravitational force in the y-direction acting on particles.
 *
 * @param particles Particles on which gravity is applied.
 * @param gravity Gravitational constant.
 */
constexpr void apply_gravity(range_of<particle> auto&& particles, double gravity) noexcept {
	for (particle& p : particles) {
		p.f += {0, p.m * gravity, 0};
	}
}
