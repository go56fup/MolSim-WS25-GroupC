#pragma once

#include <cmath>
#include <concepts>
#include <functional>
#include <stdexcept>

#include "grid/particle_container/fwd.hpp"
#include "grid/particle_container/system.hpp"
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
	requires(
		Candidate f, particle_container& container, particle_system::particle_id p1,
		particle_system::particle_id p2
	) {
		{ std::invoke(f, container, p1, p2) } -> std::same_as<void>;
	} &&
	std::is_nothrow_invocable_v<
		Candidate, particle_container&, particle_system::particle_id,
		particle_system::particle_id> &&
	std::is_trivially_copyable_v<Candidate>;

[[gnu::const]] constexpr double cubed(double x) noexcept {
	return x * x * x;
}

[[gnu::const]] constexpr double sixth_power(double x) {
	return x * x * x * x * x * x;
}

#if 0
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
#endif

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

struct lennard_jones_parameters {
	vec p1_position;
	vec p2_position;
	double sigma;
	double epsilon;
};

CONSTEXPR_IF_GCC inline vec lennard_jones_force(const lennard_jones_parameters& params) noexcept {
	const auto pos_diff = params.p1_position - params.p2_position;
	const double norm = pos_diff.euclidian_norm();
	assert(norm != 0 && "Two particles at the same position cannot interact.");
	const double normalized_sigma = sixth_power(params.sigma / norm);
	const double scaling_factor = (24 * params.epsilon / (norm * norm)) *
	                              (normalized_sigma - 2 * normalized_sigma * normalized_sigma);
	const auto result = -scaling_factor * pos_diff;
	TRACE_FORCES("Resulting force: {}", result);
	return result;
}

CONSTEXPR_IF_GCC inline void lennard_jones_force_soa(
	particle_container& container, particle_system::particle_id p1, particle_system::particle_id p2
) noexcept {
	TRACE_FORCES("Calculating Lennard-Jones forces for:\n{}\n{}\n", p1, p2);
	auto& system = container.system();
	const double pos_diff_x = system.x[p1] - system.x[p2];
	const double pos_diff_y = system.y[p1] - system.y[p2];
	const double pos_diff_z = system.z[p1] - system.z[p2];
	const double norm =
		std::sqrt(pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z);
	assert(norm != 0 && "Two particles at the same position cannot interact.");
	const auto& p1_material = container.material_for_particle(p1);
	const auto& p2_material = container.material_for_particle(p2);

	const double sigma = p1_material.sigma == p2_material.sigma
	                         ? p1_material.sigma
	                         : (p1_material.sigma + p2_material.sigma) / 2;
	const double eps = p1_material.epsilon == p2_material.epsilon
	                       ? p1_material.epsilon
	                       : std::sqrt(p1_material.epsilon * p2_material.epsilon);
	const double normalized_sigma = sixth_power(sigma / norm);
	const double scaling_factor =
		(24 * eps / (norm * norm)) * (normalized_sigma - 2 * normalized_sigma * normalized_sigma);

	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;

	system.fx[p1] += x_delta;
	system.fy[p1] += y_delta;
	system.fz[p1] += z_delta;

	system.fx[p2] -= x_delta;
	system.fy[p2] -= y_delta;
	system.fz[p2] -= z_delta;
}

constexpr void apply_gravity(particle_container& container, double gravity) noexcept {
	auto& system = container.system();
	for (particle_system::particle_id p = 0; p < system.size(); ++p) {
		system.fy[p] += container.material_for_particle(p).mass * gravity;
	}
}
