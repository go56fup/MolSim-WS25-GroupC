#pragma once

#include <cmath>
#include <concepts>
#include <functional>
#include <omp.h>

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
	TRACE_FORCES(
		"Calculating Lennard-Jones forces between {}, {} on thread {}", p1, p2, omp_get_thread_num()
	);
	auto& system = container.system();
	const double pos_diff_x = system.x[p1] - system.x[p2];
	const double pos_diff_y = system.y[p1] - system.y[p2];
	const double pos_diff_z = system.z[p1] - system.z[p2];
	const double norm =
		std::sqrt(pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z);
	assert(norm != 0 && "Two particles at the same position cannot interact.");
	const double sigma = system.sigma[p1] == system.sigma[p2]
	                         ? system.sigma[p1]
	                         : (system.sigma[p1] + system.sigma[p2]) / 2;
	const double eps = system.epsilon[p1] == system.epsilon[p2]
	                       ? system.epsilon[p1]
	                       : std::sqrt(system.epsilon[p1] * system.epsilon[p2]);
	const double normalized_sigma = sixth_power(sigma / norm);
	const double scaling_factor =
		(24 * eps / (norm * norm)) * (normalized_sigma - 2 * normalized_sigma * normalized_sigma);

	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;

#pragma omp atomic
	system.fx[p1] += x_delta;
#pragma omp atomic
	system.fy[p1] += y_delta;
#pragma omp atomic
	system.fz[p1] += z_delta;

#pragma omp atomic
	system.fx[p2] -= x_delta;
#pragma omp atomic
	system.fy[p2] -= y_delta;
#pragma omp atomic
	system.fz[p2] -= z_delta;
}

inline constexpr std::size_t BATCH_SIZE = 8;

CONSTEXPR_IF_GCC inline void lennard_jones_force_soa_batchwise(
	particle_container& container, std::array<particle_system::particle_id, BATCH_SIZE> p1_batch,
	std::array<particle_system::particle_id, BATCH_SIZE> p2_batch
) noexcept {

	auto& system = container.system();

	double pos_diff_x_batch[BATCH_SIZE];
	double pos_diff_y_batch[BATCH_SIZE];
	double pos_diff_z_batch[BATCH_SIZE];

	double scaling_factor_batch[BATCH_SIZE];

#pragma omp simd
	for (std::size_t i = 0; i < BATCH_SIZE; i++) {
		auto p1 = p1_batch[i];
		auto p2 = p2_batch[i];

		pos_diff_x_batch[i] = system.x[p1] - system.x[p2];
		pos_diff_y_batch[i] = system.y[p1] - system.y[p2];
		pos_diff_z_batch[i] = system.z[p1] - system.z[p2];
	}

	for (std::size_t i = 0; i < BATCH_SIZE; i++) {
		auto p1 = p1_batch[i];
		auto p2 = p2_batch[i];
		TRACE_FORCES("Calculating Lennard-Jones forces for:\n{}\n{}\n", p1, p2);
		const double pos_diff_x = pos_diff_x_batch[i];
		const double pos_diff_y = pos_diff_y_batch[i];
		const double pos_diff_z = pos_diff_z_batch[i];
		const double norm =
			std::sqrt(pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z);
		assert(norm != 0 && "Two particles at the same position cannot interact.");

		const double sigma = system.sigma[p1] == system.sigma[p2]
		                         ? system.sigma[p1]
		                         : (system.sigma[p1] + system.sigma[p2]) / 2;
		const double eps = system.epsilon[p1] == system.epsilon[p2]
		                       ? system.epsilon[p1]
		                       : std::sqrt(system.epsilon[p1] * system.epsilon[p2]);

		const double normalized_sigma = sixth_power(sigma / norm);
		scaling_factor_batch[i] = (24 * eps / (norm * norm)) *
		                          (normalized_sigma - 2 * normalized_sigma * normalized_sigma);
	}

	double x_delta[BATCH_SIZE];
	double y_delta[BATCH_SIZE];
	double z_delta[BATCH_SIZE];

#pragma omp simd
	for (std::size_t i = 0; i < BATCH_SIZE; i++) {
		x_delta[i] = -scaling_factor_batch[i] * pos_diff_x_batch[i];
		y_delta[i] = -scaling_factor_batch[i] * pos_diff_y_batch[i];
		z_delta[i] = -scaling_factor_batch[i] * pos_diff_z_batch[i];
	}

	for (std::size_t i = 0; i < BATCH_SIZE; i++) {
		[[maybe_unused]] auto p1 = p1_batch[i];
		[[maybe_unused]] auto p2 = p2_batch[i];

		system.fx[p1] += x_delta[i];
		system.fy[p1] += y_delta[i];
		system.fz[p1] += z_delta[i];

		// NO 3rd law of Newton opt.
	}
}

constexpr void apply_gravity(particle_container& container, double gravity) noexcept {
	auto& system = container.system();
#pragma omp parallel for simd schedule(static)
	for (particle_system::particle_id p = 0; p < system.size(); ++p) {
		system.fy[p] += system.mass[p] * gravity;
	}
}
