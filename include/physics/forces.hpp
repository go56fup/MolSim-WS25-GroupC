#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
#include "grid/particle_container/system.hpp"
#include "simulation/entities.hpp"
#include "simulation/particle.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/constants.hpp"
#include "utility/tracing/macros.hpp"

#define OVERRIDE_CALCULATE_F 1

[[gnu::const]] constexpr double cubed(double x) noexcept {
	return x * x * x;
}

[[gnu::const]] constexpr double sixth_power(double x) {
	return x * x * x * x * x * x;
}

struct lennard_jones_constants {
	double sigma{};
	double epsilon{};
};

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
 * @return A vector representing the force acted upon @p p1 by @p p2.

 **/



struct lennard_jones_parameters {
	vec p1_position;
	vec p2_position;
	lennard_jones_constants constants;
};

CONSTEXPR_IF_GCC inline vec lennard_jones_force(const lennard_jones_parameters& params) noexcept {
	const auto pos_diff = params.p1_position - params.p2_position;
	const double norm = pos_diff.euclidian_norm();
	assert(norm != 0 && "Two particles at the same position cannot interact.");
	const double normalized_sigma = sixth_power(params.constants.sigma / norm);
	const double scaling_factor = (24 * params.constants.epsilon / (norm * norm)) *
	                              (normalized_sigma - 2 * normalized_sigma * normalized_sigma);
	const auto result = -scaling_factor * pos_diff;
	TRACE_FORCES("Resulting force: {}", result);
	return result;
}

#define LENNARD_JONES_V2

constexpr double get_scaling_factor(double sigma, double eps, double r2) {
#ifdef LENNARD_JONES_V2
	const double inv_r2 = 1.0 / r2;
	// Using powers of r2 instead of sigma/norm avoids repeated divisions
	const double s_r_2 = (sigma * sigma) * inv_r2;  // (sigma/r)^2
	const double s_r_6 = s_r_2 * s_r_2 * s_r_2;     // (sigma/r)^6

	const double scaling_factor = (24.0 * eps * inv_r2) * (s_r_6 - 2.0 * s_r_6 * s_r_6);
#else
	const double norm = std::sqrt(r2);
	const double normalized_sigma = sixth_power(sigma / norm);
	const double scaling_factor =
		(24 * eps / r2) * (normalized_sigma - 2 * normalized_sigma * normalized_sigma);
#endif
	TRACE_FORCES("Scaling factor for sigma={}, eps={}, r2={}: {}", sigma, eps, r2, scaling_factor);
	return scaling_factor;
}

constexpr lennard_jones_constants
get_constants(const particle_system& system, particle_id p1, particle_id p2) noexcept {
	return {
		.sigma = system.sigma[p1] == system.sigma[p2] ? system.sigma[p1]
	                                                  : (system.sigma[p1] + system.sigma[p2]) * 0.5,
		.epsilon = system.epsilon[p1] == system.epsilon[p2]
	                   ? system.epsilon[p1]
	                   : std::sqrt(system.epsilon[p1] * system.epsilon[p2])
	};
}

#if SINGLETHREADED && !OVERRIDE_CALCULATE_F
#define NO_ATOMICS
#endif
// #define NO_ATOMICS

constexpr void apply_deltas_atomic(
	particle_system& system, particle_id p1, particle_id p2, double x_delta, double y_delta,
	double z_delta
) noexcept {
	TRACE_FORCES("Force delta for {} - {}:  {} {} {}", p1, p2, x_delta, y_delta, z_delta);

#ifndef NO_ATOMICS
#pragma omp atomic
#endif
	system.fx[p1] += x_delta;
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
	system.fy[p1] += y_delta;
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
	system.fz[p1] += z_delta;

#ifndef NO_ATOMICS
#pragma omp atomic
#endif
	system.fx[p2] -= x_delta;
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
	system.fy[p2] -= y_delta;
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
	system.fz[p2] -= z_delta;
}

CONSTEXPR_IF_GCC inline void lennard_jones_force_soa(
	particle_system& system, particle_id p1, particle_id p2, vec particle_mirror
) noexcept {
	const double pos_diff_x = system.x[p1] - system.x[p2] + particle_mirror.x;
	const double pos_diff_y = system.y[p1] - system.y[p2] + particle_mirror.y;
	const double pos_diff_z = system.z[p1] - system.z[p2] + particle_mirror.z;
	const double r2 = pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z;
	assert(r2 != 0 && "Two particles at the same position cannot interact.");

	const auto [sigma, eps] = get_constants(system, p1, p2);

	const double scaling_factor = get_scaling_factor(sigma, eps, r2);
	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;
	apply_deltas_atomic(system, p1, p2, x_delta, y_delta, z_delta);
}

CONSTEXPR_IF_GCC inline void lennard_jones_force_soa_batchwise(
	particle_system& system, particle_batch p1_batch, particle_batch p2_batch, vec particle_mirror
) noexcept {
	batch<double> pos_diff_x;
	batch<double> pos_diff_y;
	batch<double> pos_diff_z;
	batch<double> scaling_factor;

	// Position differences
#if !SINGLETHREADED
#pragma omp simd
#endif
	for (std::size_t i = 0; i < batch_size; ++i) {
		const auto p1 = p1_batch[i];
		const auto p2 = p2_batch[i];

		pos_diff_x[i] = system.x[p1] - system.x[p2] + particle_mirror.x;
		pos_diff_y[i] = system.y[p1] - system.y[p2] + particle_mirror.y;
		pos_diff_z[i] = system.z[p1] - system.z[p2] + particle_mirror.z;
	}

#if !SINGLETHREADED
#pragma omp simd
#endif
	for (std::size_t i = 0; i < batch_size; ++i) {
		const auto p1 = p1_batch[i];
		const auto p2 = p2_batch[i];

		const double dx = pos_diff_x[i];
		const double dy = pos_diff_y[i];
		const double dz = pos_diff_z[i];
		const double r2 = (dx * dx) + (dy * dy) + (dz * dz);

		assert(r2 != 0.0 && "Two particles at the same position cannot interact.");

		const auto [sigma, eps] = get_constants(system, p1, p2);
		scaling_factor[i] = get_scaling_factor(sigma, eps, r2);
	}

	batch<double> x_delta;
	batch<double> y_delta;
	batch<double> z_delta;

#if !SINGLETHREADED
#pragma omp simd
#endif
	for (std::size_t i = 0; i < batch_size; ++i) {
		x_delta[i] = -scaling_factor[i] * pos_diff_x[i];
		y_delta[i] = -scaling_factor[i] * pos_diff_y[i];
		z_delta[i] = -scaling_factor[i] * pos_diff_z[i];
	}

#if !SINGLETHREADED
#pragma omp simd
#endif
	for (std::size_t i = 0; i < batch_size; ++i) {
		const auto p1 = p1_batch[i];
		const auto p2 = p2_batch[i];
		apply_deltas_atomic(system, p1, p2, x_delta[i], y_delta[i], z_delta[i]);
	}
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
 */
CONSTEXPR_IF_GCC inline void smoothed_lennard_jones_force(
	particle_system& system, particle_id p1, particle_id p2, double cutoff_radius,
	double lower_radius, vec particle_mirror
) noexcept {
	const double pos_diff_x = system.x[p1] - system.x[p2] + particle_mirror.x;
	const double pos_diff_y = system.y[p1] - system.y[p2] + particle_mirror.y;
	const double pos_diff_z = system.z[p1] - system.z[p2] + particle_mirror.z;
	const double r2 = pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z;
	// TODO(tuna): check against cutoff * cutoff instead of sqrt
	const double norm = std::sqrt(r2);
	if (norm >= cutoff_radius) {
		return;
	}

	const auto [sigma, eps] = get_constants(system, p1, p2);

	double scaling_factor = 0;

	if (norm <= lower_radius) {
		scaling_factor = get_scaling_factor(sigma, eps, r2);
	} else {
		const double sigma_sixth_power = sixth_power(sigma);
		const double norm_sixth_power = sixth_power(norm);
		const double norm_seventh_power = norm * norm_sixth_power;
		const double norm_fourteen_power = norm_seventh_power * norm_seventh_power;
		scaling_factor =
			-24 * sigma_sixth_power * eps /
			(norm_fourteen_power * cubed(cutoff_radius - lower_radius)) * (cutoff_radius - norm) *
			((cutoff_radius * cutoff_radius) * (2 * sigma_sixth_power - norm_sixth_power) +
		     cutoff_radius * (3 * lower_radius - norm) *
		         (norm_sixth_power - 2 * sigma_sixth_power) +
		     norm * (5 * lower_radius * sigma_sixth_power - 2 * lower_radius * norm_sixth_power -
		             3 * sigma_sixth_power * norm + norm_seventh_power));
	}

	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;

#define INTERESTING_PARTICLE_LOG 0
#if INTERESTING_PARTICLE_LOG
	if (p1 == 80 || p2 == 80) {
		TRACE_FORCES(
			"interesting: {} <-> {}, delta={}, {}, {}, pos_delta={}, {}, {}; scaling_factor={}", p1,
			p2, x_delta, y_delta, z_delta, pos_diff_x, pos_diff_y, pos_diff_z, scaling_factor
		);
		TRACE_FORCES("p1: {}", system.representation(p1));
		TRACE_FORCES("p2: {}", system.representation(p2));
	}
#endif
	apply_deltas_atomic(system, p1, p2, x_delta, y_delta, z_delta);
}

/**
 * @brief Computes the harmonic force exerted on one particle by one of its horizontal or
 vertical neighbours in a membrane.
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
*/
CONSTEXPR_IF_GCC inline void harmonic_force_orthogonal(
	particle_system& system, particle_id p1, particle_id p2,
	const membrane_simulation_parameters& params
) {
	const double pos_diff_x = system.x[p1] - system.x[p2];
	const double pos_diff_y = system.y[p1] - system.y[p2];
	const double pos_diff_z = system.z[p1] - system.z[p2];
	const double norm =
		std::sqrt(pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z);

	const double scaling_factor = params.stiffness * (norm - params.average_bond_length) / norm;
	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;

	// No Newton's 3rd law, the membrane is walked in a linear and simple fashion to aid
	// parallelizability by avoiding atomics.
	system.fx[p1] += x_delta;
	system.fy[p1] += y_delta;
	system.fz[p1] += z_delta;
}

/**
 * @brief Computes the harmonic force exerted on one particle by one of its diagonal neighbours
 in a membrane.
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
 */
CONSTEXPR_IF_GCC inline void harmonic_force_diagonal(
	particle_system& system, particle_id p1, particle_id p2,
	const membrane_simulation_parameters& params
) {
	const double pos_diff_x = system.x[p1] - system.x[p2];
	const double pos_diff_y = system.y[p1] - system.y[p2];
	const double pos_diff_z = system.z[p1] - system.z[p2];
	const double norm =
		std::sqrt(pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z);

	const double scaling_factor =
		params.stiffness * (norm - std::numbers::sqrt2 * params.average_bond_length) / norm;
	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;

	system.fx[p1] += x_delta;
	system.fy[p1] += y_delta;
	system.fz[p1] += z_delta;
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
 */
CONSTEXPR_IF_GCC inline void
gravitational_force_soa(particle_system& system, particle_id p1, particle_id p2) noexcept {
	const double pos_diff_x = system.x[p1] - system.x[p2];
	const double pos_diff_y = system.y[p1] - system.y[p2];
	const double pos_diff_z = system.z[p1] - system.z[p2];
	const double norm =
		std::sqrt(pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z);

	const double reciprocal = cubed(norm);
	const double scaling_factor = system.mass[p1] * system.mass[p2] / reciprocal;

	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;
	apply_deltas_atomic(system, p1, p2, x_delta, y_delta, z_delta);
}

template <axis a>
constexpr void apply_gravity(particle_container& container, double gravity) noexcept {
	auto& system = container.system();
#if !SINGLETHREADED
#pragma omp parallel for simd schedule(static)
#endif
	for (particle_id p = 0; p < system.size(); ++p) {
		system.force_component(a)[p] += system.mass[p] * gravity;
	}
}

CONSTEXPR_IF_GCC inline void truncated_lennard_jones_force_soa(
	particle_system& system, particle_id p1, particle_id p2, vec particle_mirror
) noexcept {
	const double pos_diff_x = system.x[p1] - system.x[p2] + particle_mirror.x;
	const double pos_diff_y = system.y[p1] - system.y[p2] + particle_mirror.y;
	const double pos_diff_z = system.z[p1] - system.z[p2] + particle_mirror.z;
	const double r2 = pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z;
	const double norm = std::sqrt(r2);

	const auto [sigma, eps] = get_constants(system, p1, p2);

	if (norm >= sigma * sixth_root_of_2) return;

	const double scaling_factor = get_scaling_factor(sigma, eps, r2);

	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;
	apply_deltas_atomic(system, p1, p2, x_delta, y_delta, z_delta);
}

CONSTEXPR_IF_GCC inline vec truncated_lennard_jones_force(const lennard_jones_parameters& params
) noexcept {
	const auto pos_diff = params.p1_position - params.p2_position;
	const double norm = pos_diff.euclidian_norm();
	assert(norm != 0 && "Two particles at the same position cannot interact.");
	if (norm >= params.constants.sigma * sixth_root_of_2) return {};

	const double normalized_sigma = sixth_power(params.constants.sigma / norm);
	const double scaling_factor = (24 * params.constants.epsilon / (norm * norm)) *
	                              (normalized_sigma - 2 * normalized_sigma * normalized_sigma);
	return -scaling_factor * pos_diff;
}

CONSTEXPR_IF_GCC inline void smoothed_lennard_jones_force_pseudo_batched(
	particle_system& system, particle_batch p1_batch, particle_batch p2_batch, double cutoff_radius,
	double lower_radius, vec particle_mirror
) {
	for (std::size_t i = 0; i < batch_size; ++i) {
		smoothed_lennard_jones_force(
			system, p1_batch[i], p2_batch[i], cutoff_radius, lower_radius, particle_mirror
		);
	}
}

CONSTEXPR_IF_GCC inline void truncated_lennard_jones_force_pseudo_batched(
	particle_system& system, particle_batch p1_batch, particle_batch p2_batch, vec particle_mirror
) {
	for (std::size_t i = 0; i < batch_size; ++i) {
		truncated_lennard_jones_force_soa(system, p1_batch[i], p2_batch[i], particle_mirror);
	}
}

CONSTEXPR_IF_GCC inline void gravitational_force_soa_pseudo_batched(
	particle_system& system, particle_batch p1_batch, particle_batch p2_batch
) {
	for (std::size_t i = 0; i < batch_size; ++i) {
		gravitational_force_soa(system, p1_batch[i], p2_batch[i]);
	}
}

constexpr vec
force_calculator::single(const sim_configuration& config, const lennard_jones_parameters& params) {
	const auto& method = config.force_method;
	switch (method) {
	case lennard_jones:
		return lennard_jones_force(params);
	case truncated_lennard_jones:
		return truncated_lennard_jones_force(params);
	default:
		throw std::logic_error(
			fmt::format("Unimplemented single force calculation for method: {}", method)
		);
	}
}

constexpr void force_calculator::soa(
	particle_system& system, const sim_configuration& config, particle_id p1, particle_id p2,
	vec particle_mirror = {0, 0, 0}
) {
	switch (config.force_method) {
	case lennard_jones:
		lennard_jones_force_soa(system, p1, p2, particle_mirror);
		return;
	case smoothed_lennard_jones:
		smoothed_lennard_jones_force(
			system, p1, p2, config.cutoff_radius, config.lower_radius, particle_mirror
		);
		return;
	case truncated_lennard_jones:
		truncated_lennard_jones_force_soa(system, p1, p2, particle_mirror);
		return;
	case gravity:
		gravitational_force_soa(system, p1, p2);
		return;
	}
	std::unreachable();
}

constexpr void force_calculator::batch(
	particle_system& system, const sim_configuration& config, particle_batch p1, particle_batch p2,
	vec particle_mirror = {0, 0, 0}
) {
	switch (config.force_method) {
	case lennard_jones:
		lennard_jones_force_soa_batchwise(system, p1, p2, particle_mirror);
		return;
	case smoothed_lennard_jones:
		smoothed_lennard_jones_force_pseudo_batched(
			system, p1, p2, config.cutoff_radius, config.lower_radius, particle_mirror
		);
		return;
	case truncated_lennard_jones:
		truncated_lennard_jones_force_pseudo_batched(system, p1, p2, particle_mirror);
		return;
	case gravity:
		gravitational_force_soa_pseudo_batched(system, p1, p2);
		return;
	}
	std::unreachable();
}
