#pragma once

#include <cmath>

#include "grid/particle_container/fwd.hpp"
#include "grid/particle_container/system.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/tracing/macros.hpp"

[[gnu::const]] constexpr double cubed(double x) noexcept {
	return x * x * x;
}

[[gnu::const]] constexpr double sixth_power(double x) {
	return x * x * x * x * x * x;
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
 * @return A vector representing the force acted upon @p p1 by @p p2.

 **/

struct lennard_jones_parameters {
	vec p1_position;
	vec p2_position;
	double sigma{};
	double epsilon{};
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

// TODO(tuna): remove this
#if SINGLETHREADED
#define NO_ATOMICS
#endif
// #define NO_ATOMICS

CONSTEXPR_IF_GCC inline void lennard_jones_force_soa(
	particle_container& container, particle_system::particle_id p1, particle_system::particle_id p2,
	vec particle_mirror = {0, 0, 0}
) noexcept {
	auto& system = container.system();
	const double pos_diff_x = system.x[p1] - system.x[p2] + particle_mirror.x;
	const double pos_diff_y = system.y[p1] - system.y[p2] + particle_mirror.y;
	const double pos_diff_z = system.z[p1] - system.z[p2] + particle_mirror.z;
	const double r2 = pos_diff_x * pos_diff_x + pos_diff_y * pos_diff_y + pos_diff_z * pos_diff_z;
	assert(r2 != 0 && "Two particles at the same position cannot interact.");
	const double sigma = system.sigma[p1] == system.sigma[p2]
	                         ? system.sigma[p1]
	                         : (system.sigma[p1] + system.sigma[p2]) / 2;
	const double eps = system.epsilon[p1] == system.epsilon[p2]
	                       ? system.epsilon[p1]
	                       : std::sqrt(system.epsilon[p1] * system.epsilon[p2]);
	const double scaling_factor = get_scaling_factor(sigma, eps, r2);
	const auto x_delta = -scaling_factor * pos_diff_x;
	const auto y_delta = -scaling_factor * pos_diff_y;
	const auto z_delta = -scaling_factor * pos_diff_z;
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

CONSTEXPR_IF_GCC inline void lennard_jones_force_soa_batchwise(
	particle_container& container, particle_batch p1_batch, particle_batch p2_batch,
	vec particle_mirror = {0, 0, 0}
) noexcept {
	auto& system = container.system();

	batch<double> pos_diff_x;
	batch<double> pos_diff_y;
	batch<double> pos_diff_z;
	batch<double> scaling_factor;

	// Position differences
#ifndef NO_ATOMICS
#pragma omp simd
#endif
	for (std::size_t i = 0; i < batch_size; ++i) {
		const auto p1 = p1_batch[i];
		const auto p2 = p2_batch[i];

		pos_diff_x[i] = system.x[p1] - system.x[p2] + particle_mirror.x;
		pos_diff_y[i] = system.y[p1] - system.y[p2] + particle_mirror.y;
		pos_diff_z[i] = system.z[p1] - system.z[p2] + particle_mirror.z;
	}

	// Force scaling factor (matches scalar math exactly)
#ifndef NO_ATOMICS
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

		const double s1 = system.sigma[p1];
		const double s2 = system.sigma[p2];
		const double e1 = system.epsilon[p1];
		const double e2 = system.epsilon[p2];

		const double sigma = (s1 == s2) ? s1 : (s1 + s2) * 0.5;
		const double eps = (e1 == e2) ? e1 : std::sqrt(e1 * e2);
		scaling_factor[i] = get_scaling_factor(sigma, eps, r2);
	}

	batch<double> x_delta;
	batch<double> y_delta;
	batch<double> z_delta;

#ifndef NO_ATOMICS
#pragma omp simd
#endif
	for (std::size_t i = 0; i < batch_size; ++i) {
		x_delta[i] = -scaling_factor[i] * pos_diff_x[i];
		y_delta[i] = -scaling_factor[i] * pos_diff_y[i];
		z_delta[i] = -scaling_factor[i] * pos_diff_z[i];
	}

	for (std::size_t i = 0; i < batch_size; ++i) {
		const auto p1 = p1_batch[i];
		const auto p2 = p2_batch[i];
		TRACE_FORCES(
			"Force delta for {} - {}:  {} {} {}", p1, p2, x_delta[i], y_delta[i], z_delta[i]
		);

#ifndef NO_ATOMICS
#pragma omp atomic
#endif
		system.fx[p1] += x_delta[i];
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
		system.fy[p1] += y_delta[i];
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
		system.fz[p1] += z_delta[i];

#ifndef NO_ATOMICS
#pragma omp atomic
#endif
		system.fx[p2] -= x_delta[i];
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
		system.fy[p2] -= y_delta[i];
#ifndef NO_ATOMICS
#pragma omp atomic
#endif
		system.fz[p2] -= z_delta[i];
	}
}

#undef NO_ATOMICS

constexpr void apply_gravity(particle_container& container, double gravity) noexcept {
	auto& system = container.system();
#if !SINGLETHREADED
#pragma omp parallel for simd schedule(static)
#endif
	for (particle_system::particle_id p = 0; p < system.size(); ++p) {
		system.fy[p] += system.mass[p] * gravity;
	}
}
