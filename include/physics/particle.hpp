/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <concepts>
#include <iostream>
#include <source_location>
#include <string>
#include <string_view>
#include <utility>

#include <fmt/base.h>

#include "physics/vec_3d.hpp"
#include "utility/tracing/special_memfun.hpp"

/**
 * @brief Represents a single particle with position, velocity, force, and mass.
 */
// Discarding the argument of the special member function in flag_special_member_funcs is
// on purpose.
// NOLINTNEXTLINE(*slicing)
class particle
/// @cond DO_NOT_DOCUMENT
TRACE_SPECIAL_MEMFUNS_FOR_CLASS("particle")
/// @endcond
{
private:
	using flag = trace_special_memfuns<"particle">;
	using vec_args = std::tuple<double&, double&, double&>;

public:
	/// Position of the particle
	vec x{0, 0, 0};

	/// Velocity of the particle
	vec v{0, 0, 0};

	/// Force effective on this particle
	vec f{0, 0, 0};

	/// Force which was effective on this particle
	vec old_f{0, 0, 0};

	/// Mass of this particle
	double m{};

	double sigma{};
	double epsilon{};

	/**
	 * @brief Type of the particle.
	 * Used to separate molecules belonging to different bodies, matters; identify individual
	 * particles, etc.
	 */
	int type = 0;

	// TODO(tuna): fix docs
	/**
	 * @brief Default constructor for a particle.
	 *
	 * Creates a particle of the given type with zero-initialized attributes.
	 *
	 */
	constexpr particle() noexcept = default;

	/**
	 * @brief Construct a default particle with specific type.
	 *
	 * @param type_arg The particle type identifier.
	 */
	constexpr explicit particle(int type_arg) noexcept
		: type(type_arg) {}

	/**
	 * @brief Construct a particle with specific parameters.
	 *
	 * Initializes a particle with the given position, velocity, mass, and type.
	 *
	 * @param x_arg Initial position.
	 * @param v_arg Initial velocity.
	 * @param m_arg Mass of particle.
	 * @param type_arg Type identifier of particle (defaults to 0).
	 */

	// I would rather not change the given implementation here.
	// NOLINTBEGIN(*easily-swappable-parameters)
	constexpr particle(
		vec x_arg, vec v_arg, double m_arg, double sigma_, double epsilon_, int type_arg = 0
	) noexcept
		: x(MOVE_IF_DEBUG(x_arg))
		, v(MOVE_IF_DEBUG(v_arg))
		, m(m_arg)
		, sigma(sigma_)
		, epsilon(epsilon_)
		, type(type_arg)

	{
		flag::annotate_construction("parameterized");
	}

	constexpr particle(
		vec x_arg, vec v_arg, vec f_arg, vec old_f_arg, double m_arg, double sigma_,
		double epsilon_, int type_arg
	) noexcept
		: x(MOVE_IF_DEBUG(x_arg))
		, v(MOVE_IF_DEBUG(v_arg))
		, f(MOVE_IF_DEBUG(f_arg))
		, old_f(MOVE_IF_DEBUG(old_f_arg))
		, m(m_arg)
		, sigma(sigma_)
		, epsilon(epsilon_)
		, type(type_arg)

	{
		flag::annotate_construction("parameterized");
	}

	// TODO(tuna): see if this is even used
	/**
	 * @brief Piecewise construct a particle.
	 *
	 * @tparam VecX Deduced forwarding reference type for `x_arg`.
	 * @tparam VecV Deduced forwarding reference type for `v_arg`.
	 *
	 * @param x_arg Initial position.
	 * @param v_arg Initial velocity.
	 * @param m_arg Mass of particle.
	 * @param type_arg Type identifier of particle (defaults to 0).
	 */
	template <std::same_as<vec_args> VecX, std::same_as<vec_args> VecV>
	constexpr particle(
		// See previous constructor for rationale.
	    // NOLINTNEXTLINE(*easily-swappable-parameters)
		std::piecewise_construct_t, VecX&& x_arg, VecV&& v_arg, double m_arg, int type_arg = 0

	) noexcept
		: x(std::make_from_tuple<vec>(std::forward<VecX>(x_arg)))
		, v(std::make_from_tuple<vec>(std::forward<VecV>(v_arg)))
		, m(m_arg)
		, type(type_arg) {
		flag::annotate_construction("piecewise");
	}

	// NOLINTBEGIN(*easily-swappable-parameters)

	/**
	 * @brief Compares two particles for equality.
	 * @param other The particle to compare with.
	 * @return `true` @a iff all properties match.
	 */
	constexpr bool operator==(const particle& other) const = default;
};

/**
 * @brief Specialization of the {fmt} API for `particle`.
 *
 * Enables use of fmt::format and fmt::print with `particle` objects.
 *
 */
// NOLINTBEGIN(*convert-member-functions-to-static)
template <>
struct fmt::formatter<particle> {
	/** @brief Parses the format specification (no-op for this type). */
	constexpr auto parse(fmt::format_parse_context& ctx) {
		return ctx.begin();
	}

	/** @brief Formats `particle` objects. */
	auto format(const particle& p, fmt::format_context& ctx) const {
		return fmt::format_to(
			ctx.out(), "Particle X: {} v: {} f: {} old_f: {} epsilon: {} sigma: {} type: {}", p.x,
			p.v, p.f, p.old_f, p.epsilon, p.sigma, p.type
		);
	}
};

// NOLINTEND(*convert-member-functions-to-static)

#ifdef NDEBUG
static_assert(std::is_trivially_copyable_v<particle>);
#endif
