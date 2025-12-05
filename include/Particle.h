/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <iostream>
#include <source_location>
#include <string>
#include <string_view>

#include <fmt/base.h>

#include "Debug.h"
#include "Vector.h"

/**
 * @brief Represents a single particle with position, velocity, force, and mass.
 */
// Discarding the argument of the special member function in flag_special_member_funcs is
// on purpose.
// NOLINTNEXTLINE(*slicing)
class Particle
/// @cond DO_NOT_DOCUMENT
// TODO(tuna): without this guard GCC cannot destruct arrays of Particles at constexpr time,
// which I suspect is a bug, but I'll look into it (clang manages fine)
#ifndef IS_GCC
LOG_SPECIAL_MEMBER_FUNCS_DEBUG("Particle")
#endif
/// @endcond
{
public:
	/**
	 * Position of the particle
	 */
	vec x;

	/**
	 * Velocity of the particle
	 */
	vec v;

	/**
	 * Force effective on this particle
	 */
	vec f;

	/**
	 * Force which was effective on this particle
	 */
	vec old_f;

	/**
	 * Mass of this particle
	 */
	double m;

	/**
	 * Type of the particle. Use it for whatever you want (e.g. to separate
	 * molecules belonging to different bodies, matters, and so on)
	 */
	int type;

	/**
	 * @brief Default constructor.
	 *
	 * Creates a particle of the given type with zero-initialized attributes.
	 *
	 * @param type_arg The particle type identifier (defaults to 0).
	 */
	constexpr explicit Particle(int type_arg = 0) noexcept
		: x(0, 0, 0)
		, v(0, 0, 0)
		, m(0.)
		, type(type_arg) {}

	/**
	 * @brief Parameterized constructor.
	 *
	 * Initializes a particle with the given position, velocity, mass, and type.
	 *
	 * @param x_arg Initial position.
	 * @param v_arg Initial velocity.
	 * @param m_arg Mass of particle.
	 * @param type_arg Type identifier of particle (defaults to 0).
	 * @param loc Captured source location for logging (do not set manually).
	 */
	constexpr Particle(
		// for visualization, we need always 3 coordinates
	    // -> in case of 2d, we use only the first and the second

		// I would rather not change the given implementation here.
	    // NOLINTNEXTLINE(*easily-swappable-parameters)
		vec x_arg, vec v_arg, double m_arg, int type_arg = 0,
		const std::source_location& loc = std::source_location::current()
	) noexcept
		: x(MOVE_IF_DEBUG(x_arg))
		, v(MOVE_IF_DEBUG(v_arg))
		, m(m_arg)
		, type(type_arg)

	{
		flag_special_member_funcs<"Particle">::log_special_mem("parameterized construction", loc);
	}

	/**
	 * @brief Compares two particles for equality.
	 * @param other The particle to compare with.
	 * @return `true` @a iff all properties match.
	 */
	constexpr bool operator==(const Particle& other) const = default;
};

/**
 * @brief Stream output operator for Particle.
 *
 * Allows output of particles using `operator<<`.
 *
 * @param stream Output stream to print into.
 * @param p Particle to print.
 * @return Reference to the modified output stream.
 */
constexpr std::ostream& operator<<(std::ostream& stream, Particle& p) {
	return stream << fmt::format("Particle X: {} v: {} f: {} old_f: {} type: {}", p.x, p.v, p.f, p.old_f, p.type);
}

#ifdef NDEBUG
static_assert(std::is_trivially_copyable_v<Particle>);
#endif
