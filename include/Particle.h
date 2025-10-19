/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>

#include "Vector.h"

class Particle {
private:
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
	vec f{};

	/**
	 * Force which was effective on this particle
	 */
	vec old_f{};

	/**
	 * Mass of this particle
	 */
	double m;

	/**
	 * Type of the particle. Use it for whatever you want (e.g. to separate
	 * molecules belonging to different bodies, matters, and so on)
	 */
	int type;

public:
	explicit Particle(int type = 0);

	Particle(const Particle& other);

	Particle(
		// for visualization, we need always 3 coordinates
	    // -> in case of 2d, we use only the first and the second
		vec x_arg, vec v_arg, double m_arg, int type = 0
	);

	virtual ~Particle();

	const vec& getX() const;

	void setX(const vec&) noexcept;

	const vec& getV() const;

	void setV(const vec&) noexcept;

	const vec& getF() const;

	void setF(const vec&) noexcept;

	const vec& getOldF() const;

	void setOldF(const vec&) noexcept;

	double getM() const;

	int getType() const;

	constexpr bool operator==(const Particle& other) const = default;

	std::string toString() const;
};

std::ostream& operator<<(std::ostream& stream, Particle& p);
