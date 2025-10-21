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

	/**
	 * Destructor for this particle
	 */
	virtual ~Particle();

	/**
	 * Returns the position of this particle
	 */
	const vec& getX() const;

	/**
	 * Sets the position of this particle
	 */
   void setX(const vec&) noexcept;

	/**
	 * Returns the velocity of this particle
	 */
	const vec& getV() const;

   /**
	* Sets the velocity of this particle
	*/
	void setV(const vec&) noexcept;

	/**
	 * Returns the force effective on this particle
	 */
	const vec& getF() const;

	/**
	 * Sets the force effective on this particle
	 */
	void setF(const vec&) noexcept;

	/**
	 * Returns the previous force effective on this particle
	 */
	const vec& getOldF() const;

	/**
	 * Sets the previous force effective on this particle
	 */
	void setOldF(const vec&) noexcept;

	/**
	 * Returns the mass of this particle
	 */
	double getM() const;


	/**
	 * Returns the type of this particle
	 */
	int getType() const;

	/**
	 * @param other
	 * @return True iff other particle is equivalent to this particle
	 */
   constexpr bool operator==(const Particle& other) const = default;

	/**
	 * Returns the string representation of this particle
	 */
	std::string toString() const;
};

/**
 * @param stream
 * @param p
 * @return Output stream of the particle
 */
std::ostream& operator<<(std::ostream& stream, Particle& p);
