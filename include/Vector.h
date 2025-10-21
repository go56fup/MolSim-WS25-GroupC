#pragma once

#include <array>
#include <cmath>
#include <format>
#include <span>

namespace detail {
/**
 * Constrain type to floating point or integral
 */
template <typename T>
concept numeric = std::floating_point<T> || std::integral<T>;
}

/**
 * Three-dimensional vector having a floating point or integral type and operations on the vector
 */
template <detail::numeric Value>
struct vec_3d {
	/**
	 * Three components of the vector with the associated type
	 */
	Value x;
	Value y;
	Value z;

	/**
	 * Add another vector to this vector
	 * @param other
	 * @return Component wise addition to this vector
	 */
	constexpr auto operator+=(const vec_3d<Value>& other) noexcept {
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	/**
	 * Subtract another vector from this vector
	 * @param other
	 * @return Component wise substraction from this vector
	 */
	constexpr auto operator-=(const vec_3d<Value>& other) noexcept {
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

	/**
	 * Scale this vector
	 * @param scalar
	 * @return Component wise scaling of this vector
	 */
	constexpr auto operator*=(detail::numeric auto scalar) noexcept {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	/**
	 * Calculate the norm of this vector
	 * @return Square root of the summed squares of the components
	 */
	double euclidian_norm() const noexcept {
		return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
	}

	/**
	 * @returns A pointer to the first component of this vector
	 */
	constexpr auto begin() noexcept {
		return &x;
	}
	/**
	 * @returns A pointer to the first component of this vector for constant values
	 */
	constexpr auto begin() const noexcept {
		return &x;
	}
	/**
	 * @returns A pointer to the last component of this vector
	 */
	constexpr auto end() noexcept {
		return &x + 3;
	}
	/**
	 * @returns A pointer to the last component of this vector for constant values
	 */
	constexpr auto end() const noexcept {
		return &x + 3;
	}

	/**
	 * @returns A pointer the data of this vector for constant values
	 */
	constexpr const Value* data() const noexcept {
		return &x;
	}

	/**
	 * @returns A pointer to the data of this vector
	 */
	constexpr Value* data() noexcept {
		return &x;
	}
};


/**
 * Addition of two vectors
 * @tparam Value
 * @param lhs
 * @param rhs
 * @return Added components of two vectors
 */
template <detail::numeric Value>
constexpr auto operator+(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) noexcept {
	return vec_3d{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}
/**
 * Subtraction of two vectors
 * @tparam Value
 * @param lhs
 * @param rhs
 * @return Subtracted components of two vectors
 */
template <detail::numeric Value>
constexpr auto operator-(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) noexcept {
	return vec_3d{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

/**
 * Scalar multiplication of a vector
 * @tparam Value
 * @param scalar
 * @param vector
 * @return Scaled components of the vector
 */
template <detail::numeric Value>
constexpr auto operator*(detail::numeric auto scalar, const vec_3d<Value>& vector) noexcept {
	return vec_3d{vector.x * scalar, vector.y * scalar, vector.z * scalar};
}

/**
 * Scalar multiplication of a vector
 * @tparam Value
 * @param vector
 * @param scalar
 * @return Scaled components of the vector
 */
template <detail::numeric Value>
constexpr auto operator*(const vec_3d<Value>& vector, detail::numeric auto scalar) noexcept {
	return vec_3d{vector.x * scalar, vector.y * scalar, vector.z * scalar};
}

/**
 *
 * @tparam Value
 * @param stream
 * @param vector
 * @return Output stream
 */
template <detail::numeric Value>
std::ostream& operator<<(std::ostream& stream, const vec_3d<Value>& vector) {
	return stream << std::format("({}, {}, {})", vector.x, vector.y, vector.z);
}

/**
 * Component wise comparison of two vectors
 * @tparam Value
 * @param lhs
 * @param rhs
 * @return True iff all components of the vectors are the same
 */
template <detail::numeric Value>
constexpr bool operator==(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) {
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

/**
 *  Alias for a 3D vector of doubles.
 * @typedef vec
 */
using vec = vec_3d<double>;
static_assert(sizeof(vec) == sizeof(double[3]));
