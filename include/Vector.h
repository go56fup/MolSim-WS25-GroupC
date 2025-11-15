#pragma once

#include <cmath>
#include <ostream>

#include <fmt/base.h>
#include <fmt/format.h>

#include "CompilerTraits.h"

namespace detail {
/**
 * @brief Concept constraining a type to be arithmetic.
 *
 * This concept ensures that the given type `T` is either an integral
 * or a floating-point type (i.e., satisfies `std::is_arithmetic_v<T>`).
 *
 * @tparam T The type to check.
 */
template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;
}  // namespace detail

/**
 * @brief Three-dimensional vector with arithmetic components.
 *
 * It supports basic arithmetic operations (addition, subtraction, and scaling),
 * as well as iteration and norm calculation.
 *
 * @tparam Value The numeric type of the vector components. Must satisfy detail::arithmetic.
 */
template <detail::arithmetic Value>
// Discarding the argument of the special member function in flag_special_member_funcs is
// on purpose.
// NOLINTNEXTLINE(*slicing)
struct vec_3d
#ifndef IS_GCC
/// @cond DO_NOT_DOCUMENT
LOG_SPECIAL_MEMBER_FUNCS_DEBUG("vec_3d")
/// @endcond
#endif
{
public:
	/** @brief X component of the vector. */
	Value x;

	/** @brief Y component of the vector. */
	Value y;

	/** @brief Z component of the vector. */
	Value z;

	constexpr vec_3d(Value arg_x, Value arg_y, Value arg_z) noexcept
		: x(arg_x)
		, y(arg_y)
		, z(arg_z) {}

	constexpr vec_3d() noexcept
		: x{}
		, y{}
		, z{} {}

	/**
	 * @brief Adds another vector to this vector.
	 *
	 * Performs component-wise addition.
	 *
	 * @param other The vector to add.
	 * @return Reference to this vector after addition.
	 */
	constexpr vec_3d<Value>& operator+=(const vec_3d<Value>& other) noexcept {
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	/**
	 * @brief Subtracts another vector from this vector.
	 *
	 * Performs component-wise subtraction.
	 *
	 * @param other The vector to subtract.
	 * @return Reference to this vector after subtraction.
	 */
	constexpr vec_3d<Value>& operator-=(const vec_3d<Value>& other) noexcept {
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

	/**
	 * @brief Scales this vector by a scalar.
	 *
	 * Multiplies each component by the provided scalar.
	 *
	 * @param scalar The scalar.
	 * @return Reference to this vector after scaling.
	 */
	constexpr vec_3d<Value>& operator*=(detail::arithmetic auto scalar) noexcept {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	/**
	 * @brief Computes the Euclidean norm of this vector.
	 *
	 * @return √(x² + y² + z²).
	 */
	CONSTEXPR_IF_GCC double euclidian_norm() const noexcept {
		return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
	}

	/**
	 * @brief Returns an iterator to the first component of this vector.
	 *
	 * Enables range-based for support.
	 *
	 * @return Pointer to the first component.
	 */
	constexpr Value* begin() noexcept {
		return &x;
	}

	/**
	 * @brief Returns a const iterator to the first component of this vector.
	 *
	 * Enables range-based for support on const objects.
	 *
	 * @return Const pointer to the first component.
	 */
	constexpr const Value* begin() const noexcept {
		return &x;
	}

	// NOLINTBEGIN(*pointer-arithmetic)
	/**
	 * @brief Returns an iterator one past the last component of this vector.
	 *
	 * Enables range-based for support.
	 *
	 * @return Pointer to one past the last component.
	 */
	constexpr Value* end() noexcept {
		return &x + 3;
	}

	/**
	 * @brief Returns a const iterator one past the last component of this vector.
	 *
	 * Enables range-based for support on const objects.
	 *
	 * @return Const pointer to one past the last component.
	 */
	constexpr const Value* end() const noexcept {
		return &x + 3;
	}

	// NOLINTEND(*pointer-arithmetic)

	/**
	 * @brief Returns a const pointer to the contigious components stored within this vector.
	 *
	 * @return Const pointer to the first component.
	 */
	constexpr const Value* data() const noexcept {
		return &x;
	}

	/**
	 * @brief Returns a pointer to the contigious components stored within this vector.
	 *
	 * @return Pointer to the first component.
	 */
	constexpr Value* data() noexcept {
		return &x;
	}
	// TODO(tuna): document and see if it makes the other operator unnecessary
	constexpr bool operator==(const vec_3d&) const = default;
};

/**
 * @brief Adds two 3D vectors.
 *
 * Performs component-wise addition.
 *
 * @tparam Value The type of the components of each vector.
 * @param lhs The left-hand side vector.
 * @param rhs The right-hand side vector.
 * @return A new vector equal to the sum of the two inputs.
 */
template <detail::arithmetic Value>
constexpr vec_3d<Value> operator+(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) noexcept {
	return vec_3d{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

/**
 * @brief Subtracts one 3D vector from another.
 *
 * Performs component-wise subtraction.
 *
 * @tparam Value The type of the components of each vector.
 * @param lhs The left-hand side vector.
 * @param rhs The right-hand side vector.
 * @return A new vector equal to the difference of the two inputs.
 */
template <detail::arithmetic Value>
constexpr vec_3d<Value> operator-(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) noexcept {
	return vec_3d{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

/**
 * @brief Multiplies a 3D vector by a scalar (scalar * vector).
 *
 * Performs component-wise scaling.
 *
 * @tparam Value The type of the components of @p vector.
 * @param scalar The scalar multiplier.
 * @param vector The vector to scale.
 * @return A new scaled vector.
 */
template <detail::arithmetic Value>
constexpr vec_3d<Value> operator*(detail::arithmetic auto scalar, const vec_3d<Value>& vector) noexcept {
	return vec_3d{vector.x * scalar, vector.y * scalar, vector.z * scalar};
}

/**
 * @brief Multiplies a 3D vector by a scalar (vector * scalar).
 *
 * Performs component-wise scaling.
 *
 * @tparam Value The type of the components of @p vector.
 * @param vector The vector to scale.
 * @param scalar The scalar multiplier.
 * @return A new scaled vector.
 */
template <detail::arithmetic Value>
constexpr vec_3d<Value> operator*(const vec_3d<Value>& vector, detail::arithmetic auto scalar) noexcept {
	return scalar * vector;
}

/**
 * @brief Compares two 3D vectors for equality.
 *
 * Performs component-wise comparison.
 *
 * @tparam Value The type of the components of each vector.
 * @param lhs The left-hand side vector.
 * @param rhs The right-hand side vector.
 * @return `true` @a iff all corresponding components are equal.
 */
template <detail::arithmetic Value>
constexpr bool operator==(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) {
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

/**
 * @brief Specialization of the {fmt} API for vec_3d.
 *
 * Enables use of fmt::format and fmt::print with vec_3d objects.
 *
 * @tparam Value The type of the components of the given vector.
 */
// NOLINTBEGIN(*convert-member-functions-to-static)
template <typename T>
struct fmt::formatter<vec_3d<T>> {
	/** @brief Parses the format specification (no-op for this type). */
	constexpr auto parse(fmt::format_parse_context& ctx) {
		return ctx.begin();
	}

	/** @brief Formats @p vector as "(x, y, z)". */
	auto format(const vec_3d<T>& vector, fmt::format_context& ctx) const {
		return fmt::format_to(ctx.out(), "({}, {}, {})", vector.x, vector.y, vector.z);
	}
};

// NOLINTEND(*convert-member-functions-to-static)

/**
 * @brief Stream insertion operator for vec_3d.
 *
 * Allows output of vectors using `operator<<`.
 *
 * @tparam Value The type of the components of @p vector.
 * @param stream Output stream to print into.
 * @param vector Vector to print.
 * @return Reference to the modified output stream.
 */
template <detail::arithmetic Value>
std::ostream& operator<<(std::ostream& stream, const vec_3d<Value>& vector) {
	return stream << fmt::format("{}", vector);
}

/**
 * @brief Type alias for a 3D vector of doubles.
 *
 * Provides a shorthand for `vec_3d<double>`.
 *
 * @typedef vec
 */
using vec = vec_3d<double>;

// Ensure that a 3D vector has the expected memory layout.
// NOLINTNEXTLINE(*avoid-c-arrays)
static_assert(sizeof(vec) == sizeof(double[3]));

#ifdef NDEBUG
static_assert(std::is_trivially_copyable_v<vec>);
#endif
