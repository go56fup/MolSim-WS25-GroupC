#pragma once

#include <array>
#include <cmath>
#include <format>
#include <span>

namespace detail {
template <typename T>
concept numeric = std::floating_point<T> || std::integral<T>;
}

template <detail::numeric Value>
struct vec_3d {
	Value x;
	Value y;
	Value z;

	constexpr auto operator+=(const vec_3d<Value>& other) noexcept {
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	constexpr auto operator-=(const vec_3d<Value>& other) noexcept {
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

	constexpr auto operator*=(detail::numeric auto scalar) noexcept {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	double euclidian_norm() const noexcept {
		return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
	}

	constexpr auto begin() noexcept {
		return &x;
	}

	constexpr auto begin() const noexcept {
		return &x;
	}

	constexpr auto end() noexcept {
		return &x + 3;
	}

	constexpr auto end() const noexcept {
		return &x + 3;
	}

	constexpr const Value* data() const noexcept {
		return &x;
	}

	constexpr Value* data() noexcept {
		return &x;
	}
};

template <detail::numeric Value>
constexpr auto operator+(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) noexcept {
	return vec_3d{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

template <detail::numeric Value>
constexpr auto operator-(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) noexcept {
	return vec_3d{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

template <detail::numeric Value>
constexpr auto operator*(detail::numeric auto scalar, const vec_3d<Value>& vector) noexcept {
	return vec_3d{vector.x * scalar, vector.y * scalar, vector.z * scalar};
}

template <detail::numeric Value>
constexpr auto operator*(const vec_3d<Value>& vector, detail::numeric auto scalar) noexcept {
	return vec_3d{vector.x * scalar, vector.y * scalar, vector.z * scalar};
}

template <detail::numeric Value>
std::ostream& operator<<(std::ostream& stream, const vec_3d<Value>& vector) {
	return stream << std::format("({}, {}, {})", vector.x, vector.y, vector.z);
}

template <detail::numeric Value>
constexpr bool operator==(const vec_3d<Value>& lhs, const vec_3d<Value>& rhs) {
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

using vec = vec_3d<double>;
static_assert(sizeof(vec) == sizeof(double[3]));
