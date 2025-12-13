#pragma once

#include <array>
#include <ranges>
#include <string>
#include <string_view>

#include <fmt/format.h>
// NOLINTBEGIN(*convert-member-functions-to-static)

enum class axis : std::uint8_t { x, y, z };

/**
 * @brief Specialization of the {fmt} API for the axis enum.
 *
 * Enables use of fmt::format and fmt::print with axis enumerators.
 *
 */

template <>
struct fmt::formatter<axis> {
private:
	static constexpr std::string_view to_string(axis a) noexcept {
		using enum axis;
		switch (a) {
		case x:
			return "x";
		case y:
			return "y";
		case z:
			return "z";
		}
		std::unreachable();
	}

public:
	/** @brief Parses the format specification (no-op for this type). */
	constexpr auto parse(fmt::format_parse_context& ctx) {
		return ctx.begin();
	}

	/** @brief Formats `Axis` objects". */
	auto format(axis a, fmt::format_context& ctx) const {
		return fmt::format_to(ctx.out(), "{}", to_string(a));
	}
};

enum class boundary_type : std::uint8_t {
	x_min = 1u << 0u,  // 1
	y_min = 1u << 1u,  // 2
	z_min = 1u << 2u,  // 4
	x_max = 1u << 3u,  // 8
	y_max = 1u << 4u,  // 16
	z_max = 1u << 5u   // 32
};

constexpr boundary_type operator|(boundary_type lhs, boundary_type rhs) noexcept {
	using T = std::underlying_type_t<boundary_type>;
	return static_cast<boundary_type>(static_cast<T>(lhs) | static_cast<T>(rhs));
}

constexpr boundary_type operator&(boundary_type lhs, boundary_type rhs) noexcept {
	using T = std::underlying_type_t<boundary_type>;
	return static_cast<boundary_type>(static_cast<T>(lhs) & static_cast<T>(rhs));
}

constexpr boundary_type operator^(boundary_type lhs, boundary_type rhs) noexcept {
	using T = std::underlying_type_t<boundary_type>;
	return static_cast<boundary_type>(static_cast<T>(lhs) ^ static_cast<T>(rhs));
}

constexpr boundary_type operator~(boundary_type v) noexcept {
	using T = std::underlying_type_t<boundary_type>;
	return static_cast<boundary_type>(~static_cast<T>(v));
}

constexpr boundary_type& operator|=(boundary_type& lhs, boundary_type rhs) noexcept {
	lhs = lhs | rhs;
	return lhs;
}

constexpr boundary_type& operator&=(boundary_type& lhs, boundary_type rhs) noexcept {
	lhs = lhs & rhs;
	return lhs;
}

inline constexpr std::array<std::array<boundary_type, 2>, 3> border_pairs{
	{{boundary_type::x_min, boundary_type::x_max},
     {boundary_type::y_min, boundary_type::y_max},
     {boundary_type::z_min, boundary_type::z_max}}
};

// TODO(tuna): make this docstring an alias
/**
 * @brief Specialization of the {fmt} API for the `boundary_type` enum.
 *
 * Enables use of fmt::format and fmt::print with axis enumerators.
 *
 */
template <>
struct fmt::formatter<boundary_type> {
private:
public:
	/** @brief Parses the format specification (no-op for this type). */
	constexpr auto parse(fmt::format_parse_context& ctx) {
		return ctx.begin();
	}

	/** @brief Formats `boundary_type` enumerators". */
	constexpr auto format(boundary_type b, fmt::format_context& ctx) const {
		using enum boundary_type;
		bool first = true;
		auto out = ctx.out();

		auto emit = [&](std::string_view name) {
			if (first) {
				out = fmt::format_to(out, "{}", name);
				first = false;
			} else {
				out = fmt::format_to(out, ", {}", name);
			}
		};

		if ((b & x_min) == x_min) emit("x_min");
		if ((b & x_max) == x_max) emit("x_max");
		if ((b & y_min) == y_min) emit("y_min");
		if ((b & y_max) == y_max) emit("y_max");
		if ((b & z_min) == z_min) emit("z_min");
		if ((b & z_max) == z_max) emit("z_max");
		return out;
	}
};

// NOLINTEND(*convert-member-functions-to-static)
