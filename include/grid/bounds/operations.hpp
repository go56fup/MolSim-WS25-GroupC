#pragma once
#include <cstdint>

#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
#include "physics/forces.hpp"
#include "physics/particle.hpp"

template <boundary_type border, typename IndexT>
constexpr bool out_of_bounds(
	const vec_3d<IndexT>& pos, const vec_3d<IndexT>& bounds, const vec_3d<IndexT>& min = {0, 0, 0}
) noexcept {
	using enum boundary_type;
	static constexpr bool is_min = boundary_type_to_extremum(border) == extremum::min;
	static constexpr axis ax = boundary_type_to_axis(border);
	const auto oob = is_min ? pos[ax] <= min[ax] : pos[ax] >= bounds[ax];
	return oob;
}

// Swapping the parameters does not change the meaning of the operation.
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
template <std::unsigned_integral Unsigned>
constexpr Unsigned
apply_difference(Unsigned unsigned_, std::make_signed_t<Unsigned> signed_) noexcept {
	return static_cast<Unsigned>(static_cast<decltype(signed_)>(unsigned_) + signed_);
}

template <std::unsigned_integral Unsigned>
constexpr Unsigned
apply_difference(std::make_signed_t<Unsigned> signed_, Unsigned unsigned_) noexcept {
	return apply_difference(unsigned_, signed_);
}

// NOLINTEND(bugprone-easily-swappable-parameters)

template <std::unsigned_integral Unsigned>
constexpr vec_3d<Unsigned>
// Putting these in a struct breaks the deduction of the template parameter, and specifying
// `f(args<particle_container::size_type>{.current = i, .bounds = grid})` is way too verbose.
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
next_3d_index(const vec_3d<Unsigned>& current, const vec_3d<Unsigned>& bounds) noexcept {
	if (current.x + 1 < bounds.x) {
		return {current.x + 1, current.y, current.z};
	}

	if (current.y + 1 < bounds.y) {
		return {0, current.y + 1, current.z};
	}
	assert(current.z + 1 <= bounds.z);
	return {0, 0, current.z + 1};
}
