#pragma once

#include <ranges>

#include "grid/particle_container/fwd.hpp"


namespace detail {
class periodic_iterator {
public:
	using iterator_category = std::forward_iterator_tag;
	using value_type = particle_container::index;
	using size_type = particle_container::size_type;
	using difference_type = particle_container::difference_type;
	using pointer = value_type*;
	using reference = value_type;

private:
	using signed_index = vec_3d<difference_type>;

	particle_container* container = nullptr;
	signed_index current_virtual_idx{};
	signed_index target_cell_idx{};
	std::uint8_t displacement_idx = 0;
	signed_index signed_grid{};

	static constexpr std::uint8_t displacement_count = 13;

public:
	constexpr periodic_iterator() noexcept = default;

	// Start Constructor
	constexpr periodic_iterator(particle_container& c, const signed_index& v_idx)
		: container(&c)
		, current_virtual_idx(v_idx)
		, target_cell_idx(v_idx)
		, signed_grid{
			  static_cast<difference_type>(c.grid_size().x),
			  static_cast<difference_type>(c.grid_size().y),
			  static_cast<difference_type>(c.grid_size().z)
		  } {
		++*this;
	}

	constexpr periodic_iterator(
		particle_container& c, const signed_index& begin_idx, interactions_iterator::get_end_tag
	)
		: container(&c)
		, current_virtual_idx{begin_idx.x, 0, signed_grid.z} {}

	template <axis Axis>
	constexpr bool do_displacement(const signed_index& displacement) noexcept {
		if (displacement[Axis] == 0) return true;
		TRACE_INTERACTION_ITER("Doing displacement {} on {}", displacement, current_virtual_idx);

		// We do not need to use __builtin_add_overflow() here, because:
		assert(displacement[Axis] == -1 || displacement[Axis] == 1);
		target_cell_idx[Axis] = current_virtual_idx[Axis] + displacement[Axis];

		// If we are moving backwards, don't go beyond the first ghost layer
		const bool outside_domain_bc_negative = target_cell_idx[Axis] < 0;
		if (outside_domain_bc_negative) return false;

		const bool outside_domain_bc_max =
			out_of_bounds<axis_to_boundary_type(Axis, extremum::max)>(target_cell_idx, signed_grid);

		return !outside_domain_bc_max;
	}

	constexpr periodic_iterator& operator++() noexcept {
		// Try all displacements for the current cell
		while (displacement_idx < displacement_count) {
			const auto& displacement = displacements[displacement_idx];
			target_cell_idx = current_virtual_idx;

			const bool success = do_displacement<axis::x>(displacement) &&
			                     do_displacement<axis::y>(displacement) &&
			                     do_displacement<axis::z>(displacement);

			++displacement_idx;

			if (success) {
				TRACE_INTERACTION_ITER(
					"Found valid interaction: {} -> {}", current_cell_idx, target_cell_idx
				);
				return *this;
			}
		}

		return *this;
	}

	constexpr periodic_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	constexpr bool operator==(const periodic_iterator& other) const noexcept {
		assert(container == other.container);
		return displacement_idx == other.displacement_idx;
	}

	constexpr particle_container::index operator*() const noexcept {
		return {
			static_cast<size_type>(target_cell_idx.x), static_cast<size_type>(target_cell_idx.y),
			static_cast<size_type>(target_cell_idx.z)
		};
	}

	static constexpr std::array<signed_index, 13> displacements = {
		{{0, 0, +1},    // i,     j,     k + 1
	     {0, +1, -1},   // i,     j + 1, k - 1
	     {0, +1, 0},    // i,     j + 1, k
	     {0, +1, +1},   // i,     j + 1, k + 1
	     {1, -1, -1},   // i + 1, j - 1, k - 1
	     {+1, -1, 0},   // i + 1, j - 1, k
	     {+1, -1, +1},  // i + 1, j - 1, k + 1
	     {+1, 0, -1},   // i + 1, j,     k - 1
	     {+1, 0, 0},    // i + 1, j,     k
	     {+1, 0, +1},   // i + 1, j,     k + 1
	     {+1, +1, -1},  // i + 1, j + 1, k - 1
	     {+1, +1, 0},   // i + 1, j + 1, k
	     {+1, +1, +1}}  // i + 1, j + 1, k + 1
	};
};
}  // namespace detail

/**
 * @brief A lightweight range for periodic neighbor interactions.
 * Input: A signed virtual index.
 * Output: Yields unsigned target indices.
 */
class periodic_range {
	particle_container& container;
	vec_3d<particle_container::difference_type> begin_idx;

public:
	constexpr periodic_range(
		particle_container& c, vec_3d<particle_container::difference_type> start_v_idx
	) noexcept
		: container(c)
		, begin_idx(start_v_idx) {}

	constexpr detail::periodic_iterator begin() noexcept {
		return {container, begin_idx};
	}

	constexpr detail::periodic_iterator end() noexcept {
		return {container, begin_idx, detail::interactions_iterator::get_end_tag{}};
	}
};

static_assert(std::ranges::range<periodic_range>);
