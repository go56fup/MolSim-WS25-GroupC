#pragma once

#include <ranges>

#include "grid/particle_container/fwd.hpp"
#include "utility/tracing/macros.hpp"

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
	particle_container* container = nullptr;
	particle_container::signed_index current_virtual_idx;
	particle_container::signed_index target_cell_idx;
	std::uint8_t displacement_idx = 0;
	particle_container::signed_index signed_grid;
	// TODO(tuna): because displacement_idx is incremented unconditionally, it was never actually
	// supposed to be used to detect the .end(). However, the only state that is advanced in the
	// current implementation (which stops yielding after a single cell and never gets the next 3D
	// index) is the displacement_idx. Therefore, with an unconditional increment, there is no way
	// to tell the last element from the .end(), as if the .end() has displacement_idx + 1,
	// operator++ is not prevented from indexing into displacement_count which is 1-past-the-end.
	// The unconditional increment would need to be refactored to get rid of this extra
	// statekeeping.
	bool exhausted = false;
	static constexpr std::array<particle_container::signed_index, 13> displacements = {
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

	/// Number of registered displacements.
	static constexpr std::size_t displacement_count = std::tuple_size_v<decltype(displacements)>;

public:
	/// Tag type used to directly obtain the end iterator of a container.
	struct get_end_tag {};

	constexpr periodic_iterator() noexcept = default;

	// Start Constructor
	constexpr periodic_iterator(
		particle_container& c, const particle_container::signed_index& v_idx
	)
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

	constexpr periodic_iterator(particle_container& c, periodic_iterator::get_end_tag)
		: container(&c)
		, displacement_idx(displacement_count)
		, exhausted(true) {}

	template <axis Axis>
	constexpr bool do_displacement(const particle_container::signed_index& displacement) noexcept {
		// We do not need to use __builtin_add_overflow() here, because:
		assert(displacement[Axis] == -1 || displacement[Axis] == 1 || displacement[Axis] == 0);
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
			TRACE_PERIODIC("Doing displacement {} on {}", displacement, current_virtual_idx);

			const bool success = do_displacement<axis::x>(displacement) &&
			                     do_displacement<axis::y>(displacement) &&
			                     do_displacement<axis::z>(displacement);

			++displacement_idx;

			if (success) {
				TRACE_PERIODIC(
					"Found valid interaction: {} -> {}", current_virtual_idx, target_cell_idx
				);
				return *this;
			}
		}
		exhausted = true;
		return *this;
	}

	constexpr periodic_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	constexpr bool operator==(const periodic_iterator& other) const noexcept {
		assert(container == other.container);
		return displacement_idx == other.displacement_idx && exhausted == other.exhausted;
	}

	constexpr particle_container::index operator*() const noexcept {
		return {
			static_cast<size_type>(target_cell_idx.x), static_cast<size_type>(target_cell_idx.y),
			static_cast<size_type>(target_cell_idx.z)
		};
	}
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
		return {container, detail::periodic_iterator::get_end_tag{}};
	}
};

static_assert(std::ranges::range<periodic_range>);
