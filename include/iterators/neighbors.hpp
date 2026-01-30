#pragma once

#include "grid/bounds/operations.hpp"
#include "grid/particle_container/fwd.hpp"

namespace detail {

class neighbors_iterator {
public:
	using iterator_category = std::forward_iterator_tag;
	using size_type = particle_container::size_type;
	using difference_type = particle_container::difference_type;
	using reference = const particle_container::index&;
	using value_type = particle_container::index;

private:
	particle_container* container = nullptr;
	particle_container::index center_cell;
	particle_container::index neighbor_cell;
	std::uint8_t displacement_idx = 0;
	bool exhausted = false;

	static constexpr std::array<particle_container::signed_index, 26> displacements = {{
		{-1, -1, -1}, {-1, -1, 0},  {-1, -1, +1}, {-1, 0, -1},  {-1, 0, 0},
		{-1, 0, +1},  {-1, +1, -1}, {-1, +1, 0},  {-1, +1, +1},

		{0, -1, -1},  {0, -1, 0},   {0, -1, +1},  {0, 0, -1},   {0, 0, +1},
		{0, +1, -1},  {0, +1, 0},   {0, +1, +1},

		{+1, -1, -1}, {+1, -1, 0},  {+1, -1, +1}, {+1, 0, -1},  {+1, 0, 0},
		{+1, 0, +1},  {+1, +1, -1}, {+1, +1, 0},  {+1, +1, +1},
	}};
	static constexpr std::size_t displacement_count = std::tuple_size_v<decltype(displacements)>;

public:
	struct get_end_tag {};

	constexpr neighbors_iterator() noexcept = default;

	constexpr neighbors_iterator(
		particle_container& c, const particle_container::index& center
	) noexcept
		: container(&c)
		, center_cell(center) {
		++(*this);  // advance to first valid neighbor
	}

	constexpr neighbors_iterator(
		particle_container& c, const particle_container::index& center, get_end_tag
	) noexcept
		: container(&c)
		, center_cell(center)
		, displacement_idx(displacement_count)
		, exhausted(true) {}

	constexpr reference operator*() const noexcept {
		return neighbor_cell;
	}

	template <axis Axis>
	constexpr bool do_displacement(const particle_container::signed_index& displacement) noexcept {
		if (displacement[Axis] == 0) {
			neighbor_cell[Axis] = center_cell[Axis];
			return true;
		}

		assert(displacement[Axis] == -1 || displacement[Axis] == 1);

		if (displacement[Axis] < 0 && center_cell[Axis] == 0) return false;

		neighbor_cell[Axis] = apply_difference(center_cell[Axis], displacement[Axis]);

		return !out_of_bounds<axis_to_boundary_type(Axis, extremum::max)>(
			neighbor_cell, container->grid_size()
		);
	}

	constexpr neighbors_iterator& operator++() noexcept {
		while (displacement_idx < displacement_count) {
			const auto& d = displacements[displacement_idx++];
			neighbor_cell = center_cell;

			const bool success = do_displacement<axis::x>(d) && do_displacement<axis::y>(d) &&
			                     do_displacement<axis::z>(d);

			if (success) {
				return *this;
			}
		}

		// end iterator
		displacement_idx = displacement_count;
		exhausted = true;
		return *this;
	}

	constexpr neighbors_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	constexpr bool operator==(const neighbors_iterator& other) const noexcept {
		assert(container == other.container);
		return displacement_idx == other.displacement_idx && center_cell == other.center_cell &&
		       exhausted == other.exhausted;
	}
};

}  // namespace detail

class neighbors_range {
	particle_container& container;
	particle_container::index center;

public:
	using iterator = detail::neighbors_iterator;
	using size_type = std::size_t;

	constexpr neighbors_range(
		particle_container& c, const particle_container::index& center_cell
	) noexcept
		: container(c)
		, center(center_cell) {}

	constexpr iterator begin() noexcept {
		return iterator{container, center};
	}

	constexpr iterator end() noexcept {
		return iterator{container, center, iterator::get_end_tag{}};
	}

	/// Maximum possible neighbors (interior cell)
	static constexpr size_type max_size() noexcept {
		return 26;
	}
};

static_assert(std::ranges::range<neighbors_range>);
