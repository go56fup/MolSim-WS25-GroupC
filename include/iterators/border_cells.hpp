#pragma once

#include "grid/particle_container/fwd.hpp"
#include "utility/tracing/macros.hpp"

namespace detail {
/// Iterator for obtaining the cells at the borders of a particle container.
class border_cell_iterator {
private:
	/// Element type of the underlying container.
	using inner = typename particle_container::value_type;

public:
	/// Standard iterator category.
	using iterator_category = std::forward_iterator_tag;
	/// Type used for differences between indices.
	using difference_type = particle_container::difference_type;
	/// Type used for indexing.
	using size_type = particle_container::size_type;
	// TODO(tuna): make value_type's across the codebase const compliant
	using reference = std::pair<particle_container::index, boundary_type>;
	// TODO(tuna): see if returning index behavior is better behaved than const& to index
	/// The type returned from @ref operator*.
	using value_type = reference;

private:
	/// Pointer to container which the iterator obtains the border cells from.
	particle_container* container = nullptr;
	/// Index of the current cell of iteration.
	particle_container::index idx;
	/// Boundary type of current cell.
	boundary_type type{};

	/**
	   @brief Determines on which boundaries of the grid a cell lies.
	   @param current Index to cell to query.
	   @return Enum representing the cell's neighbourship to the grid borders.
	 */
	constexpr boundary_type
	determine_border_type(const particle_container::index& current) const noexcept {
		const auto& boundaries = container->grid_size();
		// TODO(tuna): find a way to elide these checks:
		// for the 1-cell thick case, the cell becomes both *_min and *_max; which only
		// happens there but causes the overhead of both checks for the majority of grid
		// configurations
		using enum boundary_type;
		boundary_type result{};
		if (current.x == 0) {
			result |= x_min;
		}
		if (current.x == boundaries.x - 1) {
			result |= x_max;
		}
		if (current.y == 0) {
			result |= y_min;
		}
		if (current.y == boundaries.y - 1) {
			result |= y_max;
		}
		if (current.z == 0) {
			result |= z_min;
		}
		if (current.z == boundaries.z - 1) {
			result |= z_max;
		}
		TRACE_BORDER_CELL_ITER("Border type of {}: {}", current, result);
		return result;
	};

public:
	/// Default construct a border cells iterator, does not reference any container.
	constexpr border_cell_iterator() noexcept = default;

	/**
	 * @brief Constructs a new interactions iterator.
	 *
	 * @tparam IndexT Deduced forwarding reference type for `begin_index`.
	 * @param container_ Reference to particle container to obtain border cells from.
	 * @param begin_index Beginning cell index of iteration.
	 */
	template <fwd_reference_to<particle_container::index> IndexT>
	constexpr border_cell_iterator(particle_container& c, IndexT&& begin_index) noexcept
		: container(&c)
		, idx(std::forward<IndexT>(begin_index))
		, type(determine_border_type(idx)) {}

	/**
	 * @brief Dereferences this iterator.
	 * @return Pair over the index of the current border cell and its type.
	 */
	constexpr value_type operator*() const noexcept {
		return {idx, type};
	}

	/**
	 * @brief Advances the iterator to the next border cell.
	 * @return Reference to this iterator after increment.
	 */
	constexpr border_cell_iterator& operator++() noexcept {
		const auto boundaries = container->grid_size();

		while (true) {
			idx = next_3d_index(idx, boundaries);
			type = determine_border_type(idx);
			if (std::to_underlying(type) == 0) {
				continue;
			}
			return *this;
		}
	}

	/**
	 * @brief Post-increment operator.
	 * @return Copy of this border cell iterator before increment.
	 */
	constexpr border_cell_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}


	/**
	 * @brief Compares two border cell iterators for equality.
	 * @param other Iterator to compare to.
	 * @return `true` @a iff both iterators reference the same border cell.
	 */
	constexpr bool operator==(const border_cell_iterator& other) const noexcept {
		assert(container == other.container);
		return idx == other.idx;
	}
};

}  // namespace detail

/// Range over the border cells of a particle container.
class border_cell_range {
	// This range is meant to be non-owning.
	// NOLINTNEXTLINE(*avoid-*ref-data-members)
	/// Underlying container to pull border cells from.
	particle_container& container;

public:
	/// Type of iterator returned from @ref begin and @ref end.
	using iterator = detail::border_cell_iterator;
	/// Type that represents the size of this range.
	using size_type = particle_container::size_type;

	/**
	 * @brief Constructs a new range over pairs of border cells.
	 * @param c Reference to the particle container to iterate over.
	 */
	explicit constexpr border_cell_range(particle_container& c) noexcept
		: container(c) {}

	/**
	 * @brief Returns an iterator to the first border cell.
	 * @return Iterator to the beginning of the border cell sequence.
	 */
	constexpr iterator begin() const noexcept {
		using index = particle_container::index;
		return iterator{container, index{}};
	}

	/**
	 * @brief Returns the end iterator for border cells.
	 * @return Iterator marking the end of the range.
	 */
	constexpr iterator end() const noexcept {
		return iterator{container, particle_container::index{0, 0, container.grid_size().z}};
	}

	/**
	 * @brief Returns the total number of pairs that would be generated by this range.
	 * @return The total number of pairs in range.
	 */
	constexpr size_type size() const noexcept {
		const auto& grid = container.grid_size();
		const auto x = grid.x;
		const auto y = grid.y;
		const auto z = grid.z;
		return 2 * x * z + 2 * y * z - 4 * x - 4 * y - 4 * z + 8;
	}
};
static_assert(std::ranges::range<border_cell_range>);
