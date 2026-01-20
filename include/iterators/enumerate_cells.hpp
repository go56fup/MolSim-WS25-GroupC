#pragma once

#include "grid/particle_container/fwd.hpp"

namespace detail {

	// TODO(tuna): upgrade to std::random_iterator_tag
class enumerate_cells_iterator {
private:
	using inner = typename particle_container::value_type;

public:
	/// Standard iterator category.
	using iterator_category = std::input_iterator_tag;
	/// Type used for differences between indices.
	using difference_type = particle_container::difference_type;
	/// Type used for indexing.
	using size_type = particle_container::size_type;
	/// Reference to a cell with its index.
	using reference = std::pair<particle_container::index, particle_container::cell&>;
	/// The type returned from @ref operator*.
	using value_type = reference;

private:
	/// Pointer to container to enumerate cells for.
	particle_container* container = nullptr;
	/// Current cell index.
	particle_container::index idx{};

public:
	/// Default construct an enumerate cells iterator, does not reference any container.
	constexpr enumerate_cells_iterator() noexcept = default;

	/**
	 * @brief Constructs a new cell enumeration iterator.
	 *
	 * @tparam IndexT Deduced forwarding reference type for `begin_index`.
	 * @param container_ Reference to particle container to obtain interacting cell pairs from.
	 * @param begin_index Beginning cell index of iteration.
	 */
	template <fwd_reference_to<particle_container::index> IndexT>
	constexpr enumerate_cells_iterator(
		particle_container& container_, IndexT&& begin_index
	) noexcept
		: container(&container_)
		, idx(std::forward<IndexT>(begin_index)) {}

	/**
	 * @brief Dereferences this iterator.
	 * @return Pair over index of current cell and mutable reference to cell.
	 */
	constexpr value_type operator*() const noexcept {
		return {idx, (*container)[idx]};
	}

	/**
	 * @brief Advances the iterator to the next valid cell.
	 * @return Reference to this iterator after increment.
	 */
	constexpr enumerate_cells_iterator& operator++() noexcept {
		idx = next_3d_index(idx, container->grid_size());
		return *this;
	}

	/**
	 * @brief Post-increment operator.
	 * @return Copy of this enumeration iterator before increment.
	 */
	constexpr enumerate_cells_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	/**
	 * @brief Compares two cell enumeration iterators for equality.
	 * @param other Iterator to compare to.
	 * @return `true` @a iff both iterators reference the same cell.
	 */
	constexpr bool operator==(const enumerate_cells_iterator& other) const noexcept {
		assert(container == other.container);
		return idx == other.idx;
	}

	constexpr difference_type operator-(const enumerate_cells_iterator& other) const noexcept {
		assert(container == other.container);

		return static_cast<difference_type>(container->linear_index(idx.x, idx.y, idx.z)) -
		       static_cast<difference_type>(
				   container->linear_index(other.idx.x, other.idx.y, other.idx.z)
			   );
	}

	constexpr enumerate_cells_iterator& operator+=(difference_type advance) noexcept {
		const auto linear =
			static_cast<difference_type>(container->linear_index(idx.x, idx.y, idx.z)) + advance;

		idx = container->index_from_linear(static_cast<size_type>(linear));

		return *this;
	}
};

}  // namespace detail

/// Range over the enumerated cells of a particle container.
class enumerate_cells_range {
	// This range is meant to be non-owning.
	// NOLINTNEXTLINE(*avoid-*ref-data-members)
	/// Underlying container to enumerate.
	particle_container& container;

public:
	/// Type of iterator returned from @ref begin and @ref end.
	using iterator = detail::enumerate_cells_iterator;
	/// Type that represents the size of this range.
	using size_type = particle_container::size_type;

	/**
	 * @brief Constructs a new range over enumerated cells.
	 * @param c Reference to the particle container to iterate over.
	 */
	explicit constexpr enumerate_cells_range(particle_container& c) noexcept
		: container(c) {}

	/**
	 * @brief Returns an iterator to the first cell.
	 * @return Iterator to the beginning of the cell sequence.
	 */
	constexpr iterator begin() const noexcept {
		return iterator{container, particle_container::index{}};
	}

	/**
	 * @brief Returns the end iterator for the cells in the given particle contaienr.
	 * @return Iterator marking the end of the range.
	 */
	constexpr iterator end() const noexcept {
		return iterator{container, particle_container::index{0, 0, container.grid_size().z}};
	}

	/**
	 * @brief Returns the total number of cells in this range.
	 * @return The total number of cells in the particle container.
	 */
	constexpr size_type size() const noexcept {
		const auto& grid = container.grid_size();
		const auto x = grid.x;
		const auto y = grid.y;
		const auto z = grid.z;
		return 3 * x * y * z - x * z - x * y - y * z;
	}
};

static_assert(std::ranges::range<enumerate_cells_range>);
