#pragma once
#include <type_traits>

#include "grid/particle_container/fwd.hpp"
#include "utility/tracing/macros.hpp"

namespace detail {
/// Iterator for obtaining the minimum necessary pairs of cells that are to interact via the linked
/// cell algorithm.
class interactions_iterator {
private:
	/// Element type of the underlying container.
	// It doesn't make much sense to have a const particle_container, so this iterator does not
	// concern itself with const correctness.
	// TODO(tuna): look into the above
	using inner = particle_container::value_type;

public:
	/// Standard iterator category.
	using iterator_category = std::forward_iterator_tag;
	/// Type used for indexing.
	using size_type = particle_container::size_type;
	/// Type used for differences between indices.
	using difference_type = particle_container::difference_type;
	/// Reference to the indices of a pair of interacting cells, denoting the state of this
	/// iterator.
	using reference = std::pair<const particle_container::index&, const particle_container::index&>;
	/// The type returned from @ref operator*, which is a reference.
	// TODO(tuna): consider decoupling reference from value_type and just returning the indices
	// directly.
	using value_type = reference;

private:
	/// Pointer to container which the iterator obtains the cell pairs from.
	particle_container* container = nullptr;
	/// Index of the first cell of the returned interacting pair.
	particle_container::index current_cell_idx;
	/// Index of the second cell of the returned interacting pair.
	particle_container::index target_cell_idx;
	/// Index into @ref displacements describing the progress of the scan of cells to interact with.
	std::uint8_t displacement_idx = 0;

	/// Incremental differences that cover all intercellular interactions in the least pairs by
	/// scanning to the bottom-right.
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
	static_assert(
		displacement_count * 2 + 1 == 3uz * 3 * 3,
		"Covers all adjacent cells (3x3x3 cube around current_cell_idx)"
	);

public:
	/// Tag type used to directly obtain the end iterator of a container.
	struct get_end_tag {};

	/// Default construct an interactions iterator, does not reference any container.
	constexpr interactions_iterator() noexcept = default;

	/**
	 * @brief Constructs a new interactions iterator.
	 *
	 * @tparam IndexT Deduced forwarding reference type for `begin_index`.
	 * @param container_ Reference to particle container to obtain interacting cell pairs from.
	 * @param begin_index Beginning cell index of iteration.
	 */
	template <fwd_reference_to<particle_container::index> IndexT>
	constexpr interactions_iterator(particle_container& container_, IndexT&& begin_index) noexcept
		: container(&container_)
		, current_cell_idx(std::forward<IndexT>(begin_index)) {
		// Acquire target.
		++*this;
	}

	/**
	 * @brief Constructs the end iterator for interacting cell pairs for a container.
	 * @param container_ Reference to the underlying particle container.
	 */
	constexpr interactions_iterator(particle_container& container_, get_end_tag) noexcept
		: container(&container_)
		, current_cell_idx{0, 0, container_.grid_size().z} {}

	/**
	 * @brief Dereferences this iterator.
	 * @return Pair of mutable references to the current elements.
	 */
	constexpr value_type operator*() noexcept {
		return reference{current_cell_idx, target_cell_idx};
	}

	/**
	 * @brief Performs the given displacement on @ref current_cell_idx and returns if the operation
	 * was legal.
	 *
	 * @param displacement Vector describing an incremental shift on the x, y and z axes.
	 * @return `true` @a iff the displacement resulted in a cell coordinate inside the domain.
	 */
	template <axis Axis>
	constexpr bool do_displacement(const particle_container::signed_index& displacement) noexcept {
		if (displacement[Axis] == 0) return true;
		TRACE_INTERACTION_ITER("Doing displacement {} on {}", displacement, current_cell_idx);

		// We do not need to use __builtin_add_overflow() here, because:
		assert(displacement[Axis] == -1 || displacement[Axis] == 1);
		const bool outside_domain_bc_negative =
			displacement[Axis] < 0 && current_cell_idx[Axis] == 0;
		if (outside_domain_bc_negative) return false;
		target_cell_idx[Axis] = apply_difference(current_cell_idx[Axis], displacement[Axis]);
		const bool outside_domain_bc_max =
			out_of_bounds<axis_to_boundary_type(Axis, extremum::max)>(
				target_cell_idx, container->grid_size()
			);
		return !outside_domain_bc_max;
	}

	/**
	 * @brief Advances the iterator to the next valid pair of interacting cells.
	 * @return Reference to this iterator after increment.
	 */
	constexpr interactions_iterator& operator++() noexcept {
		const auto& grid = container->grid_size();

		while (current_cell_idx.z < grid.z) {
			// Try all displacements for the current cell
			while (displacement_idx < displacement_count) {
				const auto& displacement = displacements[displacement_idx];
				target_cell_idx = current_cell_idx;

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

			// No displacement succeeded: advance to the next cell
			displacement_idx = 0;
			current_cell_idx = next_3d_index(current_cell_idx, grid);
			target_cell_idx = current_cell_idx;
		}
		// This is the end iterator
		assert(
			current_cell_idx.x == 0 && current_cell_idx.y == 0 && current_cell_idx.z == grid.z &&
			displacement_idx == 0
		);
		return *this;
	}

	/**
	 * @brief Post-increment operator.
	 * @return Copy of this interaction iterator before increment.
	 */
	constexpr interactions_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	/**
	 * @brief Compares two interaction iterators for equality.
	 * @param other Iterator to compare to.
	 * @return `true` @a iff both iterators reference the same current cell and displacement.
	 */
	constexpr bool operator==(const interactions_iterator& other) const noexcept {
		assert(container == other.container);
		return current_cell_idx == other.current_cell_idx &&
		       displacement_idx == other.displacement_idx;
	}
};

}  // namespace detail

/// Range over interacting cell indices, implementing the linked cell algorithm.
class interactions_range {
	// This range is meant to be non-owning.
	// NOLINTNEXTLINE(*avoid-*ref-data-members)
	/// Underlying container to calculate interacting cells for.
	particle_container& container;

public:
	/// Type of iterator returned from @ref begin and @ref end.
	using iterator = detail::interactions_iterator;
	/// Type that represents the size of this range.
	using size_type = std::size_t;

	/**
	 * @brief Constructs a new range over pairs of interacting cells.
	 * @param c Reference to the particle container to iterate over.
	 */
	explicit constexpr interactions_range(particle_container& c) noexcept
		: container(c) {}

	/**
	 * @brief Returns an iterator to the first pair of cells to interact.
	 * @return Iterator to the beginning of the cell pair sequence.
	 */
	constexpr iterator begin() noexcept {
		return {container, particle_container::index{}};
	}

	/**
	 * @brief Returns the end iterator for interacting cell pairs for a container.
	 * @return Iterator marking the end of the range.
	 */
	constexpr iterator end() noexcept {
		return iterator{container, iterator::get_end_tag{}};
	}

	/**
	 * @brief Returns the total number of pairs that would be generated by this range.
	 * @return The total number of pairs in range.
	 */
	constexpr size_type size() const noexcept {
		const auto& grid = container.grid_size();
		return grid.x * grid.y * (3 * grid.z - 1) - grid.z * (grid.x + grid.y);
	}
};

static_assert(std::ranges::range<interactions_range>);
