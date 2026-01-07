#pragma once

#include <ranges>
#include <type_traits>

namespace detail {
/**
 * @brief Helper class to generate pairwise iteration ranges over a container.
 *
 * Provides mechanisms for iterating over all possible pairs of elements in a container,
 * either **unique pairs** via @ref uniques (i.e., `(i, j)` where `i < j`) or **non-unique pairs**
 * via @ref pairs (i.e., all ordered pairs `(i, j)` where `i != j`).
 *
 * @tparam Container A random-access container type (e.g., `std::vector`, `std::array`).
 */
template <std::ranges::random_access_range Container>
class pairwise {
private:
	/**
	 * @brief Enumeration to distinguish between unique and non-unique pair iteration modes.
	 */
	enum class Uniqueness : std::uint8_t { unique, non_unique };

	/**
	 * @brief Iterator for generating pairs of elements in the container.
	 *
	 * Depending on the @ref Uniqueness mode, the iterator produces:
	 * - Unique pairs: `(0,1), (0,2), (1,2), ...`
	 * - Non-unique pairs: `(0,1), (0,2), (1,0), (1,2), ...`
	 *
	 * @tparam mode The uniqueness mode (unique or non-unique).
	 */
	template <Uniqueness mode>
	class p_iterator {
	private:
		/// The type of the underlying container.
		// Allowing for templating over references makes our job easier when using the
		// implementation.
		using container_t = std::remove_reference_t<Container>;
		/// Element type of the container.
		using inner_wo_quals = typename container_t::value_type;
		/// Adjusted element type respecting the const-qualification of the container.
		// Switching over the container's type here allows us to have a singular source of truth in
		// terms of const-correctness later.
		using inner =
			std::conditional_t<std::is_const_v<container_t>, const inner_wo_quals, inner_wo_quals>;

	public:
		/// Standard iterator category.
		using iterator_category = std::forward_iterator_tag;
		/// Type used for differences between indices.
		using difference_type = container_t::difference_type;
		/// Type used for indexing.
		using size_type = container_t::size_type;
		/// Reference to two elements of the container.
		using reference = std::pair<inner&, inner&>;
		/// The type returned from @ref operator*, which is a reference.
		using value_type = reference;

	private:
		/// Pointer to container which the iterator obtains the pairs from.
		container_t* container = nullptr;
		/// Index of the first element of the returned pair.
		size_type outer_idx = 0;
		/// Index of the second element of the returned pair.
		size_type inner_idx = 0;

	public:
		/// Default construct a pairwise iterator, does not reference any container.
		constexpr p_iterator() noexcept = default;

		/**
		 * @brief Constructs a new pairwise iterator.
		 *
		 * @param c Reference to the underlying container.
		 * @param first_start Initial index for the first element of the pair.
		 * @param second_start Initial index for the second element of the pair.
		 */
		constexpr p_iterator(Container& c, size_type first_start, size_type second_start) noexcept
			: container(&c)
			, outer_idx(first_start)
			, inner_idx(second_start) {}

		/**
		 * @brief Dereferences this iterator.
		 * @return Pair of references to the current elements, mutable if the underlying container
		 * is mutable.
		 */

		constexpr value_type operator*() const noexcept {
			return reference{(*container)[outer_idx], (*container)[inner_idx]};
		}

		/**
		 * @brief Advances the iterator to the next valid pair.
		 * @return Reference to this iterator after increment.
		 */
		constexpr p_iterator& operator++() noexcept {
			const auto size = std::ranges::size(*container);

			const auto check_and_do_outer_inc = [&] {
				// e.g. for {0, 1, 2}.
				// if we have reached the end of the inner loop:
				if (inner_idx >= size) {
					++outer_idx;
					if constexpr (mode == Uniqueness::unique) {
						// (0, 1); (0, 2); (1, 2) -> set inner to outer + 1
						// because we have seen the rest
						inner_idx = outer_idx + 1;
					} else {
						// (0, 1); (0, 2); (1, 0) -> set inner to 0, start over
						inner_idx = 0;
					}
				}
			};

			++inner_idx;
			check_and_do_outer_inc();

			if constexpr (mode == Uniqueness::non_unique) {
				// (outer == inner) not automatically handled here, handle that
				if (outer_idx < size && outer_idx == inner_idx) {
					++inner_idx;
					check_and_do_outer_inc();
				}
			}

			return *this;
		}

		/**
		 * @brief Post-increment operator.
		 * @return Copy of this pairwise iterator before increment.
		 */
		constexpr p_iterator operator++(int) noexcept {
			auto tmp = *this;
			++(*this);
			return tmp;
		}

		/**
		 * @brief Compares two iterators for equality.
		 * @param other Iterator to compare to.
		 * @return `true` @a iff both iterators refer to the same position in the same container.
		 */
		constexpr bool operator==(const p_iterator& other) const noexcept {
			return container == other.container && outer_idx == other.outer_idx &&
			       inner_idx == other.inner_idx;
		}
	};

	/**
	 * @brief Range that yields pairs of elements based on the uniqueness mode.
	 *
	 * Provides @ref begin and @ref end to allow range-based for over pairs.
	 *
	 * @tparam mode The uniqueness mode (unique or non-unique pairs).
	 */
	template <Uniqueness mode>
	class p_range {
		// This range is meant to be non-owning.
		// NOLINTNEXTLINE(*avoid-*ref-data-members)
		/// Underlying container to get pairs from.
		Container& container;
		using container_t = std::remove_reference_t<Container>;

	public:
		/// Type of iterator returned from @ref begin and @ref end.
		using iterator = p_iterator<mode>;
		/// Type that represents the size of this range.
		using size_type = container_t::size_type;

		/**
		 * @brief Constructs a new pair range.
		 * @param c Reference to the container to iterate over.
		 */
		explicit constexpr p_range(Container& c) noexcept
			: container(c) {}

		/**
		 * @brief Returns an iterator to the first pair in the range.
		 * @return Iterator to the beginning of the pair sequence.
		 */
		constexpr iterator begin() const noexcept {
			if (std::ranges::size(container) < 2) return end();
			return iterator(container, 0, 1);
		}

		/**
		 * @brief Returns the end iterator for non-unique pairs.
		 * @return Iterator marking the end of the range.
		 */
		constexpr iterator end() const noexcept
			requires(mode == Uniqueness::non_unique)
		{
			return iterator(container, std::ranges::size(container), 0);
		}

		/**
		 * @brief Returns the end iterator for unique pairs.
		 * @return Iterator marking the end of the range.
		 */
		constexpr iterator end() const noexcept
			requires(mode == Uniqueness::unique)
		{
			const auto size = static_cast<size_type>(std::ranges::size(container));
			return iterator(container, size - 1, size);
		}

		/**
		 * @brief Returns the total number of pairs generated by this range.
		 * @return The total number of generated pairs.
		 */
		constexpr size_type size() const noexcept {
			const auto n = std::ranges::size(container);
			if constexpr (mode == Uniqueness::unique) {
				return n * (n - 1) / 2;
			} else {
				return n * (n - 1);
			}
		}
	};

public:
	/// Type alias for a range of unique pairs.
	using uniques = p_range<Uniqueness::unique>;
	/// Type alias for a range of non-unique pairs.
	using pairs = p_range<Uniqueness::non_unique>;
};
}  // namespace detail

/**
 * @brief Creates a range for iterating over all unique pairs of elements of a range.
 *
 * Generates all `(a, b)` pairs from the provided range such that each pair
 * is unique and `a` precedes `b`.
 *
 * Example:
 * @code
 * std::vector<int> values = {1, 2, 3};
 * for (auto&& [a, b] : unique_pairs(values)) {
 *     std::cout << "(" << a << ", " << b << ")\n";
 * }
 * // Output:
 * // (1, 2)
 * // (1, 3)
 * // (2, 3)
 * @endcode
 *
 * @tparam Range A random-access range type.
 * @param range The container or range to iterate over.
 * @return Range over unique pairs.
 */
template <std::ranges::random_access_range Range>
// Ranges do not get forwarded.
// NOLINTNEXTLINE(*missing-std-forward)
constexpr auto unique_pairs(Range&& range) noexcept {
	return typename detail::pairwise<Range>::uniques{range};
}

/**
 * @brief Creates a range for iterating over all non-unique pairs of elements.
 *
 * Generates all `(a, b)` pairs from the provided range such that `a != b`.
 * This includes both `(i, j)` and `(j, i)`.
 *
 * Example:
 * @code
 * std::vector<int> values = {1, 2, 3};
 * for (auto&& [a, b] : pairs(values)) {
 *     std::cout << "(" << a << ", " << b << ")\n";
 * }
 * // Output:
 * // (1, 2)
 * // (1, 3)
 * // (2, 1)
 * // (2, 3)
 * // (3, 1)
 * // (3, 2)
 * @endcode
 *
 * @tparam Range A random-access range type.
 * @param range The container or range to iterate over.
 * @return Range over non-unique pairs.
 */
template <std::ranges::random_access_range Range>
// Ranges do not get forwarded.
// NOLINTNEXTLINE(*missing-std-forward)
constexpr auto pairs(Range&& range) noexcept {
	return typename detail::pairwise<Range>::pairs{range};
}

static_assert(std::ranges::range<decltype(pairs(std::declval<std::vector<int>>()))>);
static_assert(std::ranges::range<decltype(unique_pairs(std::declval<std::vector<int>>()))>);
