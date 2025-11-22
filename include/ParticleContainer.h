#pragma once
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "Particle.h"

// TOOD(tuna): split out into files

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
		using container_t = std::remove_reference_t<Container>;
		using inner_wo_quals = typename container_t::value_type;
		using inner = std::conditional_t<std::is_const_v<container_t>, const inner_wo_quals, inner_wo_quals>;

	public:
		/// Standard iterator category.
		using iterator_category = std::input_iterator_tag;
		/// Type used for differences between indices.
		using difference_type = std::ptrdiff_t;
		/// Type used for indexing.
		using size_type = typename container_t::size_type;
		using pointer = void;
		/// Reference to two elements of the container.
		using reference = std::pair<inner&, inner&>;
		/// The type returned from @ref operator*, which is a reference.
		using value_type = reference;

	private:
		container_t* container;
		size_type outer_idx;
		size_type inner_idx;

	public:
		constexpr p_iterator() noexcept
			: container(nullptr)
			, outer_idx(0)
			, inner_idx(0) {}

		/**
		 * @brief Constructs a new pairwise iterator.
		 *
		 * @param c Reference to the underlying container.
		 * @param i_start Initial outer index.
		 * @param j_start Initial inner index.
		 */
		constexpr p_iterator(Container& c, size_type i_start, size_type j_start) noexcept
			: container(&c)
			, outer_idx(i_start)
			, inner_idx(j_start) {}

		/**
		 * @brief Dereferences this iterator.
		 * @return Pair of mutable references to the current elements.
		 */

		constexpr reference operator*() const noexcept {
			return reference{(*container)[outer_idx], (*container)[inner_idx]};
		}

		/**
		 * @brief Advances the iterator to the next valid pair.
		 *
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
		 * @return Reference to this iterator before increment.
		 */
		constexpr p_iterator& operator++(int) noexcept {
			auto tmp = *this;
			++(*this);
			return tmp;
		}

		/**
		 * @brief Compares two iterators for equality.
		 * @param other Iterator to compare to.
		 * @return `true` @a iff both iterators reference the same position in the same container.
		 */
		constexpr bool operator==(const p_iterator& other) const noexcept {
			return container == other.container && outer_idx == other.outer_idx && inner_idx == other.inner_idx;
		}
	};

	/**
	 * @brief Range that yields pairs of elements based on the uniqueness mode.
	 *
	 * Provides @ref begin and @ref end to allow range-based for over pairs.
	 *
	 * @tparam mode The uniqueness mode (unique or non-unique).
	 */
	template <Uniqueness mode>
	class p_range {
		// NOLINTNEXTLINE(*avoid-*ref-data-members)
		Container& container;
		using container_t = std::remove_reference_t<Container>;

	public:
		/// Type of iterator returned from @ref begin and @ref end.
		using iterator = p_iterator<mode>;
		/// Type that represents the size of this range.
		using size_type = typename container_t::size_type;

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
			const auto size = std::ranges::size(container);
			return iterator(container, size - 1, size);
		}

		/**
		 * @brief Returns the total number of pairs generated by this range.
		 *
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
 * @brief Creates a range for iterating over all unique pairs of elements.
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
template <std::ranges::range Range>
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
template <std::ranges::range Range>
// NOLINTNEXTLINE(*missing-std-forward)
constexpr auto pairs(Range&& range) noexcept {
	return typename detail::pairwise<Range>::pairs{range};
}

/**
 * @brief The container type used to hold particles.
 *
 * This type is used throughout the simulation to store all particle objects.
 */
class ParticleContainer {
private:
	constexpr std::size_t linear_index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
		return (i * domain_.x * domain_.y) + (j * domain_.y) + k;
	}

public:
	using cell = std::vector<Particle>;
	static constexpr double cutoff_radius = 3.0;

	template <fwd_reference_to<index_3d> DomainT>
	constexpr explicit ParticleContainer(DomainT&& domain)
		: domain_(std::forward<DomainT>(domain)) {
		const auto x = static_cast<double>(domain_.x);
		const auto y = static_cast<double>(domain_.y);
		const auto z = static_cast<double>(domain_.z);
		grid.reserve(
			static_cast<std::size_t>(
				div_round_up(x, cutoff_radius) * div_round_up(y, cutoff_radius) * div_round_up(z, cutoff_radius)
			)
		);
	}

	constexpr explicit ParticleContainer(std::size_t x, std::size_t y, std::size_t z)
		: domain_{x, y, z} {
		const auto x_ = static_cast<double>(domain_.x);
		const auto y_ = static_cast<double>(domain_.y);
		const auto z_ = static_cast<double>(domain_.z);
		grid.reserve(
			static_cast<std::size_t>(
				div_round_up(x_, cutoff_radius) * div_round_up(y_, cutoff_radius) * div_round_up(z_, cutoff_radius)
			)
		);
	}

	template <fwd_reference_to<Particle> ParticleT>
	constexpr void place(ParticleT&& particle) noexcept {
		grid[pos_to_linear_index(particle.x)] = std::forward<ParticleT>(particle);
	}

	template <fwd_reference_to<vec> VecT, typename... Args>
		requires std::constructible_from<Particle, vec, Args...>
	constexpr void emplace(VecT&& pos, Args&&... args) {
		grid[pos_to_linear_index(pos)].emplace_back(std::forward<VecT>(pos), std::forward<Args>(args)...);
	}

	// TODO(tuna): constrain return types with concepts
	constexpr auto directional_interactions() noexcept;

	constexpr auto begin() const noexcept {
		return grid.begin();
	}

	constexpr auto end() const noexcept {
		return grid.end();
	}

	constexpr cell& operator[](const index_3d& three_d_index) noexcept {
		return grid[linear_index(three_d_index.x, three_d_index.y, three_d_index.z)];
	}

	constexpr const cell& operator[](const index_3d& three_d_index) const noexcept {
		return grid[linear_index(three_d_index.x, three_d_index.y, three_d_index.z)];
	}

	constexpr const index_3d& domain() const noexcept {
		return domain_;
	}

	constexpr range_of<Particle> auto view() noexcept {
		return grid | std::views::join;
	}

	// TODO(tuna): fix this
	constexpr range_of</* const */ Particle> auto view() const noexcept {
		return std::as_const(grid) | std::views::join;
	}

private:
	static constexpr double div_round_up(double x, double y) noexcept {
		return (x + y - 1) / y;
	}

	constexpr std::size_t pos_to_linear_index(const vec& pos) {
		const auto x = static_cast<std::size_t>(pos.x / cutoff_radius);
		const auto y = static_cast<std::size_t>(pos.y / cutoff_radius);
		const auto z = static_cast<std::size_t>(pos.z / cutoff_radius);
		return linear_index(x, y, z);
	}

	std::vector<cell> grid;
	index_3d domain_;

public:
	using size_type = std::size_t;
	using value_type = cell;
};

namespace detail {

class interactions_iterator {
private:
	using inner = typename ParticleContainer::value_type;

public:
	/// Standard iterator category.
	using iterator_category = std::input_iterator_tag;
	/// Type used for differences between indices.
	using difference_type = std::ptrdiff_t;
	/// Type used for indexing.
	using size_type = std::size_t;
	using pointer = void;
	/// Reference to two elements of the container.
	using reference = std::pair<const index_3d&, const index_3d&>;
	/// The type returned from @ref operator*, which is a reference.
	using value_type = reference;

private:
	ParticleContainer* container;
	index_3d current_cell_idx;
	std::size_t displacement_idx;

	enum class axis : std::uint8_t { x, y, z };

	template <typename T>
	static consteval auto axis_ptr(axis a) {
		switch (a) {
		case axis::x:
			return &vec_3d<T>::x;
		case axis::y:
			return &vec_3d<T>::y;
		case axis::z:
			return &vec_3d<T>::z;
		}
		std::unreachable();
	};

	template <typename T>
	static constexpr std::make_signed_t<T> add_signed(T unsigned_, std::make_signed_t<T> signed_) {
		return static_cast<decltype(signed_)>(unsigned_) + signed_;
	}

	template <axis a>
	constexpr bool out_of_bounds(const vec_3d<std::ptrdiff_t>& displacement) const noexcept {
		static constexpr auto proj_size_t = axis_ptr<std::size_t>(a);
		static constexpr auto proj_ptrdiff_t = axis_ptr<std::ptrdiff_t>(a);

		// TODO(tuna): see if caching the result of the below add_signed (which is just do_displacement)
		// is worth the increase in the iterator's size in terms of performance
		return (displacement.*proj_ptrdiff_t < 0 && current_cell_idx.*proj_size_t == 0) ||
		       (add_signed(current_cell_idx.*proj_size_t, displacement.*proj_ptrdiff_t) >
		        static_cast<std::ptrdiff_t>(container->domain().*proj_size_t));
	};

	template <axis a>
	constexpr std::size_t do_displacement(const vec_3d<std::ptrdiff_t>& displacement) const noexcept {
		static constexpr auto proj_size_t = axis_ptr<std::size_t>(a);
		static constexpr auto proj_ptrdiff_t = axis_ptr<std::ptrdiff_t>(a);
		return static_cast<std::size_t>(add_signed(current_cell_idx.*proj_size_t, displacement.*proj_ptrdiff_t));
	}

public:
	constexpr interactions_iterator() noexcept
		: container(nullptr)
		, displacement_idx(0) {}

	template <fwd_reference_to<index_3d> IndexT>
	constexpr interactions_iterator(ParticleContainer& c, IndexT&& current_cell) noexcept
		: container(&c)
		, current_cell_idx(std::forward<IndexT>(current_cell))
		, displacement_idx(0) {}

	/**
	 * @brief Dereferences this iterator.
	 * @return Pair of mutable references to the current elements.
	 */
	static constexpr std::array<vec_3d<std::ptrdiff_t>, 13> displacements = {
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

	constexpr reference operator*() const noexcept {
		const auto& displacement = displacements[displacement_idx];
		return reference{
			current_cell_idx,
			{do_displacement<axis::x>(displacement), do_displacement<axis::y>(displacement),
		     do_displacement<axis::z>(displacement)}
		};
	}

	/**
	 * @brief Advances the iterator to the next valid pair.
	 *
	 * @return Reference to this iterator after increment.
	 */
	constexpr interactions_iterator& operator++() noexcept {
		for (; displacement_idx < displacements.size(); ++displacement_idx) {
			const auto& displacement = displacements[displacement_idx];
			if (!out_of_bounds<axis::x>(displacement) && !out_of_bounds<axis::y>(displacement) &&
			    !out_of_bounds<axis::z>(displacement)) {
				break;
			}
		}
		return *this;
	}

	/**
	 * @brief Post-increment operator.
	 * @return This iterator before increment.
	 */
	constexpr interactions_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	/**
	 * @brief Compares two iterators for equality.
	 * @param other Iterator to compare to.
	 * @return `true` @a iff both iterators reference the same position.
	 */
	constexpr bool operator==(const interactions_iterator& other) const noexcept {
		// TODO(tuna): make sure there only exists one grid so this check makes sense
		return current_cell_idx == other.current_cell_idx && displacement_idx == other.displacement_idx;
	}
};

class interactions_range {
	// NOLINTNEXTLINE(*avoid-*ref-data-members)
	ParticleContainer& container;

public:
	/// Type of iterator returned from @ref begin and @ref end.
	using iterator = interactions_iterator;
	/// Type that represents the size of this range.
	using size_type = std::size_t;

	explicit constexpr interactions_range(ParticleContainer& c) noexcept
		: container(c) {}

	constexpr iterator begin() const noexcept {
		return {container, index_3d{}};
	}

	constexpr iterator end() const noexcept {
		return {container, container.domain()};
	}

	// TODO(tuna): add size
};

}  // namespace detail

// TODO(tuna): constrain return types with concepts
constexpr auto ParticleContainer::directional_interactions() noexcept {
	return detail::interactions_range(*this);
}

static_assert(std::ranges::range<decltype(unique_pairs(std::declval<ParticleContainer>()))>);
static_assert(std::ranges::range<decltype((std::declval<ParticleContainer>()).directional_interactions())>);
