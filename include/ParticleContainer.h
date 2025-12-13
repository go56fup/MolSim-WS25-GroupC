#pragma once
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <vector>

#include "CompilerTraits.h"
#include "Concepts.h"
#include "Enums.h"
#include "Particle.h"
#include "utils/MaxwellBoltzmannDistribution.h"

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
			const auto size = static_cast<size_type>(std::ranges::size(container));
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
public:
	using cell = std::vector<Particle>;
	using size_type = std::uint32_t;
	using difference_type = std::int32_t;
	using index = vec_3d<size_type>;
	using value_type = cell;

private:
	constexpr size_type linear_index(size_type i, size_type j, size_type k) const noexcept {
		auto result = (i * grid_size_.y * grid_size_.z) + (j * grid_size_.z) + k;
		assert(result < grid_size_.x * grid_size_.y * grid_size_.z && "Out of bounds cell requested");
		return result;
	}

	constexpr size_type pos_to_linear_index(const vec& pos) const noexcept {
		const auto x = static_cast<size_type>(pos.x / cutoff_radius_);
		const auto y = static_cast<size_type>(pos.y / cutoff_radius_);
		const auto z = static_cast<size_type>(pos.z / cutoff_radius_);
		return linear_index(x, y, z);
	}

	static CONSTEXPR_IF_GCC size_type div_round_up(double x, double y) noexcept {
		return static_cast<size_type>(std::ceil(x / y));
	}

	constexpr void check_if_oob(const vec& pos) const noexcept(false) {
		if (pos.x >= domain_.x || pos.y >= domain_.y || pos.z >= domain_.z) {
			throw std::domain_error(
				fmt::format(
					"Refusing to place particle at {}, which is on or outside the domain of the simulation: {}", pos,
					domain_
				)
			);
		}
	}

	std::vector<cell> grid;
	vec domain_;
	index grid_size_;
	double cutoff_radius_;

public:
	constexpr double cutoff_radius() const noexcept {
		return cutoff_radius_;
	}

	template <fwd_reference_to<vec> Vec>
	constexpr ParticleContainer(Vec domain_arg, double cutoff_radius_arg)
		: domain_{std::forward<Vec>(domain_arg)}
		, grid_size_{div_round_up(domain_.x, cutoff_radius_arg), div_round_up(domain_.y, cutoff_radius_arg), div_round_up(domain_.z, cutoff_radius_arg)}
		, cutoff_radius_{cutoff_radius_arg} {
		// Prevent implicit widening later, if x * y * z overflows for uint32, this takes care of it
		const auto grid_x = static_cast<std::size_t>(grid_size_.x);
		const auto grid_y = static_cast<std::size_t>(grid_size_.y);
		const auto grid_z = static_cast<std::size_t>(grid_size_.z);
		SPDLOG_TRACE("Setting up grid with size: {}", grid_size_);
		// TODO(tuna): see if this default constructs the particles underneath
		grid.resize(grid_x * grid_y * grid_z);
	}

	template <fwd_reference_to<Particle> ParticleT>
	constexpr void place(ParticleT&& particle) {
		check_if_oob(particle.x);
		grid[pos_to_linear_index(particle.x)].emplace_back(std::forward<ParticleT>(particle));
	}

	template <fwd_reference_to<vec> VecT, typename... Args>
		requires std::constructible_from<Particle, vec, Args...>
	constexpr void emplace(VecT&& pos, Args&&... args) {
		SPDLOG_TRACE("Emplacing particle at coords: {}", pos);
		check_if_oob(pos);
		grid[pos_to_linear_index(pos)].emplace_back(std::forward<VecT>(pos), std::forward<Args>(args)...);
	}

	// TODO(anyone): update docs
	/**
	 * @brief Create a 3D grid of particles, representing one body.
	 *
	 * @param particles Particle container to add generated particles onto
	 * @param origin Position of the lower left front-side corner
	 * @param initial_velocity Initial additional velocity of each particle
	 * @param scale Number of particles in each direction
	 * @param distance Relative distance between two particles
	 * @param mass Mass of one particle
	 * @param mean_brownian Average velocity of the Brownian Motion
	 **/
	template <std::size_t N>
	constexpr void add_cuboid(
		const vec& origin, const index& scale, double meshwidth, const vec& velocity, double mass, double brownian_mean,
		std::size_t& seq_no
	) {
		for (size_type i = 0; i < scale.x; ++i) {
			for (size_type j = 0; j < scale.y; ++j) {
				for (size_type k = 0; k < scale.z; ++k) {
					emplace(
						vec{origin.x + (i * meshwidth), origin.y + (j * meshwidth), origin.z + (k * meshwidth)},
						maxwellBoltzmannDistributedVelocity<N>(brownian_mean, seq_no) + velocity, mass, seq_no
					);
				}
			}
		}
	}

	CONSTEXPR_IF_GCC void
	add_disc(const vec& center, double radius, double meshwidth, const vec& velocity, double mass) {
		for (int i = static_cast<int>(-radius); i <= radius; i++) {
			for (int j = static_cast<int>(-radius); j <= radius; j++) {
				// TODO(tuna): switch over to manual multiplication instead of sqrt and pow
				if (std::sqrt(std::pow((i * meshwidth), 2) + std::pow((j * meshwidth), 2)) > (radius * meshwidth)) {
					continue;
				}
				emplace(vec{center.x + (i * meshwidth), center.y + (j * meshwidth), center.z}, velocity, mass);
			}
		}
	}

	// TODO(tuna): constrain return types with concepts
	constexpr auto directional_interactions() noexcept;
	constexpr auto border_cells() noexcept;
	constexpr auto enumerate_cells() noexcept;

	constexpr auto begin() const noexcept {
		return grid.begin();
	}

	constexpr auto end() const noexcept {
		return grid.end();
	}

	constexpr auto begin() noexcept {
		return grid.begin();
	}

	constexpr auto end() noexcept {
		return grid.end();
	}

	constexpr cell& operator[](const index& three_d_index) noexcept {
		return grid[linear_index(three_d_index.x, three_d_index.y, three_d_index.z)];
	}

	constexpr const cell& operator[](const index& three_d_index) const noexcept {
		return grid[linear_index(three_d_index.x, three_d_index.y, three_d_index.z)];
	}

	constexpr const cell& cell_containing(const vec& pos) const noexcept {
		return grid[pos_to_linear_index(pos)];
	}

	constexpr cell& cell_containing(const vec& pos) noexcept {
		return grid[pos_to_linear_index(pos)];
	}

	constexpr cell& operator[](size_type x, size_type y, size_type z) noexcept {
		return grid[linear_index(x, y, z)];
	}

	constexpr const cell& operator[](size_type x, size_type y, size_type z) const noexcept {
		return grid[linear_index(x, y, z)];
	}

	constexpr const vec& domain() const noexcept {
		return domain_;
	}

	constexpr const index& grid_size() const noexcept {
		return grid_size_;
	}

	constexpr range_of<Particle> auto view() noexcept {
		return grid | std::views::join;
	}

	constexpr range_of<const Particle> auto view() const noexcept {
		// TODO(tuna): i'm not sure as const is strictly necessary
		return std::as_const(grid) | std::views::join;
	}
};

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
}

namespace detail {

class interactions_iterator {
private:
	using inner = typename ParticleContainer::value_type;

public:
	/// Standard iterator category.
	using iterator_category = std::input_iterator_tag;
	/// Type used for differences between indices.
	using difference_type = ParticleContainer::difference_type;
	/// Type used for indexing.
	using size_type = ParticleContainer::size_type;
	using pointer = void;
	/// Reference to two elements of the container.
	using reference = std::pair<const ParticleContainer::index&, const ParticleContainer::index&>;
	/// The type returned from @ref operator*, which is a reference.
	using value_type = reference;

private:
	ParticleContainer* container;
	ParticleContainer::index current_cell_idx;
	ParticleContainer::index target_cell_idx;
	std::uint8_t displacement_idx;
	// TODO(tuna): maybe add 0, 0, 0 to elide first loop in calculate_forces
	static constexpr std::array<vec_3d<difference_type>, 13> displacements = {
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
	static constexpr std::size_t displacement_count = std::tuple_size_v<decltype(displacements)>;

	template <typename T>
	static constexpr std::make_signed_t<T> add_signed(T unsigned_, std::make_signed_t<T> signed_) {
		return static_cast<decltype(signed_)>(unsigned_) + signed_;
	}

	template <axis Axis>
	constexpr bool do_displacement(const vec_3d<difference_type>& displacement) noexcept {
		static constexpr auto proj_size_t = axis_ptr<size_type>(Axis);
		static constexpr auto proj_ptrdiff_t = axis_ptr<difference_type>(Axis);
		if (displacement.*proj_ptrdiff_t == 0) return true;

		// We do not need to use __builtin_add_overflow() here, as displacement can be at minimum -1.
		const auto outside_domain_bc_negative = displacement.*proj_ptrdiff_t < 0 && current_cell_idx.*proj_size_t == 0;
		if (outside_domain_bc_negative) return false;
		target_cell_idx.*proj_size_t =
			static_cast<size_type>(add_signed(current_cell_idx.*proj_size_t, displacement.*proj_ptrdiff_t));
		const auto inside_domain_max = target_cell_idx.*proj_size_t < container->grid_size().*proj_size_t;
		return inside_domain_max;
	};

public:
	struct get_end_tag {};

	constexpr interactions_iterator() noexcept
		: container(nullptr)
		, displacement_idx(0) {
		++*this;
	}

	template <fwd_reference_to<ParticleContainer::index> IndexT>
	constexpr interactions_iterator(ParticleContainer& c, IndexT&& current_cell) noexcept
		: container(&c)
		, current_cell_idx(std::forward<IndexT>(current_cell))
		, displacement_idx(0) {
		++*this;
	}

	constexpr interactions_iterator(ParticleContainer& c, get_end_tag) noexcept
		: container(&c)
		, current_cell_idx{0, 0, c.grid_size().z}
		, displacement_idx(0) {}

	/**
	 * @brief Dereferences this iterator.
	 * @return Pair of mutable references to the current elements.
	 */
	constexpr reference operator*() const noexcept {
		return reference{current_cell_idx, target_cell_idx};
	}

	/**
	 * @brief Advances the iterator to the next valid pair.
	 *
	 * @return Reference to this iterator after increment.
	 */
	constexpr interactions_iterator& operator++() noexcept {
		const auto& grid = container->grid_size();

		while (current_cell_idx.z < grid.z) {
			// Try all displacements for the current cell
			while (displacement_idx < displacement_count) {
				const auto& displacement = displacements[displacement_idx];
				// SPDLOG_TRACE("Displacing with {}, index: {}", displacement, displacement_idx);

				target_cell_idx = current_cell_idx;

				const bool success = do_displacement<axis::x>(displacement) && do_displacement<axis::y>(displacement) &&
				                     do_displacement<axis::z>(displacement);

				++displacement_idx;

				if (success) {
					return *this;
				}
			}

			// No displacement succeeded: advance to the next cell
			displacement_idx = 0;

			if (++current_cell_idx.x < grid.x) {
				target_cell_idx = current_cell_idx;
				continue;
			}
			current_cell_idx.x = 0;

			if (++current_cell_idx.y < grid.y) {
				target_cell_idx = current_cell_idx;
				continue;
			}
			current_cell_idx.y = 0;

			++current_cell_idx.z;
		}
		// This is the end iterator
		assert(
			current_cell_idx.x == 0 && current_cell_idx.y == 0 && current_cell_idx.z == grid.z && displacement_idx == 0
		);
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
		assert(container == other.container);
		/*
		SPDLOG_TRACE(
		    "Comparing: cur: {} i: {} (target: {}) ?= {} {} ({})", current_cell_idx, displacement_idx, target_cell_idx,
		    other.current_cell_idx, other.displacement_idx, other.target_cell_idx
		);
		*/
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
		return {container, ParticleContainer::index{}};
	}

	constexpr iterator end() const noexcept {
		return iterator{container, iterator::get_end_tag{}};
	}

	// TODO(tuna): add size
};

}  // namespace detail

// TODO(tuna): constrain return types with concepts
constexpr auto ParticleContainer::directional_interactions() noexcept {

	return detail::interactions_range(*this);
}

namespace detail {
class border_cell_iterator {
private:
	using inner = typename ParticleContainer::value_type;

public:
	/// Standard iterator category.
	using iterator_category = std::input_iterator_tag;
	/// Type used for differences between indices.
	using difference_type = ParticleContainer::difference_type;
	/// Type used for indexing.
	using size_type = ParticleContainer::size_type;
	using pointer = void;
	// TODO(tuna): make value_type's across the codebase const compliant
	using reference = std::pair<ParticleContainer::index, boundary_type>;
	// TODO(tuna): see if returning index behavior is better behaved than const& to index
	/// The type returned from @ref operator*.
	using value_type = reference;

private:
	ParticleContainer* container = nullptr;
	ParticleContainer::index idx;
	boundary_type type{};

	constexpr boundary_type bounds_check(const ParticleContainer::index& current) const noexcept {
		const auto& boundaries = container->grid_size();
		// TODO(tuna): find a way to elide these checks:
		// for the 1-cell thick case, the cell becomes both *_min and *_max; which only
		// happens there but causes the overhead of both checks for the majority of grid
		// configurations
		using enum boundary_type;
		boundary_type result{};
		SPDLOG_TRACE("Calculating border type of index: {}", current);
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
		SPDLOG_TRACE("got border type: {}", result);
		return result;
	};

public:
	template <fwd_reference_to<ParticleContainer::index> IndexT>
	constexpr border_cell_iterator(ParticleContainer& c, IndexT&& begin_index) noexcept
		: container(&c)
		, idx(std::forward<IndexT>(begin_index))
		, type(bounds_check(idx)) {}

	// TODO(tuna): move to defaulting default cons and memberwise default
	constexpr border_cell_iterator() noexcept = default;

	// TODO(tuna): mark return type of operator* value_type everyrwhere
	constexpr value_type operator*() const noexcept {
		return {idx, type};
	}

	constexpr border_cell_iterator& operator++() noexcept {
		const auto boundaries = container->grid_size();

		while (true) {
			if (++idx.x >= boundaries.x) {
				idx.x = 0;
				if (++idx.y >= boundaries.y) {
					idx.y = 0;
					++idx.z;
				}
			}

			type = bounds_check(idx);
			if (std::to_underlying(type) == 0) {
				continue;
			}
			return *this;
		}
	}

	/**
	 * @brief Post-increment operator.
	 * @return This iterator before increment.
	 */
	constexpr border_cell_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	/**
	 * @brief Compares two iterators for equality.
	 * @param other Iterator to compare to.
	 * @return `true` @a iff both iterators reference the same position.
	 */
	constexpr bool operator==(const border_cell_iterator& other) const noexcept {
		assert(container == other.container);
		return idx == other.idx;
	}
};

class border_cell_range {
	// NOLINTNEXTLINE(*avoid-*ref-data-members)
	ParticleContainer& container;

public:
	/// Type of iterator returned from @ref begin and @ref end.
	using iterator = border_cell_iterator;
	/// Type that represents the size of this range.
	using size_type = ParticleContainer::size_type;

	explicit constexpr border_cell_range(ParticleContainer& c) noexcept
		: container(c) {}

	constexpr iterator begin() const noexcept {
		using index = ParticleContainer::index;
		return iterator{container, index{}};
	}

	constexpr iterator end() const noexcept {

		return iterator{container, ParticleContainer::index{0, 0, container.grid_size().z}};
	}

	// TODO(tuna): add size
};

}  // namespace detail

constexpr auto ParticleContainer::border_cells() noexcept {

	return detail::border_cell_range(*this);
}

namespace detail {
class enumerate_cells_iterator {
private:
	using inner = typename ParticleContainer::value_type;

public:
	/// Standard iterator category.
	using iterator_category = std::input_iterator_tag;
	/// Type used for differences between indices.
	using difference_type = ParticleContainer::difference_type;
	/// Type used for indexing.
	using size_type = ParticleContainer::size_type;
	using pointer = void;
	// TODO(tuna): make value_type's across the codebase const compliant
	using reference = std::pair<ParticleContainer::index, ParticleContainer::cell&>;
	// TODO(tuna): see if returning index behavior is better behaved than const& to index
	/// The type returned from @ref operator*.
	using value_type = reference;

private:
	ParticleContainer* container = nullptr;
	ParticleContainer::index idx;

public:
	template <fwd_reference_to<ParticleContainer::index> IndexT>
	constexpr enumerate_cells_iterator(ParticleContainer& c, IndexT&& begin_index) noexcept
		: container(&c)
		, idx(std::forward<IndexT>(begin_index)) {}

	constexpr enumerate_cells_iterator() noexcept = default;

	constexpr value_type operator*() const noexcept {
		return {idx, (*container)[idx]};
	}

	constexpr enumerate_cells_iterator& operator++() noexcept {
		const auto boundaries = container->grid_size();

		while (true) {
			if (++idx.x >= boundaries.x) {
				idx.x = 0;
				if (++idx.y >= boundaries.y) {
					idx.y = 0;
					++idx.z;
				}
			}

			return *this;
		}
	}

	/**
	 * @brief Post-increment operator.
	 * @return This iterator before increment.
	 */
	constexpr enumerate_cells_iterator operator++(int) noexcept {
		auto tmp = *this;
		++(*this);
		return tmp;
	}

	/**
	 * @brief Compares two iterators for equality.
	 * @param other Iterator to compare to.
	 * @return `true` @a iff both iterators reference the same position.
	 */
	constexpr bool operator==(const enumerate_cells_iterator& other) const noexcept {
		assert(container == other.container);
		return idx == other.idx;
	}
};

class enumerate_cells_range {
	// NOLINTNEXTLINE(*avoid-*ref-data-members)
	ParticleContainer& container;

public:
	/// Type of iterator returned from @ref begin and @ref end.
	using iterator = enumerate_cells_iterator;
	/// Type that represents the size of this range.
	using size_type = ParticleContainer::size_type;

	explicit constexpr enumerate_cells_range(ParticleContainer& c) noexcept
		: container(c) {}

	constexpr iterator begin() const noexcept {
		using index = ParticleContainer::index;
		return iterator{container, index{}};
	}

	constexpr iterator end() const noexcept {

		return iterator{container, ParticleContainer::index{0, 0, container.grid_size().z}};
	}

	// TODO(tuna): add size
};

}  // namespace detail

constexpr auto ParticleContainer::enumerate_cells() noexcept {
	return detail::enumerate_cells_range(*this);
}

static_assert(std::ranges::range<decltype(unique_pairs(std::declval<ParticleContainer>()))>);
static_assert(std::ranges::range<decltype((std::declval<ParticleContainer>()).directional_interactions())>);
