#pragma once
#include <cstdint>
#include <ranges>
#include <vector>

namespace detail {
template <std::ranges::random_access_range Container>
class pairwise {
private:
	enum class Uniqueness : std::uint8_t { unique, non_unique };

	template <Uniqueness mode>
	class p_iterator {
	private:
		using container_t = std::remove_reference_t<Container>;
		using inner = typename container_t::value_type;

	public:
		using iterator_category = std::output_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using size_type = typename container_t::size_type;
		using pointer = void;
		using reference = std::pair<inner&, inner&>;
		using const_reference = std::pair<const inner&, const inner&>;
		using value_type = reference;

	private:
		container_t& container;
		size_type outer_idx;
		size_type inner_idx;

	public:
		constexpr p_iterator(Container& c, size_type i_start, size_type j_start) noexcept
			: container(c)
			, outer_idx(i_start)
			, inner_idx(j_start) {}

		constexpr const_reference operator*() const noexcept {
			return {container[outer_idx], container[inner_idx]};
		}

		constexpr reference operator*() noexcept {
			return {container[outer_idx], container[inner_idx]};
		}

		constexpr p_iterator& operator++() noexcept {
			const auto size = std::ranges::size(container);

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

		constexpr p_iterator operator++(int) noexcept {
			auto tmp = *this;
			++(*this);
			return tmp;
		}

		constexpr bool operator==(const p_iterator& other) const noexcept {
			return &container == &other.container && outer_idx == other.outer_idx && inner_idx == other.inner_idx;
		}
	};

	template <Uniqueness mode>
	class p_range {
		Container& container;
		using container_t = std::remove_reference_t<Container>;

	public:
		using iterator = p_iterator<mode>;
		using size_type = typename container_t::size_type;

		constexpr p_range(Container& c) noexcept
			: container(c) {}

		constexpr iterator begin() noexcept {
			if (std::ranges::size(container) < 2) return end();
			return iterator(container, 0, 1);
		}

		constexpr iterator end() noexcept
			requires(mode == Uniqueness::non_unique)
		{
			return iterator(container, std::ranges::size(container), 0);
		}

		constexpr iterator end() noexcept
			requires(mode == Uniqueness::unique)
		{
			const auto size = std::ranges::size(container);
			return iterator(container, size - 1, size);
		}

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
	using uniques = p_range<Uniqueness::unique>;
	using pairs = p_range<Uniqueness::non_unique>;
};
}  // namespace detail

template <std::ranges::range Range>
constexpr auto unique_pairs(Range&& range) noexcept {
	return typename detail::pairwise<Range>::uniques{range};
}

template <std::ranges::range Range>
constexpr auto pairs(Range&& range) noexcept {
	return typename detail::pairwise<Range>::pairs{range};
}

// Rather than having to forward through every single member function of std::vector and then being bound to it for
// underlying storage, having a (unique_)pairs adaptor requires no passthrough, is generic and still easy to use:
// for (auto&& [a, b]: /*unique_*/pairs(particles))
// for (auto& p: particles)
// for (auto&& [a, b]: pairs(std::vector{1, 2, 3})
// NOTE: The last one requires C++23 in order not to dangle.

// The only problem with a distinct ParticleContainer isn't just having to write a bunch of
// constexpr auto method(params p) { return container.method(p) }
// but also, other misc. compatibility issues like supporting std::array's inplace construction. Using aggregate
// initialization like std::array is impossible because of the user defined generic forwarding constructor. So a
// solution like std::pair's std::piecewise_construct is necessary, which is just plain painful:
// https://godbolt.org/z/4GfM9Ge38

using ParticleContainer = std::vector<Particle>;
