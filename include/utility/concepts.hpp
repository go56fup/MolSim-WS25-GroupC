#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

namespace detail {
// libc++ does not support this concept as of now (2025-11-29), so we copy the implementation from libstdc++.
#ifndef __cpp_lib_ranges_as_const
// https://github.com/gcc-mirror/gcc/blob/55049da531f5246499ddd2bd882928b57ac97519/libstdc%2B%2B-v3/include/bits/stl_iterator.h#L2616
template <std::indirectly_readable It>
using iter_const_reference_t = std::common_reference_t<const std::iter_value_t<It>&&, std::iter_reference_t<It>>;

// https://github.com/gcc-mirror/gcc/blob/55049da531f5246499ddd2bd882928b57ac97519/libstdc%2B%2B-v3/include/bits/stl_iterator.h#L2624
template <typename It>
concept constant_iterator =
	std::input_iterator<It> && std::same_as<iter_const_reference_t<It>, std::iter_reference_t<It>>;

// https://github.com/gcc-mirror/gcc/blob/55049da531f5246499ddd2bd882928b57ac97519/libstdc%2B%2B-v3/include/bits/ranges_base.h#L633
template <typename Tp>
concept constant_range = std::ranges::input_range<Tp> && constant_iterator<std::ranges::iterator_t<Tp>>;
#else
template <typename T>
concept constant_range = std::ranges::constant_range<T>;
#endif
}  // namespace detail

template <typename Range, typename Element>
concept range_of =
	std::ranges::range<Range> && std::same_as<std::ranges::range_value_t<Range>, std::remove_const_t<Element>> &&
	// if the range is const, it cannot be bound to non-const Element, but not vice versa
	(detail::constant_range<Range> ? std::is_const_v<Element> : true);

template <typename Candidate, typename Expected>
// same with the underlying type = rvalue bind case, same with lvalue ref = lvalue bind case
concept fwd_reference_to = std::same_as<std::remove_cvref_t<Candidate>, Expected>;

// TODO(tuna): actually differentiate between forwarding references and plain unqualification
template <typename Candidate, typename Referenced>
concept same_unqualified_as = std::same_as<std::remove_cvref_t<Candidate>, Referenced>;

