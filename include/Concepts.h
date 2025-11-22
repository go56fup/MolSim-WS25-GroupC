#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

template <typename Range, typename Element>
concept range_of = std::ranges::range<Range> && (std::same_as<std::ranges::range_value_t<Range>, Element> ||
                                                 std::same_as<std::ranges::range_reference_t<Range>, Element>);
template <typename Candidate, typename Expected>
// same with the underlying type = rvalue bind case, same with lvalue ref = lvalue bind case
concept fwd_reference_to = std::same_as<std::remove_cvref_t<Candidate>, Expected>;
