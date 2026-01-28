#pragma once

#include <string_view>

#include <algorithm>
#include <span>

#include <fmt/compile.h>
#include <spdlog/spdlog.h>

#include "utility/fixed_string.hpp"

template <std::size_t MaxSize, typename ValueT, typename Builder>
consteval std::span<const ValueT> constexpr_two_step(Builder&& builder) {
	struct calculated_data {
		std::array<ValueT, MaxSize> storage{};
		std::size_t size{};
	};

	constexpr auto data = std::invoke([&] {
		const auto dynamic = std::invoke(std::forward<Builder>(builder));
		calculated_data oversized{};
		auto current = oversized.storage.begin();
		for (const auto& element : dynamic) {
			*current++ = element;
			// TODO(tuna): add support for ranges or remove
			/*			if constexpr (std::ranges::range<decltype(element)>) {
			current = std::ranges::copy(element, current).out;
		} else {
			*current++ = element;
			}*/
		}
		oversized.size =
			static_cast<std::size_t>(std::distance(oversized.storage.begin(), current));
		return oversized;
	});

	static constexpr auto static_data = std::invoke([&] {
		std::array<ValueT, data.size> result{};
		std::ranges::copy(data.storage | std::views::take(data.size), result.begin());
		return result;
	});
	return static_data;
}
