#include <algorithm>
#include <array>
#include <utility>

#include <gtest/gtest.h>
#include "gtest_constexpr/macros.hpp"
#include "testing_utilities.hpp"

#include "iterators/pairwise.hpp"

// Check that unique_pairs results in the expected pairs
TEST(UniquePairTests, HappyPath) {
	GTEST_CXP std::array values = {1, 2, 3, 4, 5};
	GTEST_CXP std::array<std::pair<int, int>, 10> expected{
		{{1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}, {3, 4}, {3, 5}, {4, 5}}
	};

	GTEST_CXP bool ranges_eq = std::ranges::equal(unique_pairs(values), expected);
	STATIC_EXPECT_TRUE(ranges_eq);
}

// Check that pairs results in the expected pairs
TEST(PairTests, HappyPath) {
	GTEST_CXP std::array values = {1, 2, 3, 4, 5};
	// clang-format off
	GTEST_CXP std::array<std::pair<int, int>, 20> expected{{
			{1, 2}, {1, 3}, {1, 4}, {1, 5},
			{2, 1}, {2, 3}, {2, 4}, {2, 5},
			{3, 1}, {3, 2}, {3, 4}, {3, 5},
			{4, 1}, {4, 2}, {4, 3}, {4, 5},
			{5, 1}, {5, 2}, {5, 3}, {5, 4}
	}};
	// clang-format on
	GTEST_CXP bool ranges_eq = std::ranges::equal(pairs(values), expected);
	STATIC_EXPECT_TRUE(ranges_eq);
}

// Check that no pairs are outputted if there aren't enough elements with unique_pairs
TEST(UniquePairTests, NoPairs) {
	GTEST_CXP std::array<int, 1> values{1};
	GTEST_CXP bool unique_range_empty = std::ranges::empty(unique_pairs(values));
	STATIC_EXPECT_TRUE(unique_range_empty);
}

// Check that no pairs are outputted if there aren't enough elements with pairs
TEST(PairTests, NoPairs) {
	GTEST_CXP std::array<int, 1> values{1};
	GTEST_CXP bool nonunique_range_empty = std::ranges::empty(pairs(values));
	STATIC_EXPECT_TRUE(nonunique_range_empty);
}
