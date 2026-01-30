#include "gtest_constexpr/macros.hpp"
#include "testing_utilities.hpp"
#include <gtest/gtest.h>

#include <array>

#include "grid/particle_container/particle_container.hpp"
#include "iterators/neighbors.hpp"

// Check that all neighboring cells of a corner cell in a 2x2x2 grid are recognized correctly
TEST(Neighbors, TwoByTwoByTwoCorner) {
	using index = particle_container::index;

	static constexpr std::size_t neighbor_count = 7;

	GTEST_CXP
	std::array<index, neighbor_count> expected{{
		{0, 0, 1},
		{0, 1, 0},
		{0, 1, 1},
		{1, 0, 0},
		{1, 0, 1},
		{1, 1, 0},
		{1, 1, 1},
	}};

	GTEST_CXP auto result = CONSTEXPR_TWO_STEP<neighbor_count, index>([&] {
		particle_container container(vec{2, 2, 2}, 1);
		std::vector<index> surrogate;

		const index center{0, 0, 0};
		for (const auto& neighbor : neighbors_range(container, center)) {
			surrogate.push_back(neighbor);
		}
		return surrogate;
	});

	for (const auto& i: result) {
		SPDLOG_INFO("got: {}", i);
	}
	GTEST_CXP bool ranges_eq = std::ranges::equal(result, expected);
	STATIC_EXPECT_TRUE(ranges_eq);
}

// Check that the innermost cell of a 3x3x3 grid has 26 neigbors
TEST(Neighbors, InteriorCellHas26Neighbors) {
	GTEST_CXP auto count = std::invoke([&] {
		particle_container container(vec{3, 3, 3}, 1);
		particle_container::index center{1, 1, 1};
		std::size_t count = 0;
		for ([[maybe_unused]] const auto& n : neighbors_range(container, center)) {
			++count;
		}
		return count;
	});

	STATIC_EXPECT_EQ(count, 26);
}
