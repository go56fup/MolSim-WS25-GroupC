#include "gtest_constexpr/macros.hpp"
#include <gtest/gtest.h>

#include "testing_utilities.hpp"

#include "grid/particle_container/particle_container.hpp"

// Check that the border cells for a 2x2x2 grid are recognized correctly
TEST(BorderCellTests, TwoByTwoByTwo) {
	static constexpr std::size_t border_cell_count = 2 * 2 * 2;
	static constexpr std::array<particle_container::index, border_cell_count> expected{
		{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}
	};
	GTEST_CXP auto indices = CONSTEXPR_TWO_STEP<expected.size(), particle_container::index>([&] {
		particle_container container(vec{2, 2, 2}, 1);
		std::vector<particle_container::index> surrogate;

		for (const auto& [idx, _] : container.border_cells()) {
			surrogate.push_back(idx);
		}
		return surrogate;
	});
	GTEST_CXP bool ranges_eq = std::ranges::equal(indices, expected);
	STATIC_EXPECT_TRUE(ranges_eq);
}

// Check that the border cells for a 3x3x3 grid are recognized correctly
TEST(BorderCellTests, ThreeByThreeByThree) {
	static constexpr std::size_t border_cell_count = (3 * 3 * 3) - 1;
	static constexpr std::array<particle_container::index, border_cell_count> expected{{
		{0, 0, 0}, {1, 0, 0}, {2, 0, 0},

		{0, 1, 0}, {1, 1, 0}, {2, 1, 0},

		{0, 2, 0}, {1, 2, 0}, {2, 2, 0},

		{0, 0, 1}, {1, 0, 1}, {2, 0, 1},

		{0, 1, 1}, {2, 1, 1},

		{0, 2, 1}, {1, 2, 1}, {2, 2, 1},

		{0, 0, 2}, {1, 0, 2}, {2, 0, 2},

		{0, 1, 2}, {1, 1, 2}, {2, 1, 2},

		{0, 2, 2}, {1, 2, 2}, {2, 2, 2},
	}};
	GTEST_CXP auto indices = CONSTEXPR_TWO_STEP<expected.size(), particle_container::index>([&] {
		particle_container container(vec{3, 3, 3}, 1);
		std::vector<particle_container::index> surrogate;

		for (const auto& [idx, _] : container.border_cells()) {
			surrogate.push_back(idx);
		}
		return surrogate;
	});
	GTEST_CXP bool ranges_eq = std::ranges::equal(indices, expected);
	STATIC_EXPECT_TRUE(ranges_eq);
}
