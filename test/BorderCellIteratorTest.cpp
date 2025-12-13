#include <gtest/gtest.h>
#include <gtest_constexpr.h>

#include "ParticleContainer.h"

// TODO(tuna): use the constexpr two-step to extract the indices out and fold over them to compare
TEST(BorderCellTests, TwoByTwoByTwo) {
	static constexpr std::size_t border_cell_count = 2 * 2 * 2;
	static constexpr std::array<ParticleContainer::index, border_cell_count> expected{
		{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}
	};
	ParticleContainer container(vec{2, 2, 2}, 1);
	std::size_t i = 0;
	for (const auto& [idx, _] : container.border_cells()) {
		EXPECT_EQ(idx, expected[i++]);
	}
	EXPECT_EQ(i, border_cell_count);
}

TEST(BorderCellTests, ThreeByThreeByThree) {
	static constexpr std::size_t border_cell_count = (3 * 3 * 3) - 1;
	static constexpr std::array<ParticleContainer::index, border_cell_count> expected{{
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

	ParticleContainer container(vec{3, 3, 3}, 1);
	std::size_t i = 0;
	for (const auto& [idx, _] : container.border_cells()) {
		EXPECT_EQ(idx, expected[i++]);
	}
	EXPECT_EQ(i, border_cell_count);
}
