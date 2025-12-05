#include <gtest/gtest.h>
#include <gtest_constexpr.h>

#include "ParticleContainer.h"

// TODO(tuna): flesh out tests

// NOLINTBEGIN(*magic-numbers)
TEST(ParticleContainer, Cuboid2DPlacement) {
	ParticleContainer container(10, 10, 10, 1);
	std::size_t seq_no = 0;
	container.add_cuboid<2>({5, 5, 5}, {2, 2, 1}, 1, {}, 100, 0.1, seq_no);
	EXPECT_EQ(container.cell_containing({5, 5, 5}).at(0).m, 100);
	EXPECT_EQ(container.cell_containing({6, 6, 5}).at(0).m, 100);
	EXPECT_EQ(container.cell_containing({6, 5, 5}).at(0).m, 100);
}

TEST(ParticleContainer, Cuboid3DPlacement) {
	ParticleContainer container(10, 10, 10, 1);
	std::size_t seq_no = 0;
	container.add_cuboid<3>({5, 5, 5}, {2, 2, 2}, 1, {}, 100, 0.1, seq_no);
	EXPECT_EQ(container.cell_containing({5, 6, 5}).at(0).m, 100);
	EXPECT_EQ(container.cell_containing({6, 6, 6}).at(0).m, 100);
	EXPECT_EQ(container.cell_containing({6, 5, 6}).at(0).m, 100);
}

// NOLINTEND(*magic-numbers)
