#include "ParticleContainer.h"
#include <gtest/gtest.h>
#include <gtest_constexpr.h>

// TODO(tuna): convert to use expect_all somehow so that it can also be checked at compile time
TEST(DirectionalInteraction, BasicTest) {
	ParticleContainer container(2, 2, 2, 1);
	// Number of edges in 2^3 = K8.
	static constexpr std::size_t interactions = 8 * 7 / 2;
	[[maybe_unused]] GTEST_CXP std::array<std::pair<ParticleContainer::index, ParticleContainer::index>, interactions>
		expected{{{{0, 0, 0}, {0, 0, 1}}, {{0, 0, 0}, {0, 1, 0}}, {{0, 0, 0}, {0, 1, 1}}, {{0, 0, 0}, {1, 0, 0}},
	              {{0, 0, 0}, {1, 0, 1}}, {{0, 0, 0}, {1, 1, 0}}, {{0, 0, 0}, {1, 1, 1}}, {{1, 0, 0}, {1, 0, 1}},
	              {{1, 0, 0}, {1, 1, 0}}, {{1, 0, 0}, {1, 1, 1}}, {{0, 1, 0}, {0, 1, 1}}, {{0, 1, 0}, {1, 0, 0}},
	              {{0, 1, 0}, {1, 0, 1}}, {{0, 1, 0}, {1, 1, 0}}, {{0, 1, 0}, {1, 1, 1}}, {{1, 1, 0}, {1, 1, 1}},
	              {{0, 0, 1}, {0, 1, 0}}, {{0, 0, 1}, {0, 1, 1}}, {{0, 0, 1}, {1, 0, 0}}, {{0, 0, 1}, {1, 0, 1}},
	              {{0, 0, 1}, {1, 1, 0}}, {{0, 0, 1}, {1, 1, 1}}, {{1, 0, 1}, {1, 1, 0}}, {{1, 0, 1}, {1, 1, 1}},
	              {{0, 1, 1}, {1, 0, 0}}, {{0, 1, 1}, {1, 0, 1}}, {{0, 1, 1}, {1, 1, 0}}, {{0, 1, 1}, {1, 1, 1}}}};
	std::size_t i = 0;
	for (const auto& pair : container.directional_interactions()) {
		EXPECT_EQ(pair, expected[i++]);
	}
	EXPECT_EQ(i, interactions);
}

// TODO(tuna): add tests for boxes of weird shapes like 2*3*1, 1*1*1 etc.
