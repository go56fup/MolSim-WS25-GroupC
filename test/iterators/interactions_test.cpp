#include "gtest_constexpr/macros.hpp"
#include "testing_utilities.hpp"
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

#include "grid/particle_container/particle_container.hpp"

// Check that the directionally interacting cells for a 2x2x2 grid are recognized correctly
TEST(DirectionalInteraction, TwoByTwoByTwo) {
	// Number of edges in 2^3 = K8.
	static constexpr std::size_t interaction_count = 8 * 7 / 2;
	using interaction = std::pair<particle_container::index, particle_container::index>;
	GTEST_CXP
	std::array<interaction, interaction_count> expected{
		{{{0, 0, 0}, {0, 0, 1}}, {{0, 0, 0}, {0, 1, 0}}, {{0, 0, 0}, {0, 1, 1}},
	     {{0, 0, 0}, {1, 0, 0}}, {{0, 0, 0}, {1, 0, 1}}, {{0, 0, 0}, {1, 1, 0}},
	     {{0, 0, 0}, {1, 1, 1}}, {{1, 0, 0}, {1, 0, 1}}, {{1, 0, 0}, {1, 1, 0}},
	     {{1, 0, 0}, {1, 1, 1}}, {{0, 1, 0}, {0, 1, 1}}, {{0, 1, 0}, {1, 0, 0}},
	     {{0, 1, 0}, {1, 0, 1}}, {{0, 1, 0}, {1, 1, 0}}, {{0, 1, 0}, {1, 1, 1}},
	     {{1, 1, 0}, {1, 1, 1}}, {{0, 0, 1}, {0, 1, 0}}, {{0, 0, 1}, {0, 1, 1}},
	     {{0, 0, 1}, {1, 0, 0}}, {{0, 0, 1}, {1, 0, 1}}, {{0, 0, 1}, {1, 1, 0}},
	     {{0, 0, 1}, {1, 1, 1}}, {{1, 0, 1}, {1, 1, 0}}, {{1, 0, 1}, {1, 1, 1}},
	     {{0, 1, 1}, {1, 0, 0}}, {{0, 1, 1}, {1, 0, 1}}, {{0, 1, 1}, {1, 1, 0}},
	     {{0, 1, 1}, {1, 1, 1}}}
	};

	GTEST_CXP auto result = CONSTEXPR_TWO_STEP<interaction_count, interaction>([&] {
		particle_container container(vec{2, 2, 2}, 1);
		std::vector<interaction> surrogate;

		for (const auto& pair : container.directional_interactions()) {
			surrogate.push_back(pair);
		}
		return surrogate;
	});
	GTEST_CXP bool ranges_eq = std::ranges::equal(result, expected);
	STATIC_EXPECT_TRUE(ranges_eq);
}

// TODO(tuna): add tests for boxes of weird shapes like 2*3*1, 1*1*1 etc.
