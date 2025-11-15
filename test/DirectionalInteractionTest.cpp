#include "ParticleContainer.h"
#include <gtest/gtest.h>
#include <gtest_constexpr.h>

TEST(DirectionalInteraction, BasicTest) {
	ParticleContainer container(2, 2, 2, 1);
	for (const auto& [current_cell_idx, target_cell_idx] : container.directional_interactions()) {
		SPDLOG_TRACE("Got: {} -> {}", current_cell_idx, target_cell_idx);
	}
}
