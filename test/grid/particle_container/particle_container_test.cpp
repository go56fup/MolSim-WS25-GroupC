#include <algorithm>

#include "gtest_constexpr/macros.hpp"
#include <gtest/gtest.h>

#include "grid/particle_container/particle_container.hpp"

#include "testing_utilities.hpp"

// TODO(tuna): flesh out tests

// NOLINTBEGIN(*magic-numbers)
// Check that a square is placed into the correct position in a domain
TEST(particle_container, SquarePlacement) {
	GTEST_CXP_GCC auto ok = std::invoke([] {
		particle_container container(vec{10, 10, 10}, 1);
		std::size_t seq_no = 0;
		container.add_cuboid<2>({5, 5, 5}, {2, 2, 1}, 1, {}, 100, 0.1, seq_no);
		return std::array{
			container.cell_containing({5, 5, 5}).at(0).m == 100,
			container.cell_containing({6, 6, 5}).at(0).m == 100,
			container.cell_containing({6, 5, 5}).at(0).m == 100
		};
	});
	GCC_STATIC_EXPECT_ALL(ok);
}

// Check that a cuboid is placed into the correct position in a domain
TEST(particle_container, Cuboid3DPlacement) {
	GTEST_CXP_GCC auto ok = std::invoke([] {
		particle_container container(vec{10, 10, 10}, 1);
		std::size_t seq_no = 0;
		container.add_cuboid<3>({5, 5, 5}, {2, 2, 2}, 1, {}, 100, 0.1, seq_no);
		return std::array{
			container.cell_containing({5, 6, 5}).at(0).m == 100,
			container.cell_containing({6, 6, 6}).at(0).m == 100,
			container.cell_containing({6, 5, 6}).at(0).m == 100
		};
	});
	STATIC_EXPECT_ALL(ok);
}

// TODO(tuna): this test doesn't really make sense, either put in a special debug version of
// add_cuboid that has this behavior or just remove entirely
#if 0
TEST(particle_container, ParticleIdentifiers) {
	particle_container container(vec{100, 100, 100}, 1);
	std::size_t seq_no = 0;
	static constexpr auto particle_count = 50 * 50;
	container.add_cuboid<2>({2, 2, 2}, {50, 50, 1}, 1, {}, 100, 0.1, seq_no);
	std::vector<decltype(particle::type)> seen_types;
	for (auto& p : container.particles()) {
		const auto type_repeats = std::ranges::contains(seen_types, p.type);
		EXPECT_FALSE(type_repeats);
		seen_types.push_back(p.type);
	}
	EXPECT_EQ(seen_types.size(), particle_count);
}
#endif

// NOLINTEND(*magic-numbers)
