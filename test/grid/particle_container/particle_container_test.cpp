#include <algorithm>

#include "gtest_constexpr/macros.hpp"
#include <gtest/gtest.h>

#include "grid/particle_container/particle_container.hpp"

#include "simulation/config/entities.hpp"
#include "testing_utilities.hpp"

// TODO(tuna): flesh out tests

// NOLINTBEGIN(*magic-numbers)
// Check that a square is placed into the correct position in a domain
TEST(particle_container, SquarePlacement) {
	GTEST_CXP_GCC auto ok = std::invoke([] {
		const vec domain{10, 10, 10};
		particle_container container(domain, 1);
		std::size_t seq_no = 0;
		container.add_cuboid<2>(
			cuboid_parameters<2>{.origin{5, 5}, .scale{2, 2}}.extend_to_3d(domain),
			{.meshwidth = 1, .brownian_mean = 5.0}, {}, {.mass = 100, .sigma = 1.0, .epsilon = 5.0},
			seq_no
		);
		const auto first = container.cell_containing({5, 5, 5}).at(0);
		const auto second = container.cell_containing({6, 6, 5}).at(0);
		const auto third = container.cell_containing({6, 5, 5}).at(0);
		return std::array{
			container.material_for_particle(first).mass == 100,
			container.material_for_particle(second).mass == 100,
			container.material_for_particle(third).mass == 100
		};
	});
	GCC_STATIC_EXPECT_ALL(ok);
}

// Check that a cuboid is placed into the correct position in a domain
TEST(particle_container, Cuboid3DPlacement) {
	GTEST_CXP_GCC auto ok = std::invoke([] {
		particle_container container(vec{10, 10, 10}, 1);
		std::size_t seq_no = 0;
		container.add_cuboid<3>(
			{.origin{5, 5, 5}, .scale{2, 2, 2}}, {.meshwidth = 1, .brownian_mean = 5.0}, {},
			{.mass = 100, .sigma = 1.0, .epsilon = 5.0}, seq_no
		);
		const auto first = container.cell_containing({5, 6, 5}).at(0);
		const auto second = container.cell_containing({6, 6, 6}).at(0);
		const auto third = container.cell_containing({6, 5, 6}).at(0);

		return std::array{
			container.material_for_particle(first).mass == 100,
			container.material_for_particle(second).mass == 100,
			container.material_for_particle(third).mass == 100
		};
	});
	STATIC_EXPECT_ALL(ok);
}

// NOLINTEND(*magic-numbers)
