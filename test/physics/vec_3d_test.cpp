#include <functional>

#include <gtest/gtest.h>
#include "gtest_constexpr/macros.hpp"
#include "testing_utilities.hpp"

#include "grid/particle_container/fwd.hpp"
#include "physics/thermostat.hpp"
#include "physics/vec_3d.hpp"
#include "simulation/molsim.hpp"

// NOLINTBEGIN(*magic-numbers)

// Check that vector addition works
TEST(VectorTests, Addition) {
	GTEST_CXP vec a{1, 2, 3};
	GTEST_CXP vec b{4, 5, 6};
	GTEST_CXP vec result{5, 7, 9};
	STATIC_EXPECT_EQ(a + b, result);
}

// Check that vector subtraction works
TEST(VectorTests, Subtraction) {
	GTEST_CXP vec a{5, 7, 9};
	GTEST_CXP vec b{4, 5, 6};
	GTEST_CXP vec result{1, 2, 3};
	STATIC_EXPECT_EQ(a - b, result);
}

// Check that vector plus-equals works
TEST(VectorTests, AdditionAssignment) {
	GTEST_CXP vec result = std::invoke([] {
		vec a{1, 2, 3};
		vec b{4, 5, 6};
		a += b;
		return a;
	});
	GTEST_CXP vec expected{5, 7, 9};
	STATIC_EXPECT_EQ(result, expected);
}

// Check that vector minus-equals works
TEST(VectorTests, SubtractionAssignment) {
	GTEST_CXP vec result = std::invoke([] {
		vec a{4, 5, 6};
		vec b{1, 2, 3};
		a -= b;
		return a;
	});
	GTEST_CXP vec expected{3, 3, 3};
	STATIC_EXPECT_EQ(result, expected);
}

// Check that vector scaling works
TEST(VectorTests, Scalar) {
	GTEST_CXP vec a{5, 7, 9};
	GTEST_CXP vec result{10, 14, 18};
	STATIC_EXPECT_EQ(a * 2, result);
	STATIC_EXPECT_EQ(2 * a, result);
}

// Check that vector times-equals works
TEST(VectorTests, ScalarAssignment) {
	GTEST_CXP vec result = std::invoke([] {
		vec a{5, 7, 9};
		a *= 2;
		return a;
	});

	GTEST_CXP vec expected{10, 14, 18};
	STATIC_EXPECT_EQ(result, expected);
}

// Check that vector normalization works
TEST(VectorTests, Norm) {
	GTEST_CXP vec a{2, 3, 6};
	GTEST_CXP double result = 7;
	GCC_STATIC_EXPECT_EQ(a.euclidian_norm(), result);
}

// Check that vector equality works
TEST(VectorTests, Equal) {
	GTEST_CXP vec a{1, 2, 3};
	GTEST_CXP vec b{1, 2, 3};
	STATIC_EXPECT_EQ(a, b);
}

TEST(ThermostatTests, Cooling) {
	GTEST_CXP double init_temp = 40;
	GTEST_CXP double target_temp = 30;
	GTEST_CXP double mass = 2.0;
	GTEST_CXP double brownian_mean = std::sqrt(init_temp/ mass);
	GTEST_CXP particle_container container(vec{2, 2, 2}, 1);
	GTEST_CXP thermostat_parameters therm = {init_temp, 2, target_temp, std::numeric_limits<double>::infinity()};
	GTEST_CXP size_t seq = 2;

	container.add_cuboid<2>({0,0,0}, {1,1,0}, 1.1225, {0,0,0}, 2.0, 1.0, 5.0, brownian_mean, seq);
	run_thermostat(container, therm, 2);

	double resulting_temp = 0;
	for (const particle& p : container.particles()) {
		const double squared_norm = (p.v.x * p.v.x) + (p.v.y * p.v.y) + (p.v.z * p.v.z);
		resulting_temp += p.m * squared_norm;
	}
	resulting_temp /= 2 * static_cast<double>(container.size());
	STATIC_EXPECT_EQ(resulting_temp, target_temp);
}

TEST(ThermostatTests, Heating) {
	GTEST_CXP double init_temp = 30;
	GTEST_CXP double target_temp = 40;
	GTEST_CXP double mass = 2.0;
	GTEST_CXP double brownian_mean = std::sqrt(init_temp/ mass);
	GTEST_CXP particle_container container(vec{2, 2, 2}, 1);
	GTEST_CXP thermostat_parameters therm = {init_temp, 2, target_temp, std::numeric_limits<double>::infinity()};
	GTEST_CXP size_t seq = 2;

	container.add_cuboid<2>({0,0,0}, {1,1,0}, 1.1225, {0,0,0}, 2.0, 1.0, 5.0, brownian_mean, seq);
	run_thermostat(container, therm, 2);

	double resulting_temp = 0;
	for (const particle& p : container.particles()) {
		const double squared_norm = (p.v.x * p.v.x) + (p.v.y * p.v.y) + (p.v.z * p.v.z);
		resulting_temp += p.m * squared_norm;
	}
	resulting_temp /= 2 * static_cast<double>(container.size());
	STATIC_EXPECT_EQ(resulting_temp, target_temp);
}

TEST(ThermostatTests, Holding) {
	GTEST_CXP double init_temp = 40;
	GTEST_CXP double mass = 2.0;
	GTEST_CXP double brownian_mean = std::sqrt(init_temp/ mass);
	GTEST_CXP particle_container container(vec{2, 2, 2}, 1);
	GTEST_CXP thermostat_parameters therm = {init_temp, 2, init_temp, std::numeric_limits<double>::infinity()};
	GTEST_CXP size_t seq = 2;
	GTEST_CXP auto config = sim_configuration{
		.delta_t = 0.1,
		.cutoff_radius = 1,
		.boundary_behavior = boundary_conditions_descriptor::all(boundary_condition::outflow),
		.end_time = 10,
		.write_frequency = 1,
		.base_name{std::from_range, "unused"},
		.domain{2, 2, 2},
		.create_checkpoint = false,
		.dimensions = 2
	};

	container.add_cuboid<2>({0,0,0}, {1,1,0}, 1.1225, {0,0,0}, 2.0, 1.0, 5.0, brownian_mean, seq);
	for (int i = 0; i< 10; i++) {
		run_thermostat(container, therm, 2);
		run_sim_iteration(lennard_jones_force, container, config);
	}

	double resulting_temp = 0;
	for (const particle& p : container.particles()) {
		const double squared_norm = (p.v.x * p.v.x) + (p.v.y * p.v.y) + (p.v.z * p.v.z);
		resulting_temp += p.m * squared_norm;
	}
	resulting_temp /= 2 * static_cast<double>(container.size());
	STATIC_EXPECT_EQ(resulting_temp, init_temp);
}


// NOLINTEND(*magic-numbers)
