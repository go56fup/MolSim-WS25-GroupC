#include "gtest_constexpr/macros.hpp"

#include "testing_utilities.hpp"
#include <gtest/gtest.h>

#include "grid/particle_container/fwd.hpp"
#include "simulation/entities.hpp"
#include "simulation/molsim.hpp"
#include "simulation/thermostat.hpp"

// Check that the thermostat can cool the simulation down
TEST(ThermostatTests, Cooling) {
	static constexpr double init_temp = 40;
	static constexpr double mass = 2.0;
	GTEST_CXP_GCC double brownian_mean = std::sqrt(init_temp / mass);
	static constexpr double target_temp = 30;
	static constexpr auto dims = 2;
	static constexpr vec domain{2, 2, 2};
	static constexpr double cutoff_radius = 1;

	static constexpr thermostat_parameters therm{
		.initial_temperature = init_temp,
		.application_frequency = 1,
		.target_temperature = target_temp,
		.max_temperature_difference = std::numeric_limits<double>::infinity(),
		.enforce_initial_temperature = true
	};
	GTEST_CXP_GCC double resulting_temp = std::invoke([&] {
		particle_container container(domain, cutoff_radius);
		std::size_t seq = 0;

		container.add_cuboid<2>(
			cuboid_parameters<2>{.origin{0, 0}, .scale{1, 1}}.extend_to_3d(domain),
			{.meshwidth = 1.1225, .brownian_mean = brownian_mean}, {0, 0, 0},
			{.mass = mass, .sigma = 1, .epsilon = 5}, seq
		);

		run_thermostat(container, therm, dims);
		return get_temperature(container, dims);
	});
	STATIC_EXPECT_DOUBLE_EQ(resulting_temp, target_temp);
}

// Check that the thermostat can heat the simulation up
TEST(ThermostatTests, Heating) {
	static constexpr double init_temp = 30;
	static constexpr double target_temp = 40;
	static constexpr double mass = 2.0;
	GTEST_CXP_GCC double brownian_mean = std::sqrt(init_temp / mass);
	static constexpr auto dims = 2;
	static constexpr vec domain{2, 2, 2};
	static constexpr double cutoff_radius = 1;

	static constexpr thermostat_parameters therm{
		.initial_temperature = init_temp,
		.application_frequency = 1,
		.target_temperature = target_temp,
		.enforce_initial_temperature = true
	};
	GTEST_CXP_GCC double resulting_temp = std::invoke([&] {
		particle_container container(domain, cutoff_radius);
		std::size_t seq = 0;

		container.add_cuboid<2>(
			cuboid_parameters<2>{.origin{0, 0}, .scale{1, 1}}.extend_to_3d(domain),
			{.meshwidth = 1.1225, .brownian_mean = brownian_mean}, {},
			{.mass = mass, .sigma = 1.0, .epsilon = 5.0}, seq
		);

		run_thermostat(container, therm, dims);
		return get_temperature(container, dims);
	});

	STATIC_EXPECT_DOUBLE_EQ(resulting_temp, target_temp);
}

// Check that the thermostat is able to hold the temperature of the simulation constant over
// multiple iterations
TEST(ThermostatTests, Holding) {
	static constexpr double init_temp = 2;
	static constexpr double mass = 2.0;
	GTEST_CXP_GCC double brownian_mean = std::sqrt(init_temp / mass);

	GTEST_CXP thermostat_parameters therm{
		.initial_temperature = init_temp,
		.application_frequency = 1,
		.target_temperature = init_temp,
		.enforce_initial_temperature = true
	};

	GTEST_CXP auto config = sim_configuration{
		.delta_t = 0.1,
		.cutoff_radius = 1,
		.boundary_behavior = boundary_conditions_descriptor::all(boundary_condition::reflecting),
		.thermostat = therm,
		.end_time = 10,
		.domain{5, 5, 5},
		.dimensions = 3,
	};

	double resulting_temp = std::invoke([&] {
		particle_container container(config.domain, config.cutoff_radius);
		std::size_t seq = 0;
		container.add_cuboid<3>(
			{.origin{3, 3, 3}, .scale{1, 1, 1}},
			{.meshwidth = 1.1225, .brownian_mean = brownian_mean}, {},
			{.mass = 2.0, .sigma = 1.0, .epsilon = 5.0}, seq
		);

		static constexpr sim_iteration_t hold_temperature_over_n_iterations = 10;
		for (sim_iteration_t i = 0; i < hold_temperature_over_n_iterations; ++i) {
			run_sim_iteration(container, config, i);
		}
		return get_temperature(container, config.dimensions);
	});

	// TODO(tuna): this test has interesting behavior. it returns false at compile time,
	// and fails when run with other tests; but passes on its own.
	EXPECT_DOUBLE_EQ(resulting_temp, init_temp);
}
