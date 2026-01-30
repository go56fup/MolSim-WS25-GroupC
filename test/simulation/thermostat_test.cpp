#include "gtest_constexpr/macros.hpp"

#include "testing_utilities.hpp"
#include <gtest/gtest.h>

#include "grid/particle_container/fwd.hpp"
#include "simulation/config/entities.hpp"
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
			{0, 0, 0}, {1, 1, 1}, 1.1225, {0, 0, 0}, mass, 1.0, 5.0, brownian_mean, seq
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
		.max_temperature_difference = std::numeric_limits<double>::infinity(),
		.enforce_initial_temperature = true
	};
	GTEST_CXP_GCC double resulting_temp = std::invoke([&] {
		particle_container container(domain, cutoff_radius);
		std::size_t seq = 0;

		container.add_cuboid<2>(
			{0, 0, 0}, {1, 1, 1}, 1.1225, {0, 0, 0}, mass, 1.0, 5.0, brownian_mean, seq
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
		.max_temperature_difference = std::numeric_limits<double>::infinity(),
		.enforce_initial_temperature = true
	};

	GTEST_CXP auto config = sim_configuration{
		.delta_t = 0.1,
		.cutoff_radius = 1,
		.boundary_behavior = boundary_conditions_descriptor::all(boundary_condition::reflecting),
		.thermostat = therm,
		.end_time = 10,
		.write_frequency = 1,
		.base_name{std::from_range, "unused"},
		.domain{5, 5, 5},
		.create_checkpoint = false,
		.dimensions = 3,
		.gravitational_constant = 0
	};

	GTEST_CXP_GCC double resulting_temp = std::invoke([&] {
		particle_container container(config.domain, config.cutoff_radius);
		std::size_t seq = 0;
		container.add_cuboid<3>(
			{3, 3, 3}, {1, 1, 1}, 1.1225, {0, 0, 0}, 2.0, 1.0, 5.0, brownian_mean, seq
		);
		for (sim_iteration_t i = 0; i < 10; i++) {
			run_sim_iteration(lennard_jones_force, container, config, i);
		}
		return get_temperature(container, config.dimensions);
	});

	STATIC_EXPECT_DOUBLE_EQ(resulting_temp, init_temp);
}
