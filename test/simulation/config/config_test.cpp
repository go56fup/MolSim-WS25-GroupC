#include <limits>
#include <utility>

#include "daw/json/daw_json_link.h"
#include "grid/particle_container/system.hpp"
#include "gtest_constexpr/macros.hpp"
#include <gtest/gtest.h>

#include "grid/particle_container/particle_container.hpp"
#include "physics/particle.hpp"
#include "simulation/config/entities.hpp"
#include "simulation/config/json_schema.hpp"
#include "simulation/config/parse.hpp"

#include "testing_utilities.hpp"

// NOLINTBEGIN(*magic-numbers*)

// Check that a 3D vector description parses to the correct vector
TEST(ConfigTest, Vec3DComponent) {
	GTEST_CXP auto doubles = daw::json::from_json<vec>("[0.1, 0.2, 0.3]");
	STATIC_EXPECT_EQ(doubles, vec(0.1, 0.2, 0.3));
	using int_vec = vec_3d<particle_container::size_type>;
	GTEST_CXP auto integers = daw::json::from_json<int_vec>("[1, 2, 3]");
	STATIC_EXPECT_EQ(integers, int_vec(1, 2, 3));
}

// Check that a vector's components can be specified as ints but result in doubles
TEST(ConfigTest, VecIntSpecifiable) {
	GTEST_CXP auto doubles = daw::json::from_json<vec>("[1, 2, 3]");
	STATIC_EXPECT_EQ(doubles, vec(1.0, 2.0, 3.0));
}

// Check that a 2D vector description parses to the correct vector
TEST(ConfigTest, Vec2DComponent) {
	GTEST_CXP auto doubles = daw::json::from_json<vec_2d<double>>("[0.1, 0.2]");
	STATIC_EXPECT_EQ(doubles, vec_2d(0.1, 0.2));
	using size_type = particle_container::size_type;
	GTEST_CXP auto integers = daw::json::from_json<vec_2d<size_type>>("[1, 2]");
	STATIC_EXPECT_EQ(integers, vec_2d<size_type>(1, 2));
}

// Check that a boundary condition description parses to the correct boundaries
TEST(ConfigTest, BoundaryConditionsDescriptorComponent) {
	GTEST_CXP auto params = daw::json::from_json<boundary_conditions_descriptor>(R"({
"x_min": "reflecting",
"y_min": "outflow",
"z_min": "reflecting",
"x_max": "outflow",
"y_max": "reflecting",
"z_max": "outflow"
})");
	using enum boundary_condition;
	STATIC_EXPECT_EQ(params.x_min, reflecting);
	STATIC_EXPECT_EQ(params.y_min, outflow);
	STATIC_EXPECT_EQ(params.z_min, reflecting);
	STATIC_EXPECT_EQ(params.x_max, outflow);
	STATIC_EXPECT_EQ(params.y_max, reflecting);
	STATIC_EXPECT_EQ(params.z_max, outflow);
}

// Check that a simulation description parses to the correct simulation parameters
TEST(ConfigTest, SimConfigurationComponent) {
	GTEST_CXP auto sim_config = daw::json::from_json<sim_configuration>(R"({
"delta_t": 0.1,
"cutoff_radius": 0.2,
"boundary_conditions": {
"x_min": "reflecting",
"y_min": "outflow",
"z_min": "reflecting",
"x_max": "outflow",
"y_max": "reflecting",
"z_max": "outflow"
},
"end_time": 0.3,
"write_frequency": 10,
"base_name": "base_name",
"domain": [0.4, 0.5, 0.6],
"create_checkpoint": true
})");
	using enum boundary_condition;
	using enum boundary_type;
	STATIC_EXPECT_EQ(sim_config.delta_t, 0.1);
	STATIC_EXPECT_EQ(sim_config.cutoff_radius, 0.2);
	STATIC_EXPECT_EQ(sim_config.boundary_behavior[x_min], reflecting);
	STATIC_EXPECT_EQ(sim_config.boundary_behavior.y_min, outflow);
	STATIC_EXPECT_EQ(sim_config.boundary_behavior[z_min], reflecting);
	STATIC_EXPECT_EQ(sim_config.boundary_behavior.x_max, outflow);
	STATIC_EXPECT_EQ(sim_config.boundary_behavior[y_max], reflecting);
	STATIC_EXPECT_EQ(sim_config.boundary_behavior.z_max, outflow);
	STATIC_EXPECT_EQ(sim_config.end_time, 0.3);
	STATIC_EXPECT_EQ(sim_config.write_frequency, 10);
	STATIC_EXPECT_STREQ(sim_config.base_name.c_str(), "base_name");
	STATIC_EXPECT_EQ(sim_config.domain, vec(0.4, 0.5, 0.6));
	STATIC_EXPECT_TRUE(sim_config.create_checkpoint);
	STATIC_EXPECT_DOUBLE_EQ(sim_config.gravitational_constant, 0);
}

// Check that a cuboid description parses to the correct cuboid
TEST(ConfigTest, Cuboid3DParametersComponent) {
	GTEST_CXP auto params = daw::json::from_json<cuboid_parameters<3>>(R"({
"origin": [0.1, 0.2, 0.3],
"scale": [4, 5, 6],
})");
	STATIC_EXPECT_EQ(params.origin, vec(0.1, 0.2, 0.3));
	STATIC_EXPECT_EQ(params.scale, particle_container::index(4, 5, 6));
}

// Check that a square description parses to the correct rectangle
TEST(ConfigTest, RectangleParametersComponent) {
	GTEST_CXP auto params = daw::json::from_json<cuboid_parameters<2>>(R"({
"origin": [0.1, 0.2],
"scale": [4, 5]
})");
	STATIC_EXPECT_EQ(params.origin, vec_2d<double>(0.1, 0.2));
	STATIC_EXPECT_EQ(params.scale, vec_2d<particle_container::size_type>(4, 5));
}

// Check that a particle description parses to the correct particle
TEST(ConfigTest, ParticleComponent) {
	GTEST_CXP auto p = daw::json::from_json<particle_parameters>(R"({
"position": [0.1, 0.2, 0.3]
})");
	STATIC_EXPECT_EQ(p.position, vec(0.1, 0.2, 0.3));
}

// Check that a thermostat description parses to the correct parameters to the simulation
TEST(ConfigTest, ThermostatComponent) {
	GTEST_CXP auto p = daw::json::from_json<thermostat_parameters>(R"({
"initial_temperature": 0.1,
"application_frequency": 2,
"target_temperature": 0.2,
"max_temperature_difference": 0.3,
"enforce_initial_temperature": false
})");
	STATIC_EXPECT_DOUBLE_EQ(p.initial_temperature, 0.1);
	STATIC_EXPECT_EQ(p.application_frequency, 2);
	STATIC_EXPECT_DOUBLE_EQ(p.target_temperature, 0.2);
	STATIC_EXPECT_DOUBLE_EQ(p.max_temperature_difference, 0.3);
	STATIC_EXPECT_FALSE(p.enforce_initial_temperature);
}

// Check that user not specifying the temperature difference defaults to +inf
TEST(ConfigTest, OptionalTemperatureDifference) {
	GTEST_CXP auto p = daw::json::from_json<thermostat_parameters>(R"({
"initial_temperature": 0.1,
"application_frequency": 2,
"target_temperature": 0.2,
"enforce_initial_temperature": false
})");
	STATIC_EXPECT_DOUBLE_EQ(p.initial_temperature, 0.1);
	STATIC_EXPECT_EQ(p.application_frequency, 2);
	STATIC_EXPECT_DOUBLE_EQ(p.target_temperature, 0.2);
	STATIC_EXPECT_DOUBLE_EQ(p.max_temperature_difference, std::numeric_limits<double>::infinity());
	STATIC_EXPECT_FALSE(p.enforce_initial_temperature);
}

// Check that user not specifying the target temperature defaults to initial temperature
TEST(ConfigTest, OptionalTargetTemperature) {
	GTEST_CXP auto p = daw::json::from_json<thermostat_parameters>(R"({
"initial_temperature": 0.1,
"application_frequency": 2,
"enforce_initial_temperature": true
})");
	STATIC_EXPECT_DOUBLE_EQ(p.initial_temperature, 0.1);
	STATIC_EXPECT_EQ(p.application_frequency, 2);
	STATIC_EXPECT_DOUBLE_EQ(p.target_temperature, 0.1);
	STATIC_EXPECT_TRUE(p.enforce_initial_temperature);
}

// Check that user can specify null explicitly to get the default temperature difference
TEST(ConfigTest, TemperatureDifferenceNullDefault) {
	GTEST_CXP auto p = daw::json::from_json<thermostat_parameters>(R"({
"initial_temperature": 0.1,
"application_frequency": 2,
"target_temperature": 0.2,
"max_temperature_difference": null,
"enforce_initial_temperature": true
})");
	STATIC_EXPECT_DOUBLE_EQ(p.initial_temperature, 0.1);
	STATIC_EXPECT_EQ(p.application_frequency, 2);
	STATIC_EXPECT_DOUBLE_EQ(p.target_temperature, 0.2);
	STATIC_EXPECT_DOUBLE_EQ(p.max_temperature_difference, std::numeric_limits<double>::infinity());
	STATIC_EXPECT_TRUE(p.enforce_initial_temperature);
}

// Check that parsing a configuration file results in the correct parameters to the simulation
TEST(ConfigTest, BasicConfig) {
	GTEST_CXP std::string_view config_json_data = R"({
  "delta_t": 0.1,
  "cutoff_radius": 1.0,
  "boundary_conditions": {
    "x_min": "reflecting",
    "y_min": "reflecting",
    "z_min": "reflecting",
    "x_max": "reflecting",
    "y_max": "reflecting",
    "z_max": "reflecting"
  },
  "end_time": 0.3,
  "write_frequency": 10,
  "base_name": "MD_vtk",
  "domain": [
    13.1,
    13.2,
    13.3
  ],
  "create_checkpoint": true
})";
	using enum boundary_condition;

	GTEST_CXP auto cfg = config::parse_config(config_json_data);
	STATIC_EXPECT_EQ(cfg.delta_t, 0.1);
	STATIC_EXPECT_EQ(cfg.cutoff_radius, 1.0);
	STATIC_EXPECT_EQ(
		cfg.boundary_behavior, boundary_conditions_descriptor::all(boundary_condition::reflecting)
	);
	STATIC_EXPECT_EQ(cfg.end_time, 0.3);
	STATIC_EXPECT_STREQ(cfg.base_name.c_str(), "MD_vtk");
	STATIC_EXPECT_EQ(cfg.write_frequency, 10);
	STATIC_EXPECT_VEC_DOUBLE_EQ(cfg.domain, vec(13.1, 13.2, 13.3));
	STATIC_EXPECT_FALSE(cfg.thermostat.has_value());
	STATIC_EXPECT_TRUE(cfg.create_checkpoint);

	GTEST_CXP std::string_view body_json_data = R"([
{
      "type": "cuboid",
      "geometry": {
        "origin": [0, 0, 0],
        "scale": [2, 2, 2]
      },
      "velocity": [0.1, 0.1, 0.1],
      "parameters": {
        "meshwidth": 1.12,
        "brownian_mean": 0.1
       },
       "material": {
        "mass": 10,
        "sigma": 1,
        "epsilon": 5
       }
    },
    {
      "type": "particle",
      "geometry": {
        "position": [2.5, 2.5, 2.5]
      },
      "velocity": [0.2, 0.3, 0.4],
      "material": {
        "mass": 0.1,
        "sigma": 1,
        "epsilon": 5
      }
    }
])";

	GTEST_CXP_GCC auto ok = std::invoke([&] {
		particle_container container(cfg.domain, cfg.cutoff_radius);
		auto bodies = config::parse_bodies(body_json_data);
		config::populate_simulation(container, cfg, bodies);
		const particle_system::particle_id rectangle =
			container.cell_containing({2.5, 2.5, 2.5}).at(0);
		const particle_system::particle_id particle =
			container.cell_containing({0.2, 0.2, 0.2}).at(0);

		return std::array{
			container.system().serialize_velocity(particle) == vec{0.2, 0.3, 0.4},
			container.material_for_particle(particle).mass == 0.1,
			container.system().serialize_velocity(rectangle) == vec{0.1, 0.1, 0.1},
			container.material_for_particle(rectangle).mass == 10
		};
	});
	GCC_STATIC_EXPECT_ALL(ok);
}

// Check that supplying the optional thermostat key correctly parses relevant parameters
TEST(ConfigTest, ThermostatOptional) {
	GTEST_CXP std::string_view json_data = R"({
    "delta_t": 0.1,
    "cutoff_radius": 1.0,
    "boundary_conditions": {
      "x_min": "reflecting",
      "y_min": "reflecting",
      "z_min": "reflecting",
      "x_max": "reflecting",
      "y_max": "reflecting",
      "z_max": "reflecting"
    },
    "thermostat": {
      "initial_temperature": 1.1,
      "application_frequency": 2,
      "target_temperature": 1.3,
      "enforce_initial_temperature": true
    },

    "end_time": 0.3,
    "write_frequency": 10,
    "base_name": "MD_vtk",
    "domain": [
      13.1,
      13.2,
      13.3
    ],
    "create_checkpoint": true
})";
	using enum boundary_condition;

	GTEST_CXP auto cfg = config::parse_config(json_data);
	STATIC_EXPECT_TRUE(cfg.thermostat.has_value());
	STATIC_EXPECT_DOUBLE_EQ(cfg.thermostat->initial_temperature, 1.1);
	STATIC_EXPECT_EQ(cfg.thermostat->application_frequency, 2);
	STATIC_EXPECT_DOUBLE_EQ(cfg.thermostat->target_temperature, 1.3);
	STATIC_EXPECT_DOUBLE_EQ(
		cfg.thermostat->max_temperature_difference, std::numeric_limits<double>::infinity()
	);
	STATIC_EXPECT_TRUE(cfg.thermostat->enforce_initial_temperature);
}

#if HAS_EMBED
// Check that configuration parsing via #embed works
TEST(ConfigTest, EmbedTest) {
	// NOLINTNEXTLINE(*avoid-c-arrays)
	static constexpr char data[]{
#embed "assignment5/debug/simple/3_collision_of_two_bodies/config.json"
		, 0
	};
	constexpr auto json_data = std::string_view(data);
	GTEST_CXP auto parse_result = config::parse(json_data);
	GTEST_CXP const auto& cfg = parse_result.config;
	STATIC_EXPECT_DOUBLE_EQ(cfg.delta_t, 0.0005);
	STATIC_EXPECT_DOUBLE_EQ(cfg.cutoff_radius, 3);
	STATIC_EXPECT_EQ(cfg.domain, vec(180, 90, 3));
}
#endif

// NOLINTEND(*magic-numbers*)
