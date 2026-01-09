#include "daw/json/daw_json_link.h"
#include "gtest_constexpr.h"
#include "ParticleContainer.h"
#include "gtest/gtest.h"
#include <utility>

#include "FileReader.h"

#include "SharedTools.h"

// NOLINTBEGIN(*magic-numbers*)

TEST(ConfigTest, Vec3DComponent) {
	GTEST_CXP auto doubles = daw::json::from_json<vec>("[0.1, 0.2, 0.3]");
	STATIC_EXPECT_EQ(doubles, vec(0.1, 0.2, 0.3));
	using int_vec = vec_3d<ParticleContainer::size_type>;
	GTEST_CXP auto integers = daw::json::from_json<int_vec>("[1, 2, 3]");
	STATIC_EXPECT_EQ(integers, int_vec(1, 2, 3));
}

// TODO(tuna): runtime test for catching mismatch of double <-> int
TEST(ConfigTest, Vec2DComponent) {
	GTEST_CXP auto doubles = daw::json::from_json<vec_2d<double>>("[0.1, 0.2]");
	STATIC_EXPECT_EQ(doubles, vec(0.1, 0.2, 0));
	using size_type = ParticleContainer::size_type;
	GTEST_CXP auto integers = daw::json::from_json<vec_2d<size_type>>("[1, 2]");
	STATIC_EXPECT_EQ(integers, vec_3d<size_type>(1, 2, 0));
}

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

// TODO(tuna): move config out of the ok-paradigm when we move to fixed-string for its base_name storage
TEST(ConfigTest, SimConfigurationComponent) {
	GTEST_CXP auto ok = std::invoke([&] {
		sim_configuration sim_config = daw::json::from_json<sim_configuration>(R"({
"delta_t": 0.1,
"cutoff_radius": 0.2,
"sigma": 0.7,
"epsilon": 5,
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
"domain": [0.4, 0.5, 0.6]
})");
		using enum boundary_condition;
		using enum boundary_type;
		return std::array{
			sim_config.delta_t == 0.1,
			sim_config.cutoff_radius == 0.2,
			sim_config.sigma == 0.7,
			sim_config.epsilon == 5,
			sim_config.boundary_behavior[x_min] == reflecting,
			sim_config.boundary_behavior.y_min == outflow,
			sim_config.boundary_behavior[z_min] == reflecting,
			sim_config.boundary_behavior.x_max == outflow,
			sim_config.boundary_behavior[y_max] == reflecting,
			sim_config.boundary_behavior.z_max == outflow,
			sim_config.end_time == 0.3,
			sim_config.write_frequency == 10,
			sim_config.base_name == "base_name",
			sim_config.domain == vec{0.4, 0.5, 0.6}
		};
	});
	STATIC_EXPECT_ALL(ok);
}

TEST(ConfigTest, Cuboid3DParametersComponent) {
	GTEST_CXP auto params = daw::json::from_json<cuboid_parameters<3>>(R"({
"origin": [0.1, 0.2, 0.3],
"scale": [4, 5, 6],
"velocity": [0.7, 0.8, 0.9],
"brownian_mean": 1.1,
"particle_mass": 1.2,
})");
	STATIC_EXPECT_EQ(params.origin, vec(0.1, 0.2, 0.3));
	STATIC_EXPECT_EQ(params.scale, ParticleContainer::index(4, 5, 6));
	STATIC_EXPECT_EQ(params.velocity, vec(0.7, 0.8, 0.9));
	STATIC_EXPECT_DOUBLE_EQ(params.brownian_mean, 1.1);
	STATIC_EXPECT_DOUBLE_EQ(params.particle_mass, 1.2);
}

TEST(ConfigTest, Cuboid2DParametersComponent) {
	GTEST_CXP auto params = daw::json::from_json<cuboid_parameters<2>>(R"({
"origin": [0.1, 0.2],
"scale": [4, 5],
"velocity": [0.7, 0.8],
"brownian_mean": 1.1,
"particle_mass": 1.2
})");
	STATIC_EXPECT_EQ(params.origin, vec(0.1, 0.2, 0));
	STATIC_EXPECT_EQ(params.scale, ParticleContainer::index(4, 5, 0));
	STATIC_EXPECT_EQ(params.velocity, vec(0.7, 0.8, 0));
	STATIC_EXPECT_DOUBLE_EQ(params.brownian_mean, 1.1);
	STATIC_EXPECT_DOUBLE_EQ(params.particle_mass, 1.2);
}

TEST(ConfigTest, ParticleComponent) {
	GTEST_CXP auto particle = daw::json::from_json<Particle>(R"({
"position": [0.1, 0.2, 0.3],
"velocity": [0.4, 0.5, 0.6],
"mass": 0.7,
"type": 8
})");
	STATIC_EXPECT_EQ(particle.x, vec(0.1, 0.2, 0.3));
	STATIC_EXPECT_EQ(particle.v, vec(0.4, 0.5, 0.6));
	STATIC_EXPECT_DOUBLE_EQ(particle.m, 0.7);
	STATIC_EXPECT_EQ(particle.type, 8);
}

TEST(ConfigTest, BasicConfig) {
	GTEST_CXP std::string_view json_data = R"({
  "configuration": {
    "delta_t": 0.1,
    "cutoff_radius": 1.0,
    "sigma": 3.4,
    "epsilon": 5,
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
    ]
  },
  "bodies": [
    {
      "type": "cuboid_3d",
      "parameters": {
        "origin": [0, 0, 0],
        "scale": [2, 2, 2],
        "velocity": [0.1, 0.1, 0.1],
        "brownian_mean": 0.1,
        "particle_mass": 10
      }
    },
    {
      "type": "particle",
      "parameters": {
        "position": [2.5, 2.5, 2.5],
        "velocity": [0.2, 0.3, 0.4],
        "mass": 0.1,
        "type": 1
      }
    }
  ]
})";
	GTEST_CXP_GCC auto ok = std::invoke([&] {
		using enum boundary_condition;
		auto [cfg, particles] = FileReader::parse(json_data);
		return std::array{
			cfg.delta_t == 0.1,
			cfg.cutoff_radius == 1.0,
			cfg.sigma == 3.4,
			cfg.epsilon == 5,
			cfg.boundary_behavior == boundary_conditions_descriptor::all(boundary_condition::reflecting),
			cfg.end_time == 0.3,
			cfg.base_name == "MD_vtk",
			cfg.write_frequency == 10,
			cfg.domain == vec{13.1, 13.2, 13.3},
			particles.cell_containing({2.5, 2.5, 2.5})[0].v == vec{0.2, 0.3, 0.4},
			particles.cell_containing({0.2, 0.2, 0.2})[0].m == 10
		};
	});
	GCC_STATIC_EXPECT_ALL(ok);
}

/*
TEST(ConfigTest, EmbedTest) {
    // NOLINTNEXTLINE(*avoid-c-arrays)
    static constexpr char data[] {
#embed "../../input/assignment3/task2/collision_of_two_bodies.json"
        , 0
    };
    constexpr auto json_data = std::string_view(data);
    GTEST_CXP_GCC auto ok = std::invoke([&] {
        auto [cfg, particles] = FileReader::parse(json_data);
    });
}
*/
// NOLINTEND(*magic-numbers*)
