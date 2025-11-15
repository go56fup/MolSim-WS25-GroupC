#include "daw/json/daw_json_link.h"
#include "gtest_constexpr.h"
#include "ParticleContainer.h"
#include "gtest/gtest.h"
#include <utility>

#include "FileReader.h"

template class flag_special_member_funcs<"Particle">;

#define INTERNAL_COMPONENT_CMP(i, comparison, a, b, ...) comparison((a).i, (b).i __VA_OPT__(, ) __VA_ARGS__)
#define INTERNAL_VEC_COMPONENTWISE_CMP(...)                                                                            \
	INTERNAL_COMPONENT_CMP(x, __VA_ARGS__);                                                                            \
	INTERNAL_COMPONENT_CMP(y, __VA_ARGS__);                                                                            \
	INTERNAL_COMPONENT_CMP(z, __VA_ARGS__)

#define GCC_STATIC_EXPECT_VEC_NEAR(a, b, abs_error)                                                                    \
	INTERNAL_VEC_COMPONENTWISE_CMP(GCC_STATIC_EXPECT_NEAR, a, b, abs_error)
#define GCC_STATIC_EXPECT_VEC_DOUBLE_EQ(a, b) INTERNAL_VEC_COMPONENTWISE_CMP(GCC_STATIC_EXPECT_DOUBLE_EQ, a, b)

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

// TODO(tuna): maybe not only bool but also implicitly convertible to bool, also maybe constrain the auto in
// Results
template <auto Results>
void expect_all() {
	[&]<std::size_t... I>(std::index_sequence<I...>) {
		(..., ([&] { STATIC_EXPECT_TRUE(Results[I]); }()));
	}(std::make_index_sequence<std::tuple_size_v<decltype(Results)>>{});
}

template <std::size_t N>
void expect_all(const std::array<bool, N>& results) {
	for (std::size_t i = 0; i < results.size(); ++i) {
		EXPECT_TRUE(results[i]) << "testing at index: " << i;
	}
}

#ifndef GTEST_CONFIG_RUNTIME_STATIC_EXPECT
#define STATIC_EXPECT_ALL(arr) expect_all<arr>()
#else
#define STATIC_EXPECT_ALL(arr) expect_all(arr)
#endif

TEST(ConfigTest, SimConfigurationComponent) {
	GTEST_CXP auto ok = std::invoke([&] {
		sim_configuration sim_config = daw::json::from_json<sim_configuration>(R"({
"delta_t": 0.1,
"cutoff_radius": 0.2,
"end_time": 0.3,
"write_frequency": 10,
"base_name": "base_name",
"domain": [0.4, 0.5, 0.6]
})");
		return std::array{sim_config.delta_t == 0.1,           sim_config.cutoff_radius == 0.2,
		                  sim_config.end_time == 0.3,          sim_config.write_frequency == 10,
		                  sim_config.base_name == "base_name", sim_config.domain == vec{0.4, 0.5, 0.6}};
	});
	STATIC_EXPECT_ALL(ok);
}

TEST(ConfigTest, CuboidParticleParametersComponent) {
	GTEST_CXP auto params = daw::json::from_json<cuboid_particle_parameters>(R"({
"meshwidth": 0.1,
"mass": 0.2
})");
	STATIC_EXPECT_DOUBLE_EQ(params.meshwidth, 0.1);
	STATIC_EXPECT_DOUBLE_EQ(params.mass, 0.2);
}

TEST(ConfigTest, Cuboid3DParametersComponent) {
	GTEST_CXP auto params = daw::json::from_json<cuboid_parameters<3>>(R"({
"origin": [0.1, 0.2, 0.3],
"scale": [4, 5, 6],
"velocity": [0.7, 0.8, 0.9],
"brownian_mean": 1.1,
"particles": {
"meshwidth": 1.2,
"mass": 1.3
}
})");
	STATIC_EXPECT_EQ(params.origin, vec(0.1, 0.2, 0.3));
	STATIC_EXPECT_EQ(params.scale, ParticleContainer::index(4, 5, 6));
	STATIC_EXPECT_EQ(params.velocity, vec(0.7, 0.8, 0.9));
	STATIC_EXPECT_DOUBLE_EQ(params.brownian_mean, 1.1);
	GTEST_CXP auto expected_params = cuboid_particle_parameters{.meshwidth = 1.2, .mass = 1.3};
	STATIC_EXPECT_EQ(params.particle_params, expected_params);
}

TEST(ConfigTest, Cuboid2DParametersComponent) {
	GTEST_CXP auto params = daw::json::from_json<cuboid_parameters<2>>(R"({
"origin": [0.1, 0.2],
"scale": [4, 5],
"velocity": [0.7, 0.8],
"brownian_mean": 1.1,
"particles": {
"meshwidth": 1.2,
"mass": 1.3
}
})");
	STATIC_EXPECT_EQ(params.origin, vec(0.1, 0.2, 0));
	STATIC_EXPECT_EQ(params.scale, ParticleContainer::index(4, 5, 0));
	STATIC_EXPECT_EQ(params.velocity, vec(0.7, 0.8, 0));
	STATIC_EXPECT_DOUBLE_EQ(params.brownian_mean, 1.1);
	GTEST_CXP auto expected_params = cuboid_particle_parameters{.meshwidth = 1.2, .mass = 1.3};
	STATIC_EXPECT_EQ(params.particle_params, expected_params);
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
    "end_time": 0.3,
    "write_frequency": 10,
    "base_name": "MD_vtk",
    "domain": [
      3,
      3,
      3
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
        "particles": {
          "meshwidth": 0.2,
          "mass": 10
        }
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
	GTEST_CXP auto ok = std::invoke([&] {
		auto [cfg, particles] = FileReader::parse(json_data);
		return std::array{
			cfg.delta_t == 0.1,
			cfg.cutoff_radius == 1.0,
			cfg.end_time == 0.3,
			cfg.base_name == "MD_vtk",
			cfg.write_frequency == 10,
			cfg.domain == vec{3, 3, 3},
			particles.cell_containing({2.5, 2.5, 2.5})[0].v == vec{0.2, 0.3, 0.4},
			particles.cell_containing({0.2, 0.2, 0.2})[0].m == 10
		};
	});
	STATIC_EXPECT_ALL(ok);
}
