#include <array>

#include "grid/particle_container/fwd.hpp"
#include "gtest_constexpr/macros.hpp"
#include "testing_utilities.hpp"

#include "grid/particle_container/particle_container.hpp"
#include "physics/forces.hpp"
#include "physics/particle.hpp"
#include "simulation/molsim.hpp"
#include "utility/concepts.hpp"

// NOLINTBEGIN(*magic-numbers)
// TODO(tuna): uncomment after figuring out what to do with gravitation + boundaries

#if 0
TEST(ForceTests, BasicGravitation) {
	GTEST_CXP auto ok = std::invoke([] {
		auto config = sim_configuration{
			.delta_t = 0.0014,
			.cutoff_radius = 1.5,
			.sigma = 1,
			.end_time = 0.0014,
			.write_frequency = 1000,
			.base_name = "unused",
			.domain{2, 2, 2},
			.meshwidth = 1.225,
		};
		// TODO(tuna): make the constructor one arged, and the tests always gotta use this config first idiom
		ParticleContainer particles(config.domain.x, config.domain.y, config.domain.z, config.cutoff_radius);
		particles.emplace(vec{}, vec{}, 1, 10);
		particles.emplace(vec{1, 0, 0}, vec{}, 1, 15);
		run_sim_iteration(gravitational_force, particles, config);
		auto moving_par = particles.view() | std::views::filter([](const Particle& p) { return p.type == 15; });
		return moving_par.begin()->f == vec{1, 0, 0};
	});
	GCC_STATIC_EXPECT_TRUE(ok);
}
#endif

template <fwd_reference_to<particle_container> Container>
constexpr auto& particle_from_type(Container&& container, decltype(particle::type) type) {
	auto flat = std::forward<Container>(container).particles();
	return *std::ranges::find_if(flat, [type](const particle& p) { return p.type == type; });
}

// Check that the Lennard-Jones force is calculated correctly
TEST(ForceTests, LennardJones) {
	GTEST_CXP particle first(vec{4, 4, 4}, vec{}, 1, 10);
	GTEST_CXP particle second(vec{5, 6, 4}, vec{}, 1, 15);
	GTEST_CXP auto config = sim_configuration{
		.delta_t = 0.0014,
		.cutoff_radius = 5,
		.sigma = 1,
		.epsilon = 5,
		.boundary_behavior = boundary_conditions_descriptor::all(boundary_condition::reflecting),
		.end_time = 0.0014,
		.write_frequency = 1000,
		.base_name{std::from_range, "unused"},
		.domain{10, 10, 10},
		.meshwidth = 1.125
	};
	GTEST_CXP auto calc = [sigma = config.sigma,
	                       eps = config.epsilon](const particle& p1, const particle& p2) noexcept {
		return lennard_jones_force(p1, p2, sigma, eps);
	};
	GTEST_CXP auto first_force = std::invoke([&] -> vec {
		particle_container particles(config.domain, config.cutoff_radius);
		particles.place(first);
		particles.place(second);
		run_simulation(particles, config, calc, "unused");
		return particle_from_type(particles, 10).f;
	});

	GTEST_CXP auto direct_calc = lennard_jones_force(first, second, config.sigma, config.epsilon);
	STATIC_EXPECT_VEC_DOUBLE_EQ(first_force, direct_calc);
	static constexpr auto hand_calculated_result = vec(2952, 5904, 0);
	GTEST_CXP auto scaled_force = 15625 * first_force;
	STATIC_EXPECT_VEC_DOUBLE_EQ(hand_calculated_result, scaled_force);
}

#if HAS_EMBED
// Check that two particles colliding works
TEST(ConfigTest, GridTest) {
	// NOLINTNEXTLINE(*avoid-c-arrays)
	static constexpr char data[]{
#embed "assignment3/debug/two_particles_moving_towards_each_other.json"
		, 0
	};
	constexpr auto json_data = std::string_view(data);
	GTEST_CXP auto parse_result = config::parse(json_data);
	GTEST_CXP const auto& config = parse_result.config;
	GTEST_CXP_SIM auto particle_result = std::invoke([&] -> std::pair<particle, particle> {
		particle_container particles(config.domain, config.cutoff_radius);
		config::populate_simulation(container, parse_result);
		run_simulation(particles, config, calc, "unused");
		return {particle_from_type(particles, 10), particle_from_type(particles, 15)};
	});
	// TODO(tuna): test that the two particles actually collide at iteration ~400
}
#endif

/*
// Check that the Halley's comet simulation from worksheet 1 results in the correct simulation
// (specifically, that the comet returns to its original position modulo some deviation, and
// the comet passes by Earth)
TEST(SimTests, HalleysComet) {
    GTEST_CXP Particle sun{vec{}, vec{}, 1.0};
    GTEST_CXP Particle earth{vec{0., 1., 0.}, vec{-1., 0., 0.}, 3.0e-6};
    GTEST_CXP Particle jupiter{vec{0., 5.36, 0.}, vec{-0.45, 0., 0.}, 9.55e-4};
    GTEST_CXP Particle comet{vec{34.75, 0., 0.}, vec{0., 0.0296, 0.}, 1e-14};

    enum planet_ids : std::uint8_t { sun_i = 0, earth_i = 1, jupiter_i = 2, comet_i = 3 };

    GTEST_CXP_GCC auto run_sim = [&](std::size_t iteration) {
        std::array<Particle, 4> planets{sun, earth, jupiter, comet};
        GTEST_CXP double delta_t = 0.014;

        for (std::size_t i = 0; i < iteration; ++i) {
            run_sim_iteration(gravitational_force, planets, delta_t);
        }

        return planets;
    };
    GTEST_CXP auto until_earth_sighting = 16600;
    GTEST_CXP_GCC auto earth_sighting = run_sim(until_earth_sighting);
    GTEST_CXP_GCC const vec& earth_coords = earth_sighting[earth_i].x;
    GTEST_CXP_GCC const vec& comet_coords = earth_sighting[comet_i].x;

    // "During its 1986 appearance, Halley's nearest approach to Earth occurred (...) at a distance
of
    // 0.42 AU (...)"
    // source: https://science.nasa.gov/solar-system/comets/1p-halley/
    GCC_STATIC_EXPECT_VEC_NEAR(earth_coords, comet_coords, 0.16);

    GTEST_CXP auto until_return_to_initial_comet_pos = 33690;
    GTEST_CXP const vec& comet_initial = comet.x;
    GTEST_CXP_GCC vec comet_after_one_revolution =
run_sim(until_return_to_initial_comet_pos)[comet_i].x;

    GTEST_CXP double deviation_x = 1.4064506667406391;
    GTEST_CXP double deviation_y = -0.004078428002306117;
    GCC_STATIC_EXPECT_NEAR(comet_initial.x, comet_after_one_revolution.x - deviation_x, 0.001);
    GCC_STATIC_EXPECT_NEAR(comet_initial.y, comet_after_one_revolution.y - deviation_y, 0.001);
}

constexpr void nothing_io(std::span<const Particle>, std::string_view, int) {}

TEST(SimTests, PlottingTest) {
    GTEST_CXP_GCC auto result = std::invoke([] {
        std::array<Particle, 2> particles = {Particle{vec{}, vec{}, 1}, Particle{vec{1, 0, 0},
vec{}, 1}}; GTEST_CXP double delta_t = 0.014; run_simulation<{.force = gravitational_force, .io =
nothing_io}>( particles, {.delta_t = delta_t, .end_time = delta_t}, "unused"
        );
        return particles;
    });
    GCC_STATIC_EXPECT_EQ(result[0].f, vec(1, 0, 0));
}
*/
// NOLINTEND(*magic-numbers)
