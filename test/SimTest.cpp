#include "Debug.h"
#include "FileReader.h"
#include "gtest_constexpr.h"
#include "ParticleContainer.h"
#include "gtest/gtest.h"

#include <array>

#include "Debug.h"
#include "ForceCalculators.h"
#include "MolSim.h"
#include "Particle.h"

#include "SharedTools.h"

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

// TODO(tuna): do not copy the particle
constexpr Particle get_particle_from_type(ParticleContainer& particles, decltype(Particle::type) type) {
	auto view = particles.view() | std::views::filter([&](const Particle& p) { return p.type == type; });
	return *view.begin();
}

#if 0
// Check that the Lennard-Jones force is calculated correctly
// TODO(tuna): make sim_configuration constexpr capable by migrating to fixed_string of 256,
// and then make the test actually constexpr
TEST(ForceTests, LennardJones) {
	GTEST_CXP Particle first(vec{4, 4, 4}, vec{}, 1, 10);
	GTEST_CXP Particle second(vec{5, 6, 4}, vec{}, 1, 15);
	auto config = sim_configuration{
		.delta_t = 0.0014,
		.cutoff_radius = 5,
		.sigma = 1,
		.epsilon = 5,
		.boundary_behavior = boundary_conditions_descriptor::all(boundary_condition::reflecting),
		.end_time = 0.0014,
		.write_frequency = 1000,
		.base_name = "unused",
		.domain{10, 10, 10},
		.meshwidth = 1.125
	};
	auto calc = [sigma = config.sigma, eps = config.epsilon](const Particle& p1, const Particle& p2) noexcept {
		return lennard_jones_force(p1, p2, sigma, eps);
	};
	auto first_force = std::invoke([&] -> vec {
		ParticleContainer particles(config.domain, config.cutoff_radius);
		particles.place(first);
		particles.place(second);
		run_sim_iteration(calc, particles, config);
		return get_particle_from_type(particles, 10).f;
	});

	auto direct_calc = lennard_jones_force(first, second, config.sigma, config.epsilon);
	EXPECT_VEC_DOUBLE_EQ(first_force, direct_calc);
	static constexpr auto hand_calculated_result = vec(2952, 5904, 0);
	auto scaled_force = 15625 * first_force;
	EXPECT_VEC_DOUBLE_EQ(hand_calculated_result, scaled_force);
}

TEST(ForceTests, GridTest) {
	Particle first{vec{2, 1, 1}, vec{1, 0, 0}, 1, 10};
	Particle second{vec{5, 1, 1}, vec{-1, 0, 0}, 1, 15};
		auto config = sim_configuration{
		.delta_t = 0.001,
		.cutoff_radius = 2,
		.sigma = 1,
		.epsilon = 5,
		.boundary_behavior = boundary_conditions_descriptor::all(boundary_condition::outflow),
		.end_time = 2,
		.write_frequency = 1000,
		.base_name = "unused",
		.domain{21, 21, 21},
		.meshwidth = 1.125
	};
	auto calc = [sigma = config.sigma, eps = config.epsilon](const Particle& p1, const Particle& p2) noexcept {
		return lennard_jones_force(p1, p2, sigma, eps);
	};
	auto particle_result = std::invoke([&] -> std::pair<Particle, Particle> {
		ParticleContainer particles(config.domain, config.cutoff_radius);
		particles.place(first);
		particles.place(second);
		run_simulation(particles, config, calc, "/tmp/paraview");
		return {get_particle_from_type(particles, 10), get_particle_from_type(particles, 15)};
	});
	SPDLOG_DEBUG("got particles at the end of sim: {}\n{}", particle_result.first, particle_result.second);
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

    // "During its 1986 appearance, Halley's nearest approach to Earth occurred (...) at a distance of
    // 0.42 AU (...)"
    // source: https://science.nasa.gov/solar-system/comets/1p-halley/
    GCC_STATIC_EXPECT_VEC_NEAR(earth_coords, comet_coords, 0.16);

    GTEST_CXP auto until_return_to_initial_comet_pos = 33690;
    GTEST_CXP const vec& comet_initial = comet.x;
    GTEST_CXP_GCC vec comet_after_one_revolution = run_sim(until_return_to_initial_comet_pos)[comet_i].x;

    GTEST_CXP double deviation_x = 1.4064506667406391;
    GTEST_CXP double deviation_y = -0.004078428002306117;
    GCC_STATIC_EXPECT_NEAR(comet_initial.x, comet_after_one_revolution.x - deviation_x, 0.001);
    GCC_STATIC_EXPECT_NEAR(comet_initial.y, comet_after_one_revolution.y - deviation_y, 0.001);
}

constexpr void nothing_io(std::span<const Particle>, std::string_view, int) {}

TEST(SimTests, PlottingTest) {
    GTEST_CXP_GCC auto result = std::invoke([] {
        std::array<Particle, 2> particles = {Particle{vec{}, vec{}, 1}, Particle{vec{1, 0, 0}, vec{}, 1}};
        GTEST_CXP double delta_t = 0.014;
        run_simulation<{.force = gravitational_force, .io = nothing_io}>(
            particles, {.delta_t = delta_t, .end_time = delta_t}, "unused"
        );
        return particles;
    });
    GCC_STATIC_EXPECT_EQ(result[0].f, vec(1, 0, 0));
}
*/
// NOLINTEND(*magic-numbers)
