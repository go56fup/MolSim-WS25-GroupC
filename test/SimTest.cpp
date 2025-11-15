#include "Debug.h"
#include "gtest_constexpr.h"
#include "gtest/gtest.h"

#include <array>

#include "Debug.h"
#include "ForceCalculators.h"
#include "MolSim.h"
#include "Particle.h"

// TODO(tuna): uncomment

TEST(ForceTests, BasicGravitation) {
	GTEST_CXP auto ok = std::invoke([] {
		ParticleContainer particles(2, 2, 2, 1.5);
		particles.emplace(vec{}, vec{}, 1, 10);
		particles.emplace(vec{1, 0, 0}, vec{}, 1, 15);
		run_sim_iteration(gravitational_force, particles, 0.0014);
		auto moving_par = particles.view() | std::views::filter([](const Particle& p) { return p.type == 15; });
		return moving_par.begin()->f == vec{1, 0, 0};
	});
	GCC_STATIC_EXPECT_TRUE(ok);
}

/*
// Check that the Lennard-Jones force is calculated correctly
TEST(ForceTests, LennardJones) {
    GTEST_CXP_GCC auto particles_f = std::invoke([] {
        std::array<Particle, 2> particles = {Particle{vec{}, vec{}, 1}, Particle{vec{1, 2, 0}, vec{}, 1}};
        calculateF(lennard_jones_force, particles);
        return particles;
    });

    // NOLINTNEXTLINE(*magic-numbers)
    GTEST_CXP_GCC auto scaled_force = 15625 * particles_f[0].f;
    static constexpr auto result = vec(2952, 5904, 0);
    GCC_STATIC_EXPECT_VEC_DOUBLE_EQ(scaled_force, result);
}

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
    GCC_STATIC_EXPECT_DOUBLE_EQ(comet_initial.x, comet_after_one_revolution.x - deviation_x);
    GCC_STATIC_EXPECT_DOUBLE_EQ(comet_initial.y, comet_after_one_revolution.y - deviation_y);
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
