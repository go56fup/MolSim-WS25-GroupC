#include "Debug.h"
#include "gtest_constexpr.h"
#include "gtest/gtest.h"

#include <array>
#include <cstdint>
#include <functional>

#include "ForceCalculators.h"
#include "MolSim.h"
#include "Particle.h"

template class vec_3d<double>;
template class flag_special_member_funcs<"Particle">;
template class flag_special_member_funcs<"Vector">;

#define INTERNAL_VEC_EQ_AT(a, b, i, abs_error, macro) macro((a).i, (b).i, (abs_error))
#define INTERNAL_VEC_NEAR(a, b, abs_error, macro)                                                                      \
	INTERNAL_VEC_EQ_AT(a, b, x, abs_error, macro);                                                                     \
	INTERNAL_VEC_EQ_AT(a, b, y, abs_error, macro);                                                                     \
	INTERNAL_VEC_EQ_AT(a, b, z, abs_error, macro)

#define GCC_STATIC_EXPECT_VEC_NEAR(a, b, abs_error) INTERNAL_VEC_NEAR(a, b, abs_error, GCC_STATIC_EXPECT_NEAR)
#define EXPECT_VEC_NEAR(a, b, abs_error) INTERNAL_VEC_NEAR(a, b, abs_error, EXPECT_NEAR)

TEST(ForceTests, BasicGravitation) {
	static GTEST_CXP_GCC auto particles_f = std::invoke([] {
		std::array<Particle, 2> particles = {{{vec{}, vec{}, 1}, {vec{1, 0, 0}, vec{}, 1}}};
		calculateF(gravitational_force, particles);
		return particles;
	});

	GCC_STATIC_EXPECT_EQ(particles_f[0].f, vec(1, 0, 0));
}

TEST(SimTests, HalleysComet) {
	static constexpr Particle sun{vec{}, vec{}, 1.0};
	static constexpr Particle earth{vec{0., 1., 0.}, vec{-1., 0., 0.}, 3.0e-6};
	static constexpr Particle jupiter{vec{0., 5.36, 0.}, vec{-0.45, 0., 0.}, 9.55e-4};
	static constexpr Particle comet{vec{34.75, 0., 0.}, vec{0., 0.0296, 0.}, 1e-14};

	enum planet_ids : std::uint8_t { sun_i = 0, earth_i = 1, jupiter_i = 2, comet_i = 3 };

	static GTEST_CXP_GCC auto run_sim = [](std::size_t iteration) {
		std::array<Particle, 4> planets{sun, earth, jupiter, comet};
		static constexpr double delta_t = 0.014;

		for (std::size_t i = 0; i < iteration; ++i) {
			calculateX(planets, delta_t);
			calculateF(gravitational_force, planets);
			calculateV(planets, delta_t);
			update_values(planets);
		}

		// run_simulation(planets, 0.0014, iteration);
		return planets;
	};
	static constexpr auto until_earth_sighting = 16600;
	static GTEST_CXP_GCC auto earth_sighting = run_sim(until_earth_sighting);
	static GTEST_CXP_GCC const vec& earth_coords = earth_sighting[earth_i].x;
	static GTEST_CXP_GCC const vec& comet_coords = earth_sighting[comet_i].x;

	// "During its 1986 appearance, Halley's nearest approach to Earth occurred (...) at a distance of
	// 0.42 AU (...)"
	// source: https://science.nasa.gov/solar-system/comets/1p-halley/
	GCC_STATIC_EXPECT_VEC_NEAR(earth_coords, comet_coords, 0.16);

	static constexpr auto until_return_to_initial_comet_pos = 33690;
	static constexpr const vec& comet_initial = comet.x;
	static GTEST_CXP_GCC vec comet_after_one_revolution = run_sim(until_return_to_initial_comet_pos)[comet_i].x;

	static constexpr double deviation_x = 1.4064506667406391;
	static constexpr double deviation_y = -0.004078428002306117;
	GCC_STATIC_EXPECT_DOUBLE_EQ(comet_initial.x, comet_after_one_revolution.x - deviation_x);
	GCC_STATIC_EXPECT_DOUBLE_EQ(comet_initial.y, comet_after_one_revolution.y - deviation_y);
}

// NOLINTBEGIN(*magic-numbers)

TEST(VectorTests, Addition) {
	static constexpr vec a{1, 2, 3};
	static constexpr vec b{4, 5, 6};
	static constexpr vec result{5, 7, 9};
	STATIC_EXPECT_EQ(a + b, result);
}

TEST(VectorTests, Subtraction) {
	static constexpr vec a{5, 7, 9};
	static constexpr vec b{4, 5, 6};
	static constexpr vec result{1, 2, 3};
	STATIC_EXPECT_EQ(a - b, result);
}

TEST(VectorTests, AdditionAssignment) {
	static GTEST_CXP vec result = std::invoke([] {
		vec a{1, 2, 3};
		vec b{4, 5, 6};
		a += b;
		return a;
	});
	static constexpr vec expected{5, 7, 9};
	STATIC_EXPECT_EQ(result, expected);
}

TEST(VectorTests, SubtractionAssignment) {
	static GTEST_CXP vec result = std::invoke([] {
		vec a{4, 5, 6};
		vec b{1, 2, 3};
		a -= b;
		return a;
	});
	static constexpr vec expected{3, 3, 3};
	STATIC_EXPECT_EQ(result, expected);
}

TEST(VectorTests, Scalar) {
	static constexpr vec a{5, 7, 9};
	static constexpr vec result{10, 14, 18};
	STATIC_EXPECT_EQ(a * 2, result);
	STATIC_EXPECT_EQ(2 * a, result);
}

TEST(VectorTests, ScalarAssignment) {
	static GTEST_CXP vec result = std::invoke([] {
		vec a{5, 7, 9};
		a *= 2;
		return a;
	});

	static constexpr vec expected{10, 14, 18};
	STATIC_EXPECT_EQ(result, expected);
}

TEST(VectorTests, Norm) {
	static constexpr vec a{2, 3, 6};
	static constexpr double result = 7;
	GCC_STATIC_EXPECT_EQ(a.euclidian_norm(), result);
}

TEST(VectorTests, Equal) {
	static constexpr vec a{1, 2, 3};
	static constexpr vec b{1, 2, 3};
	STATIC_EXPECT_EQ(a, b);
}

// NOLINTEND(*magic-numbers)
