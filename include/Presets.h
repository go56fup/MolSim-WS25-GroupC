#pragma once

#include <string_view>

#include "MolSim.h"
#include "Vector.h"

namespace presets {
template <sim_traits Traits = {}>
void assignment2(std::string_view filename) {
	static constexpr vec x1{0, 0, 0};
	static constexpr vec v1{0, 0, 0};
	static constexpr vec_3d<int> N1{40, 8, 1};

	static constexpr vec x2{15, 15, 0};
	static constexpr vec v2{0, -10, 0};
	static constexpr vec_3d<int> N2{8, 8, 1};

	static constexpr double h = 1.1225;
	static constexpr double m = 1;
	static constexpr double average_brownian = 0.1;

	ParticleContainer particles;
	particles.reserve((N1.x * N1.y * N1.z) + (N2.x * N2.y * N2.z));
	cuboid_generator(particles, x1, N1, h, v1, m, average_brownian);
	cuboid_generator(particles, x2, N2, h, v2, m, average_brownian);
	// NOLINTNEXTLINE(*magic-numbers)
	run_simulation<Traits>(particles, {.delta_t = 0.0002, .end_time = 5}, filename);
}
}  // namespace presets
