#pragma once

#include "grid/bounds/operations.hpp"
#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
#include "physics/forces.hpp"
#include "simulation/config/entities.hpp"
#include "iterators/periodic.hpp"

// The wrapping in do-while increases this metric dramatically. Otherwise, the rest of the logic is
// impossible to express in more granular functions without having to recompute information. The
// template-for facility from C++26 would have come in useful here.
// NOLINTBEGIN(*cognitive-complexity)
// NOLINTBEGIN(*avoid-do-while)

namespace reflect {
constexpr void macro(
	particle_container::cell& current_cell, const vec& domain, boundary_type border_type,
	force_calculator auto calculator
) noexcept {
#define REFLECT_IF(check, x_mod, y_mod, z_mod)                                                     \
	do {                                                                                           \
		for (auto& p : current_cell) {                                                             \
			const double ghost_particle_threshold = p.sigma * sixth_root_of_2;                     \
			if (p.x.check) {                                                                       \
				const particle ghost{                                                              \
					{p.x.x_mod, p.x.y_mod, p.x.z_mod}, {}, p.m, p.sigma, p.epsilon                 \
				};                                                                                 \
				TRACE_GRID("Reflecting {} with ghost: {}", p, ghost);                              \
				p.f += std::invoke(calculator, p, ghost);                                          \
			}                                                                                      \
		}                                                                                          \
	} while (0)

	using enum boundary_type;

	switch (border_type) {
	case x_min:
		REFLECT_IF(x <= ghost_particle_threshold, x * -1, y, z);
		return;
	case y_min:
		REFLECT_IF(y <= ghost_particle_threshold, x, y * -1, z);
		return;
	case z_min:
		REFLECT_IF(z <= ghost_particle_threshold, x, y, z * -1);
		return;
	case x_max:
		REFLECT_IF(x >= domain.x - ghost_particle_threshold, x * -1 + (2 * domain.x), y, z);
		return;
	case y_max:
		REFLECT_IF(y >= domain.y - ghost_particle_threshold, x, y * -1 + (2 * domain.y), z);
		return;
	case z_max:
		REFLECT_IF(z >= domain.z - ghost_particle_threshold, x, y, z * -1 + (2 * domain.z));
		return;
	}
	assert(!"This should never be called with a non-border, the iterator is wrong.");
#undef REFLECT_IF
}

// NOLINTEND(*avoid-do-while)
// NOLINTEND(*cognitive-complexity)

constexpr void lambda(
	std::vector<particle>& current_cell, const vec& domain, boundary_type border_type,
	force_calculator auto calculator
) noexcept {
	using enum boundary_type;

	auto reflect = [&]<boundary_type b>(const auto& make_ghost_pos) {
		for (auto& p : current_cell) {
			const double ghost_particle_threshold = p.sigma * sixth_root_of_2;
			const vec min{
				ghost_particle_threshold, ghost_particle_threshold, ghost_particle_threshold
			};
			const auto max = domain - min;
			const bool check = out_of_bounds<b>(p.x, max, min);
			if (check) {
				const particle ghost{make_ghost_pos(p.x), {}, p.m};
				p.f += std::invoke(calculator, p, ghost);
			}
		}
	};

	switch (border_type) {
	case x_min:
		reflect.template operator()<x_min>([](const vec& x) static { return vec{-x.x, x.y, x.z}; });
		return;

	case y_min:
		reflect.template operator()<y_min>([](const vec& x) static { return vec{x.x, -x.y, x.z}; });
		return;

	case z_min:
		reflect.template operator()<z_min>([](const vec& x) static { return vec{x.x, x.y, -x.z}; });
		return;

	case x_max:
		reflect.template operator(
		)<x_max>([&](const vec& x) { return vec{-x.x + 2 * domain.x, x.y, x.z}; });
		return;

	case y_max:
		reflect.template operator(
		)<y_max>([&](const vec& x) { return vec{x.x, -x.y + 2 * domain.y, x.z}; });
		return;

	case z_max:
		reflect.template operator(
		)<z_max>([&](const vec& x) { return vec{x.x, x.y, -x.z + 2 * domain.z}; });
		return;
	}
}

}  // namespace reflect

constexpr void delete_ouflowing_particles(
	particle_container::cell& cell, const vec& domain, boundary_type border
) {
	using enum boundary_type;
	auto erase_if = [&](bool oob, std::size_t i) {
		if (oob)
			cell.erase(
				std::next(cell.begin(), static_cast<particle_container::cell::difference_type>(i))
			);
	};
	for (std::size_t i = 0; i < cell.size(); ++i) {
		const auto& particle_coords = cell[i].x;
		switch (border) {
		case x_min:
			erase_if(out_of_bounds<x_min>(particle_coords, domain), i);
			break;
		case x_max:
			erase_if(out_of_bounds<x_max>(particle_coords, domain), i);
			break;
		case y_min:
			erase_if(out_of_bounds<y_min>(particle_coords, domain), i);
			break;
		case y_max:
			erase_if(out_of_bounds<y_max>(particle_coords, domain), i);
			break;
		case z_min:
			erase_if(out_of_bounds<z_min>(particle_coords, domain), i);
			break;
		case z_max:
			erase_if(out_of_bounds<z_max>(particle_coords, domain), i);
			break;
		}
	}
}

#define reflect_via_ghost_particle reflect::macro

constexpr void periodic(
		std::vector<particle>& current_cell, boundary_type border_type,
		force_calculator auto calculator, const particle_container::index& idx,
		particle_container particles
		){

	using enum boundary_type;
	using difference_type = particle_container::difference_type;
	using signed_index = vec_3d<difference_type>;
	signed_index current_virtual_idx{};

	const auto& grid = particles.grid_size();
	signed_index signed_grid ={
			static_cast<difference_type>(grid.x),
			static_cast<difference_type>(grid.y),
			static_cast<difference_type>(grid.z)
	};
	vec_3d<difference_type> const signed_idx = {static_cast<difference_type>(idx.x), static_cast<difference_type>(idx.y),
							   static_cast<difference_type>(idx.z)};

	current_virtual_idx = signed_idx;

	if ((border_type & x_min) == x_min) {
		current_virtual_idx.x = signed_grid.x;
	} if ((border_type & x_max) == x_max) {
		current_virtual_idx.x = -1;
	} if ((border_type & y_min) == y_min) {
		current_virtual_idx.y = signed_grid.y;
	} if ((border_type & y_max) == y_max) {
		current_virtual_idx.y = -1;
	} if ((border_type & z_min) == z_min) {
		current_virtual_idx.z = signed_grid.z;
	} if ((border_type & z_max) == z_max) {
		current_virtual_idx.z = -1;
	}

	for (const auto& i : periodic_range(particles, current_virtual_idx)) {

		for (auto& p1 : current_cell) {
			for (auto& p2 : particles[i]) {
				apply_force_interaction(calculator, p1, p2);
			}
		}
	}
}


constexpr void handle_boundary_condition(
	particle_container::cell& cell, boundary_type border, const vec& domain,
	force_calculator auto force_calc, const sim_configuration& config
) noexcept {
	using enum boundary_type;
	auto exercise_boundary_condition = [&](boundary_type border) {
		switch (config.boundary_behavior[border]) {
		case boundary_condition::reflecting:
			reflect_via_ghost_particle(cell, domain, border, force_calc);
			break;
		case boundary_condition::outflow:
			delete_ouflowing_particles(cell, domain, border);
			break;
		case boundary_condition::periodic:
			// TODO(gabriel): implement this
			throw std::logic_error("Not implemented");
		}
	};

	for (auto [min, max] : border_pairs) {
		const bool on_min = (border & min) == min;
		const bool on_max = (border & max) == max;
		if (on_min) {
			exercise_boundary_condition(min);
		}
		if (on_max) {
			exercise_boundary_condition(max);
		}
	}
}
