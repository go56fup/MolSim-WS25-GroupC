#pragma once

#include "grid/bounds/operations.hpp"
#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
#include "iterators/periodic.hpp"
#include "physics/forces.hpp"
#include "utility/constants.hpp"

namespace detail {
constexpr vec x_min_f(const vec& x) noexcept {
	return {-x.x, x.y, x.z};
}

constexpr vec y_min_f(const vec& x) noexcept {
	return {x.x, -x.y, x.z};
}

constexpr vec z_min_f(const vec& x) {
	return {x.x, x.y, -x.z};
}
}  // namespace detail

constexpr void reflect_via_ghost_particle(
	particle_container::cell& current_cell, particle_container& container, boundary_type border_type
) noexcept {
	using enum boundary_type;

	auto& system = container.system();
	const auto& domain = container.domain();

	auto reflect = [&]<boundary_type b>(const auto& make_ghost_pos) {
		for (auto p : current_cell) {
			const auto& material = container.material_for_particle(p);
			const double ghost_particle_threshold = material.sigma * sixth_root_of_2;
			const vec min{
				ghost_particle_threshold, ghost_particle_threshold, ghost_particle_threshold
			};
			const auto max = domain - min;
			const bool check = out_of_bounds_soa<b>(system, p, max, min);

			if (check) {
				const vec current_pos{system.x[p], system.y[p], system.z[p]};
				const auto& material = container.material_for_particle(p);
				const auto resulting_force = lennard_jones_force(
					{.p1_position = current_pos,
				     .p2_position = make_ghost_pos(current_pos),
				     .sigma = material.sigma,
				     .epsilon = material.epsilon}
				);
				system.fx[p] += resulting_force.x;
				system.fy[p] += resulting_force.y;
				system.fz[p] += resulting_force.z;
			}
		}
	};

	const auto x_max_f = [&](const vec& x) { return vec{-x.x + 2 * domain.x, x.y, x.z}; };
	const auto y_max_f = [&](const vec& x) { return vec{x.x, -x.y + 2 * domain.y, x.z}; };
	const auto z_max_f = [&](const vec& x) { return vec{x.x, x.y, -x.z + 2 * domain.z}; };

	switch (border_type) {
	case x_min:
		reflect.template operator()<x_min>(detail::x_min_f);
		return;

	case y_min:
		reflect.template operator()<y_min>(detail::y_min_f);
		return;

	case z_min:
		reflect.template operator()<z_min>(detail::z_min_f);
		return;

	case x_max:
		reflect.template operator()<x_max>(x_max_f);
		return;

	case y_max:
		reflect.template operator()<y_max>(y_max_f);
		return;

	case z_max:
		reflect.template operator()<z_max>(z_max_f);
		return;
	}
}

constexpr void delete_ouflowing_particles(
	particle_container::cell& cell, particle_container& container, boundary_type border
) {
	using enum boundary_type;
	auto& system = container.system();
	const auto& domain = container.domain();

	auto erase_if = [&]<boundary_type b>(std::size_t cell_idx) {
		if (out_of_bounds_soa<b>(system, cell[cell_idx], domain)) {
			cell.erase(std::next(
				cell.begin(), static_cast<particle_container::cell::difference_type>(cell_idx)
			));
		}
	};
	for (std::size_t cell_idx = 0; cell_idx < cell.size(); ++cell_idx) {
		switch (border) {
		case x_min:
			erase_if.template operator()<x_min>(cell_idx);
			break;
		case x_max:
			erase_if.template operator()<x_max>(cell_idx);
			break;
		case y_min:
			erase_if.template operator()<y_min>(cell_idx);
			break;
		case y_max:
			erase_if.template operator()<y_max>(cell_idx);
			break;
		case z_min:
			erase_if.template operator()<z_min>(cell_idx);
			break;
		case z_max:
			erase_if.template operator()<z_max>(cell_idx);
			break;
		}
	}
}

constexpr void periodic(
	particle_container::cell& current_cell, boundary_type border_type,
	force_calculator auto calculator, const particle_container::index& idx,
	particle_container& container
) {

	using enum boundary_type;
	using difference_type = particle_container::difference_type;
	using signed_index = vec_3d<difference_type>;
	signed_index current_virtual_idx{};

	const auto& grid = container.grid_size();
	signed_index signed_grid = {
		static_cast<difference_type>(grid.x), static_cast<difference_type>(grid.y),
		static_cast<difference_type>(grid.z)
	};
	const vec_3d<difference_type> signed_idx = {
		static_cast<difference_type>(idx.x), static_cast<difference_type>(idx.y),
		static_cast<difference_type>(idx.z)
	};

	current_virtual_idx = signed_idx;
	vec particle_mirror{};
	const double cell_width = container.cutoff_radius();

	// TODO(tuna): fold into axis and lambda
	if ((border_type & x_min) == x_min) {
		particle_mirror.x = signed_grid.x * cell_width;
		current_virtual_idx.x = signed_grid.x;
	}
	if ((border_type & x_max) == x_max) {
		particle_mirror.x = -signed_grid.x * cell_width;
		current_virtual_idx.x = -1;
	}
	if ((border_type & y_min) == y_min) {
		particle_mirror.y = signed_grid.y * cell_width;
		current_virtual_idx.y = signed_grid.y;
	}
	if ((border_type & y_max) == y_max) {
		particle_mirror.y = -signed_grid.y * cell_width;
		current_virtual_idx.y = -1;
	}
	if ((border_type & z_min) == z_min) {
		particle_mirror.z = signed_grid.z * cell_width;
		current_virtual_idx.z = signed_grid.z;
	}
	if ((border_type & z_max) == z_max) {
		particle_mirror.z = -signed_grid.z * cell_width;
		current_virtual_idx.z = -1;
	}
	auto& system = container.system();
	TRACE_GRID("Calculating periodic interactions via virtual index: {}", current_virtual_idx);
	for (const auto& periodic_target : periodic_range(container, current_virtual_idx)) {
		for (auto current_p : current_cell) {
			for (auto periodic_p : container[periodic_target]) {
				// TODO(tuna): is this actually the best way of doing this?
				system.x[current_p] += particle_mirror.x;
				system.y[current_p] += particle_mirror.y;
				system.z[current_p] += particle_mirror.z;
				std::invoke(calculator, container, current_p, periodic_p);
				system.x[current_p] -= particle_mirror.x;
				system.y[current_p] -= particle_mirror.y;
				system.z[current_p] -= particle_mirror.z;
			}
		}
	}
}
