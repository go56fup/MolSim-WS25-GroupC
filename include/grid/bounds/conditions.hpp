#pragma once

#include "grid/bounds/operations.hpp"
#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
#include "iterators/periodic.hpp"
#include "physics/forces.hpp"
#include "utility/constants.hpp"
#include "utility/tracing/macros.hpp"

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

template <boundary_type Boundary>
constexpr void reflect_via_ghost_particle(
	particle_container& container, const particle_container::index& cell_idx
) noexcept {
	using enum boundary_type;

	auto& current_cell = container[cell_idx];
	auto& system = container.system();
	const auto& domain = container.domain();

	auto reflect = [&]<boundary_type b>(const auto& make_ghost_pos) {
		for (auto p : current_cell) {
			const double ghost_particle_threshold = system.sigma[p] * sixth_root_of_2;
			const vec min{
				ghost_particle_threshold, ghost_particle_threshold, ghost_particle_threshold
			};
			const auto max = domain - min;
			if (out_of_bounds_soa<b>(system, p, max, min)) {
				const vec current_pos{system.x[p], system.y[p], system.z[p]};
				const auto resulting_force = lennard_jones_force(
					{.p1_position = current_pos,
				     .p2_position = make_ghost_pos(current_pos),
				     .sigma = system.sigma[p],
				     .epsilon = system.epsilon[p]}
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

	switch (Boundary) {
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

template <boundary_type Boundary>
constexpr void delete_ouflowing_particles(
	particle_container& container, const particle_container::index& cell_idx
) {
	using enum boundary_type;
	auto& cell = container[cell_idx];
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
		switch (Boundary) {
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

template <boundary_type Boundary>
constexpr void periodic_particle_interactions(
	particle_container& container, const particle_container::index& cell_idx
) {
	using enum boundary_type;
	using signed_ = particle_container::difference_type;

	particle_container::signed_index current_virtual_idx{};
	const auto& grid = container.grid_size();

	particle_container::signed_index signed_grid{
		static_cast<signed_>(grid.x), static_cast<signed_>(grid.y), static_cast<signed_>(grid.z)
	};

	const particle_container::signed_index signed_idx{
		static_cast<signed_>(cell_idx.x), static_cast<signed_>(cell_idx.y),
		static_cast<signed_>(cell_idx.z)
	};

	current_virtual_idx = signed_idx;

	vec particle_mirror{};
	const double cell_width = container.cutoff_radius();

	auto handle_boundary = [&]<boundary_type BoundaryChecked> {
		if constexpr ((Boundary & BoundaryChecked) != BoundaryChecked) {
			return;
		}
		static constexpr auto axis_ = boundary_type_to_axis(BoundaryChecked);
		static constexpr bool is_min = boundary_type_to_extremum(BoundaryChecked) == extremum::min;
		static constexpr int sign = is_min ? 1 : -1;

		particle_mirror[axis_] = signed_grid[axis_] * cell_width * sign;
		current_virtual_idx[axis_] = is_min ? signed_grid[axis_] : -1;
	};

	handle_boundary.template operator()<x_min>();
	handle_boundary.template operator()<x_max>();
	handle_boundary.template operator()<y_min>();
	handle_boundary.template operator()<y_max>();
	handle_boundary.template operator()<z_min>();
	handle_boundary.template operator()<z_max>();

	auto use_batch = [&](particle_batch batch_p1, particle_batch batch_p2) {
		lennard_jones_force_soa_batchwise(container, batch_p1, batch_p2, particle_mirror);
	};
	auto use_batch_piecewise = [&](particle_batch batch_p1, particle_batch batch_p2,
	                               std::size_t up_to) {
		for (std::size_t i = 0; i < up_to; ++i) {
			lennard_jones_force_soa(container, batch_p1[i], batch_p2[i], particle_mirror);
		}
	};

	for (const auto& periodic_target : periodic_range(container, current_virtual_idx)) {
		TRACE_PERIODIC(
			"Calculating periodic interactions for {} at {} via virtual index: {}", cell_idx,
			Boundary, current_virtual_idx
		);

#ifndef SINGLETHREADED
#pragma omp parallel for schedule(dynamic)
#endif
		for (auto current_p : container[cell_idx]) {
			std::size_t count = 0;
			particle_batch batch_p1{};
			particle_batch batch_p2{};

			for (auto periodic_p : container[periodic_target]) {
				batch_p1[count] = current_p;
				batch_p2[count] = periodic_p;
				++count;

				if (count == batch_size) {
					use_batch(batch_p1, batch_p2);
					count = 0;
				}
			}
			use_batch_piecewise(batch_p1, batch_p2, count);
		}
	}
}

namespace detail {
constexpr bool has_periodic(const sim_configuration& config) {
	using enum boundary_type;
	static const bool result = config.boundary_behavior[x_min] == boundary_condition::periodic ||
	                           config.boundary_behavior[y_min] == boundary_condition::periodic ||
	                           config.boundary_behavior[z_min] == boundary_condition::periodic;
	return result;
}

template <boundary_type... Boundaries>
constexpr bool is_periodic(const sim_configuration& config) {
	// TODO(tuna): see if making the result of this check static actually results in different
	// static bools across the templates
	return (... && (config.boundary_behavior[Boundaries] == boundary_condition::periodic));
}
}  // namespace detail

template <boundary_type Boundary>
constexpr void handle_boundary_condition(
	const sim_configuration& config, particle_container& container,
	const particle_container::index& cell_idx
) {
	using enum boundary_type;

	TRACE_BORDER_CELL_ITER("Checking boundary conditions for index: {}", cell_idx);
	const auto behavior = config.boundary_behavior[Boundary];
	switch (behavior) {
	case boundary_condition::reflecting:
		reflect_via_ghost_particle<Boundary>(container, cell_idx);
		break;
	case boundary_condition::outflow:
		delete_ouflowing_particles<Boundary>(container, cell_idx);
		break;
	case boundary_condition::periodic:
		periodic_particle_interactions<Boundary>(container, cell_idx);
		break;
	}
}

constexpr void handle_boundaries(particle_container& container, const sim_configuration& config) {
	const auto& grid = container.grid_size();
	using enum boundary_type;

	// TODO(gabriel):parallelis this (Attention corners need to be done atomically)
	for (particle_container::size_type i = 0; i < grid.y; ++i) {
		for (particle_container::size_type j = 0; j < grid.z; ++j) {
			handle_boundary_condition<x_min>(config, container, {0, i, j});
			handle_boundary_condition<x_max>(config, container, {grid.x - 1, i, j});
		}
	}

	for (particle_container::size_type i = 0; i < grid.x; ++i) {
		for (particle_container::size_type j = 0; j < grid.z; ++j) {
			handle_boundary_condition<y_min>(config, container, {i, 0, j});
			handle_boundary_condition<y_max>(config, container, {i, grid.y - 1, j});
		}
	}

	for (particle_container::size_type i = 0; i < grid.x; ++i) {
		for (particle_container::size_type j = 0; j < grid.y; ++j) {
			handle_boundary_condition<z_min>(config, container, {i, j, 0});
			handle_boundary_condition<z_max>(config, container, {i, j, grid.z - 1});
		}
	}

	if (!detail::has_periodic(config)) return;

	if (detail::is_periodic<x_min, y_min>(config)) {
		for (particle_container::size_type i = 0; i < grid.z; ++i) {
			periodic_particle_interactions<x_min | y_min>(container, {0, 0, i});
			periodic_particle_interactions<x_min | y_max>(container, {0, grid.y, i});
			periodic_particle_interactions<x_max | y_min>(container, {grid.x, 0, i});
			periodic_particle_interactions<x_max | y_max>(container, {grid.x, grid.y, i});
		}
	}

	if (detail::is_periodic<x_min, z_min>(config)) {
		for (particle_container::size_type i = 0; i < grid.y; ++i) {
			periodic_particle_interactions<x_min | z_min>(container, {0, i, 0});
			periodic_particle_interactions<x_min | z_max>(container, {0, i, grid.z});
			periodic_particle_interactions<x_max | z_min>(container, {grid.x, i, 0});
			periodic_particle_interactions<x_max | z_max>(container, {grid.x, i, grid.z});
		}
	}

	if (detail::is_periodic<y_min, z_min>(config)) {
		for (particle_container::size_type i = 0; i < grid.y; ++i) {
			periodic_particle_interactions<y_min | z_min>(container, {i, 0, 0});
			periodic_particle_interactions<y_min | z_max>(container, {i, 0, grid.z});
			periodic_particle_interactions<y_max | z_min>(container, {i, grid.y, 0});
			periodic_particle_interactions<y_max | z_max>(container, {i, grid.y, grid.z});
		}
	}

	if (detail::is_periodic<x_min, y_min, z_min>(config)) {
		periodic_particle_interactions<x_min | y_min | z_min>(container, {0, 0, 0});

		periodic_particle_interactions<x_max | y_min | z_min>(container, {grid.x, 0, 0});
		periodic_particle_interactions<x_min | y_max | z_min>(container, {0, grid.y, 0});
		periodic_particle_interactions<x_min | y_min | z_max>(container, {0, 0, grid.z});

		periodic_particle_interactions<x_max | y_max | z_min>(container, {grid.x, grid.y, 0});
		periodic_particle_interactions<x_min | y_max | z_max>(container, {0, grid.y, grid.z});
		periodic_particle_interactions<x_max | y_min | z_max>(container, {grid.x, 0, grid.z});

		periodic_particle_interactions<x_max | y_max | z_max>(container, {grid.x, grid.y, grid.z});
	}
}
