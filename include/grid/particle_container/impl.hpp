#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <vector>

#include "grid/bounds/operations.hpp"
#include "grid/enums.hpp"
#include "iterators/border_cells.hpp"
#include "iterators/enumerate_cells.hpp"
#include "iterators/interactions.hpp"
#include "physics/maxwell_boltzmann.hpp"
#include "physics/particle.hpp"
#include "physics/vec_3d.hpp"
#include "simulation/config/entities.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/concepts.hpp"
#include "utility/tracing/config.hpp"
#include "utility/tracing/macros.hpp"

constexpr particle_container::size_type particle_container::linear_index(
	particle_container::size_type x, particle_container::size_type y,
	particle_container::size_type z
) const {
	const auto result = (x * grid_size_.y * grid_size_.z) + (y * grid_size_.z) + z;
#if LOG_GRID
	if (result >= grid_size_.x * grid_size_.y * grid_size_.z) {
		throw std::out_of_range(
			fmt::format("Out of bounds cell requested for grid {}: {}, {}, {}", grid_size_, x, y, z)
		);
	}
#endif
	return result;
}

constexpr particle_container::index particle_container::index_from_linear(size_type linear) const {
	const size_type yz = grid_size_.y * grid_size_.z;

	// TODO(tuna): move this check to its own macro that is empty when !LOG_GRID
#if LOG_GRID
	if (linear >= grid_size_.x * yz) {
		throw std::out_of_range(
			fmt::format("Out of bounds linear index {} for grid {}", linear, grid_size_)
		);
	}
#endif

	const size_type x = linear / yz;
	const size_type rem = linear % yz;
	const size_type y = rem / grid_size_.z;
	const size_type z = rem % grid_size_.z;

	return {x, y, z};
}

constexpr particle_container::size_type particle_container::pos_to_linear_index(const vec& pos
) const noexcept {
	const auto x = static_cast<size_type>(pos.x / cutoff_radius_);
	const auto y = static_cast<size_type>(pos.y / cutoff_radius_);
	const auto z = static_cast<size_type>(pos.z / cutoff_radius_);
	return linear_index(x, y, z);
}

namespace detail {
/**
 * @brief Divides two doubles and rounds the result up to @ref size_type.
 * @param dividend Dividend.
 * @param divisor Divisor.
 * @return \f[ \left\lceil \frac{dividend}{divisor} \right\rceil \f]
 */

CONSTEXPR_IF_GCC inline particle_container::size_type
div_round_up(double dividend, double divisor) noexcept {
	return static_cast<particle_container::size_type>(std::ceil(dividend / divisor));
}
}  // namespace detail

constexpr void particle_container::check_if_out_of_domain_max(const vec& pos) const
	noexcept(false) {
	using enum boundary_type;
	if (out_of_bounds<x_max>(pos, domain_) || out_of_bounds<y_max>(pos, domain_) ||
	    out_of_bounds<z_max>(pos, domain_)) {
		throw std::domain_error(fmt::format(
			"Refusing to place particle at {}, which is on or outside the domain of the "
			"simulation: {}",
			pos, domain_
		));
	}
}

template <fwd_reference_to<vec> Vec>
constexpr particle_container::particle_container(Vec&& domain_arg, double cutoff_radius_arg)
		: domain_{std::forward<Vec>(domain_arg)}
		, grid_size_{
				detail::div_round_up(domain_.x, cutoff_radius_arg),
				detail::div_round_up(domain_.y, cutoff_radius_arg),
				detail::div_round_up(domain_.z, cutoff_radius_arg)}
		, cutoff_radius_{cutoff_radius_arg} {
	// Prevent implicit widening later, if x * y * z overflows for uint32, this takes care of it
	const auto grid_x = static_cast<std::size_t>(grid_size_.x);
	const auto grid_y = static_cast<std::size_t>(grid_size_.y);
	const auto grid_z = static_cast<std::size_t>(grid_size_.z);
	TRACE_PARTICLE_CONTAINER("Setting up grid with size: {}", grid_size_);
	grid.resize(grid_x * grid_y * grid_z);
}

constexpr void particle_container::add_particle(
	const vec& position, const vec& velocity, particle_system::particle_type_t type
) {
	system_.add_particle(position, velocity, type);
	grid[pos_to_linear_index(position)].push_back(system_.size() - 1);
}

template <fwd_reference_to<material_description> Desc>
constexpr std::uint8_t particle_container::register_material(Desc&& desc) noexcept {
	auto it = std::ranges::find(materials_, desc);
	if (it == materials_.end()) {
		++particle_type_count;
		assert(particle_type_count < materials_.size());
		materials_[particle_type_count] = std::forward<Desc>(desc);
		return particle_type_count;
	}
	return static_cast<std::uint8_t>(std::distance(materials_.begin(), it));
}

constexpr void particle_container::reload_particle_state(
	const particle_state_parameters& parameters, const vec& velocity,
	particle_system::particle_type_t type
) {
	system_.add_particle(
		parameters.position, velocity, parameters.force, parameters.old_force, type
	);
	grid[pos_to_linear_index(parameters.position)].push_back(system_.size() - 1);
}

template <std::size_t N>
constexpr void particle_container::add_cuboid(
	const cuboid_parameters<3>& cuboid, const body_common_parameters& body_parameters,
	const vec& velocity, const material_description& material, std::size_t& seq_no
) {
	const auto type = register_material(material);
	for (size_type x = 0; x < cuboid.scale.x; ++x) {
		for (size_type y = 0; y < cuboid.scale.y; ++y) {
			for (size_type z = 0; z < cuboid.scale.z; ++z) {
				add_particle(
					{cuboid.origin.x + (x * body_parameters.meshwidth),
				     cuboid.origin.y + (y * body_parameters.meshwidth),
				     cuboid.origin.z + (z * body_parameters.meshwidth)},
					maxwell_boltzmann_distributed_velocity<N>(
						body_parameters.brownian_mean, seq_no
					) + velocity,
					type
				);
			}
		}
	}
}

constexpr void particle_container::add_disc(
	const disc_parameters<3>& disc, const body_common_parameters& body_parameters,
	const vec& velocity, const material_description& material, std::size_t& seq_no
) {
	const auto& r = disc.radius;
	const auto type = register_material(material);
	for (int i = static_cast<int>(-r); i <= r; i++) {
		for (int j = static_cast<int>(-r); j <= r; j++) {
			if (i * i + j * j > r * r) {
				continue;
			}
			add_particle(
				{disc.center.x + (i * body_parameters.meshwidth),
			     disc.center.y + (j * body_parameters.meshwidth), disc.center.z},
				maxwell_boltzmann_distributed_velocity<2>(body_parameters.brownian_mean, seq_no) +
					velocity,
				type
			);
		}
	}
}

// TODO(tuna): see if the return type specified conflicts when the container is const
constexpr range_of<
	std::pair<const particle_container::index&, const particle_container::index&>> auto
particle_container::directional_interactions() noexcept {
	return interactions_range(*this);
}

// TODO(tuna): change this once the return type of the new border cell is fixed
constexpr /* range_of<border_cell_range::iterator::value_type> */ auto
particle_container::border_cells() noexcept {
	return border_cell_range(*this);
}

// The return type of this contains references to cells, which are conditionally const, which
// cannot be expressed in here: so we leave it as auto.
constexpr auto particle_container::enumerate_cells() noexcept {
	return enumerate_cells_range(*this);
}

constexpr range_of<particle_container::cell> auto& particle_container::cells() noexcept {
	return grid;
}

constexpr const range_of<const particle_container::cell> auto&
particle_container::cells() const noexcept {
	return grid;
}

constexpr particle_container::cell&
particle_container::operator[](const particle_container::index& three_d_index) noexcept {
	return grid[linear_index(three_d_index.x, three_d_index.y, three_d_index.z)];
}

constexpr const particle_container::cell&
particle_container::operator[](const particle_container::index& three_d_index) const noexcept {
	return grid[linear_index(three_d_index.x, three_d_index.y, three_d_index.z)];
}

constexpr particle_container::cell& particle_container::operator[](
	particle_container::size_type x, particle_container::size_type y,
	particle_container::size_type z
) noexcept {
	return grid[linear_index(x, y, z)];
}

constexpr const particle_container::cell& particle_container::operator[](
	particle_container::size_type x, particle_container::size_type y,
	particle_container::size_type z
) const noexcept {
	return grid[linear_index(x, y, z)];
}

constexpr particle_container::cell& particle_container::cell_containing(const vec& position
) noexcept {
	return grid[pos_to_linear_index(position)];
}

constexpr const particle_container::cell& particle_container::cell_containing(const vec& pos
) const noexcept {
	return grid[pos_to_linear_index(pos)];
}

constexpr const vec& particle_container::domain() const noexcept {
	return domain_;
}

constexpr const particle_container::index& particle_container::grid_size() const noexcept {
	return grid_size_;
}

constexpr double particle_container::cutoff_radius() const noexcept {
	return cutoff_radius_;
}

constexpr const material_description&
particle_container::material_for_particle(particle_system::particle_id i) const noexcept {
	return materials_[system_.type[i]];
}

constexpr const particle_system& particle_container::system() const noexcept {
	return system_;
}

constexpr particle_system& particle_container::system() noexcept {
	return system_;
}
