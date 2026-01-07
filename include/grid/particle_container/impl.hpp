#pragma once
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <vector>

#include "fwd.hpp"
#include "grid/bounds/operations.hpp"
#include "grid/enums.hpp"
#include "iterators/border_cells.hpp"
#include "iterators/enumerate_cells.hpp"
#include "iterators/interactions.hpp"
#include "physics/maxwell_boltzmann.hpp"
#include "physics/particle.hpp"
#include "physics/vec_3d.hpp"
#include "utility/compiler_traits.hpp"
#include "utility/concepts.hpp"
#include "utility/tracing/config.hpp"

constexpr particle_container::size_type particle_container::linear_index(
	particle_container::size_type x, particle_container::size_type y,
	particle_container::size_type z
) const noexcept {
	const auto result = (x * grid_size_.y * grid_size_.z) + (y * grid_size_.z) + z;
	assert(result < grid_size_.x * grid_size_.y * grid_size_.z && "Out of bounds cell requested");
	return result;
}

constexpr particle_container::size_type
particle_container::pos_to_linear_index(const vec& pos) const noexcept {
	const auto x = static_cast<size_type>(pos.x / cutoff_radius_);
	const auto y = static_cast<size_type>(pos.y / cutoff_radius_);
	const auto z = static_cast<size_type>(pos.z / cutoff_radius_);
	return linear_index(x, y, z);
}

CONSTEXPR_IF_GCC particle_container::size_type
particle_container::div_round_up(double dividend, double divisor) noexcept {
	return static_cast<particle_container::size_type>(std::ceil(dividend / divisor));
}

constexpr void
particle_container::check_if_out_of_domain_max(const vec& pos) const noexcept(false) {
	using enum boundary_type;
	if (out_of_bounds<x_max>(pos, domain_) || out_of_bounds<y_max>(pos, domain_) ||
	    out_of_bounds<z_max>(pos, domain_)) {
		throw std::domain_error(
			fmt::format(
				"Refusing to place particle at {}, which is on or outside the domain of the "
				"simulation: {}",
				pos, domain_
			)
		);
	}
}

template <fwd_reference_to<vec> Vec>
constexpr particle_container::particle_container(Vec&& domain_arg, double cutoff_radius_arg)
		: domain_{std::forward<Vec>(domain_arg)}
		, grid_size_{div_round_up(domain_.x, cutoff_radius_arg), div_round_up(domain_.y, cutoff_radius_arg), div_round_up(domain_.z, cutoff_radius_arg)}
		, cutoff_radius_{cutoff_radius_arg} {
	// Prevent implicit widening later, if x * y * z overflows for uint32, this takes care of it
	const auto grid_x = static_cast<std::size_t>(grid_size_.x);
	const auto grid_y = static_cast<std::size_t>(grid_size_.y);
	const auto grid_z = static_cast<std::size_t>(grid_size_.z);
	TRACE_PARTICLE_CONTAINER("Setting up grid with size: {}", grid_size_);
	grid.resize(grid_x * grid_y * grid_z);
}

template <fwd_reference_to<particle> ParticleT>
constexpr void particle_container::place(ParticleT&& particle) {
	TRACE_PARTICLE_CONTAINER("Placing particle at coordinates: {}", particle.x);
	check_if_out_of_domain_max(particle.x);
	grid[pos_to_linear_index(particle.x)].emplace_back(std::forward<ParticleT>(particle));
}

template <fwd_reference_to<vec> Vec, typename... Args>
	requires std::constructible_from<particle, vec, Args...>
constexpr void particle_container::emplace(Vec&& position, Args&&... args) {
	TRACE_PARTICLE_CONTAINER("Emplacing particle at coordinates: {}", position);
	check_if_out_of_domain_max(position);
	grid[pos_to_linear_index(position)].emplace_back(
		std::forward<Vec>(position), std::forward<Args>(args)...
	);
}

template <std::size_t N>
constexpr void particle_container::add_cuboid(
	const vec& origin, const index& scale, double meshwidth, const vec& velocity, double mass,
	double brownian_mean, std::size_t& seq_no
) {
	for (size_type i = 0; i < scale.x; ++i) {
		for (size_type j = 0; j < scale.y; ++j) {
			for (size_type k = 0; k < scale.z; ++k) {
				emplace(
					vec{origin.x + (i * meshwidth), origin.y + (j * meshwidth),
				        origin.z + (k * meshwidth)},
					maxwell_boltzmann_distributed_velocity<N>(brownian_mean, seq_no) + velocity, mass,
					LOG_PARTICLE_TYPE ? seq_no : 0
				);
			}
		}
	}
}

constexpr void particle_container::add_disc(
	const vec& center, double radius, double meshwidth, const vec& velocity, double mass,
	double brownian_mean, std::size_t& seq_no
) {
	for (int i = static_cast<int>(-radius); i <= radius; i++) {
		for (int j = static_cast<int>(-radius); j <= radius; j++) {
			// TODO(tuna) check that this still works after the following simplification
			// if (std::sqrt(std::pow((i * meshwidth), 2) + std::pow((j * meshwidth), 2)) >
			// (radius * meshwidth)) {
			if (i * i + j * j > radius * radius) {
				continue;
			}
			emplace(
				vec{center.x + (i * meshwidth), center.y + (j * meshwidth), center.z},
				maxwell_boltzmann_distributed_velocity<2>(brownian_mean, seq_no) + velocity, mass,
				LOG_PARTICLE_TYPE ? seq_no : 0
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

constexpr range_of<particle> auto particle_container::particles() noexcept {
	return grid | std::views::join;
}

constexpr range_of<const particle> auto particle_container::particles() const noexcept {
	// TODO(tuna): i'm not sure as const is strictly necessary
	return std::as_const(grid) | std::views::join;
}

constexpr range_of<particle_container::cell> auto& particle_container::cells() noexcept {
	return grid;
}

constexpr const range_of<const particle_container::cell> auto& particle_container::cells() const noexcept {
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

constexpr particle_container::cell&
particle_container::cell_containing(const vec& position) noexcept {
	return grid[pos_to_linear_index(position)];
}

constexpr const particle_container::cell&
particle_container::cell_containing(const vec& pos) const noexcept {
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
