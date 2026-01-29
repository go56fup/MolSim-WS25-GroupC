#pragma once

#include <stdexcept>
#include <string_view>

#include <daw/json/daw_json_link.h>
#include <spdlog/spdlog.h>

#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
#include "simulation/config/json_schema.hpp"
#include "simulation/entities.hpp"
#include "utility/macros.hpp"

namespace detail {
constexpr void check_for_2d_domain(const sim_configuration& config) noexcept(false) {
	if (config.dimensions == 3) {
		throw std::domain_error(
			"Refusing to create a 3D grid with a 2D object inside. Set "
			"the Z-size of the domain to the cutoff radius to obtain a 2D grid, "
			"or replace the 2D body with a 3D equivalent."
		);
	}
}

}  // namespace detail

namespace config {
constexpr sim_configuration parse_config(std::string_view json_data) noexcept(false) {
	TRACE_INPUT_PARSING("Using config: {}", json_data);
	const auto cfg = daw::json::from_json<
		sim_configuration>(json_data, daw::json::options::parse_flags<daw::json::options::UseExactMappingsByDefault::yes>);
	auto is_periodic = [&](boundary_type border) {
		return cfg.boundary_behavior[border] == boundary_condition::periodic;
	};
	for (auto [min, max] : border_pairs) {
		// TODO(tuna): add test that checks for this being covered
		if (is_periodic(min) != is_periodic(max)) {
			throw std::invalid_argument(fmt::format(
				"Mismatched periodic boundary condition between borders {} and {}.", min, max
			));
		}
	}
	return cfg;
}

constexpr std::vector<body_entry> parse_bodies(std::string_view json_data) {
	TRACE_INPUT_PARSING("Using bodies: {}", json_data);
	return daw::json::from_json_array<body_entry>(json_data);
}

constexpr void populate_simulation(
	particle_container& particles, const sim_configuration& config, std::span<body_entry> bodies
) noexcept(false) {
	std::once_flag two_d_domain_check_flag;
	std::size_t seq_no = 0;
	auto two_d_domain_check = [&] {
		std::call_once(two_d_domain_check_flag, detail::check_for_2d_domain, config);
	};
	auto decide_brownian = [&](body_entry& body) {
		STATIC_IF_NOT_TESTING const bool initial_temp_enforced =
			config.thermostat.has_value() && config.thermostat->enforce_initial_temperature;
		if (initial_temp_enforced) {
			assert(body.parameters->brownian_mean == 0.0);
			body.parameters->brownian_mean =
				std::sqrt(config.thermostat->initial_temperature / body.material.mass);
		}
	};
	for (auto& body : bodies) {
		TRACE_INPUT_PARSING("Parsing body with geometry: {}", daw::json::to_json(body.geometry));
		if (body.type == "particle") {
			particles.add_particle(
				body.geometry.as<particle_parameters>().position, body.velocity, body.material
			);
			continue;
		}
		if (body.type == "particle_state") {
			particles.reload_particle_state(
				body.geometry.as<particle_state_parameters>(), body.velocity, body.material
			);
			continue;
		}

		decide_brownian(body);

		if (body.type == "cuboid") {
			particles.add_cuboid<3>(
				body.geometry.as<cuboid_parameters<3>>(), body.parameters.value(), body.velocity,
				body.material, seq_no
			);
			continue;
		}
		if (body.type == "rectangle") {
			two_d_domain_check();
			particles.add_cuboid<2>(
				body.geometry.as<cuboid_parameters<2>>().extend_to_3d(config.domain),
				body.parameters.value(), body.velocity, body.material, seq_no
			);
			continue;
		}
		if (body.type == "disc") {
			two_d_domain_check();
			particles.add_disc(
				body.geometry.as<disc_parameters<2>>().extend_to_3d(config.domain),
				body.parameters.value(), body.velocity, body.material, seq_no
			);
			continue;
		}
		if (body.type == "membrane") {
			particles.add_membrane(
				body.geometry.as<cuboid_parameters<2>>().extend_to_3d(config.domain),
				body.parameters.value(), body.velocity, body.material, seq_no
			);
			continue;
		}
		throw std::invalid_argument(fmt::format("Unknown body type: {}", body.type));
	}
}

}  // namespace config
