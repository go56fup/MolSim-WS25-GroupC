/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include <daw/json/daw_json_link.h>
#include <spdlog/spdlog.h>

#include "grid/particle_container/fwd.hpp"
#include "physics/vec_3d.hpp"
#include "simulation/config/entities.hpp"
#include "simulation/config/json_schema.hpp"
#include "utility/concepts.hpp"

namespace detail {
constexpr void check_for_2d_domain(const vec& domain, double cutoff) noexcept(false) {
	if (domain.z > cutoff) {
		throw std::domain_error(
			"Refusing to create a 3D grid with a 2D object inside. Set "
			"the Z-size of the domain to the cutoff radius to obtain a 2D grid, "
			"or replace the 2D body with a 3D equivalent."
		);
	}
}

}  // namespace detail

namespace config {
constexpr unprocessed_config parse(std::string_view json_data) {
	TRACE_INPUT_PARSING("Running from input content: {}", json_data);
	return daw::json::from_json<unprocessed_config>(
		json_data,
		daw::json::options::parse_flags<daw::json::options::UseExactMappingsByDefault::yes>
	);
}

constexpr void populate_simulation(
	particle_container& particles, const unprocessed_config& desc
) noexcept(false) {
	std::once_flag two_d_domain_check;
	std::size_t seq_no = 0;
	for (const auto& [_, body] : desc.bodies) {
		const std::string type = body["type"].as<std::string>();
		if (type == "cuboid") {
			const auto params = body["parameters"].as<cuboid_parameters<3>>();
			particles.add_cuboid<3>(
				params.origin, params.scale, desc.config.meshwidth, params.velocity,
				params.particle_mass, params.brownian_mean, seq_no
			);
		} else if (type == "square") {
			std::call_once(
				two_d_domain_check, detail::check_for_2d_domain, desc.config.domain,
				desc.config.cutoff_radius
			);
			const auto params = body["parameters"].as<cuboid_parameters<2>>();
			particles.add_cuboid<2>(
				{params.origin.x, params.origin.y, desc.config.domain.z / 2},
				{params.scale.x, params.scale.y, 1}, desc.config.meshwidth, params.velocity,
				params.particle_mass, params.brownian_mean, seq_no
			);
		} else if (type == "particle") {
			const auto params = body["parameters"].as<particle>();
			particles.place(MOVE_IF_DEBUG(params));
		} else if (type == "disc") {
			std::call_once(
				two_d_domain_check, detail::check_for_2d_domain, desc.config.domain,
				desc.config.cutoff_radius
			);
			const auto params = body["parameters"].as<disc_parameters>();
			particles.add_disc(
				{params.center.x, params.center.y, desc.config.domain.z / 2}, params.radius,
				desc.config.meshwidth, params.velocity, params.particle_mass, params.brownian_mean,
				seq_no
			);
		} else {
			throw std::invalid_argument(fmt::format("Unknown body type: {}", type));
		}
	}
}
}  // namespace config
