#pragma once

#include <fstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include <daw/json/daw_json_link.h>
#include <spdlog/spdlog.h>

#include "entities.hpp"
#include "grid/enums.hpp"
#include "grid/particle_container/particle_container.hpp"
#include "physics/vec_3d.hpp"
#include "simulation/config/entities.hpp"
#include "simulation/config/json_schema.hpp"
#include "utility/concepts.hpp"
#include "utility/constants.hpp"

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
constexpr unprocessed_config parse(std::string_view json_data) noexcept(false) {
	TRACE_INPUT_PARSING("Running from input content: {}", json_data);
	auto cfg = daw::json::from_json<
		unprocessed_config>(json_data, daw::json::options::parse_flags<daw::json::options::UseExactMappingsByDefault::yes>);
	auto is_periodic = [&](boundary_type border) {
		return cfg.config.boundary_behavior[border] == boundary_condition::periodic;
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

// TODO(tuna): split the configuration into global and body files, so the checkpointing basically
// becomes a body file generator
constexpr void populate_simulation(
	particle_container& particles, const sim_configuration& config,
	const daw::json::json_value& bodies
) noexcept(false) {
	std::once_flag two_d_domain_check;
	std::size_t seq_no = 0;
	auto decide_brownian = [&](const auto& params) {
		static const bool initial_temp_enforced = config.thermostat.has_value() && config.thermostat->enforce_initial_temperature;
		if (initial_temp_enforced) {
			return std::sqrt(config.thermostat->initial_temperature / params.particle_mass);
		}
		return params.brownian_mean;
	};

	for (const auto& [_, body] : bodies) {
		const std::string type = body["type"].as<std::string>();
		if (type == "cuboid") {
			const auto params = body["parameters"].as<cuboid_parameters<3>>();
			particles.add_cuboid<3>(
				params.origin, params.scale, params.meshwidth, params.velocity,
				params.particle_mass, params.sigma, params.epsilon, decide_brownian(params), seq_no
			);
		} else if (type == "rectangle") {
			std::call_once(two_d_domain_check, detail::check_for_2d_domain, config);
			const auto params = body["parameters"].as<cuboid_parameters<2>>();
			particles.add_cuboid<2>(
				{params.origin.x, params.origin.y, config.domain.z / 2},
				{params.scale.x, params.scale.y, 1}, params.meshwidth, params.velocity,
				params.particle_mass, params.sigma, params.epsilon, decide_brownian(params), seq_no
			);
		} else if (type == "particle") {
			const auto params = body["parameters"].as<particle_parameters>();
			particles.emplace(
				params.position, params.velocity, params.mass, params.sigma, params.epsilon
			);
		} else if (type == "disc") {
			std::call_once(two_d_domain_check, detail::check_for_2d_domain, config);
			const auto params = body["parameters"].as<disc_parameters>();
			particles.add_disc(
				{params.center.x, params.center.y, config.domain.z / 2}, params.radius,
				params.meshwidth, params.velocity, params.particle_mass, params.sigma,
				params.epsilon, decide_brownian(params), seq_no
			);
		} else if (type == "particle_state") {
			particles.place(body["parameters"].as<particle>());
		} else {
			throw std::invalid_argument(fmt::format("Unknown body type: {}", type));
		}
	}
}

// TODO(tuna): move to output part of the codebase
constexpr std::string dump_state(particle_container& container) {
	std::string out = "[\n";
	for (const particle& p : container.particles()) {
		const auto view = serialization_view<particle>{.type = "particle_state", .parameters = p};
		const auto json = daw::json::
			to_json(view, daw::json::options::output_flags<daw::json::options::SerializationFormat::Pretty>);
		fmt::format_to(std::back_inserter(out), "{},\n", json);
	}
	out.pop_back();
	out.pop_back();
	out += ']';
	return out;
};

inline void write_state_to_file(std::string_view state, std::string_view output_path) {
	std::ofstream out(output_path.data());
	out << state;
}

}  // namespace config
