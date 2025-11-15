/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <spdlog/spdlog.h>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include <daw/json/daw_json_link.h>

#include "Concepts.h"
#include "ParticleContainer.h"
#include "Vector.h"

enum class boundary_condition : std::uint8_t { outflow, reflecting };

struct boundary_conditions_descriptor {
	boundary_condition x_min;
	boundary_condition y_min;
	boundary_condition z_min;
	boundary_condition x_max;
	boundary_condition y_max;
	boundary_condition z_max;

	constexpr boundary_condition operator[](boundary_type type) const noexcept {
		switch (type) {
		case boundary_type::x_min:
			return x_min;
		case boundary_type::y_min:
			return y_min;
		case boundary_type::z_min:
			return z_min;
		case boundary_type::x_max:
			return x_max;
		case boundary_type::y_max:
			return y_max;
		case boundary_type::z_max:
			return z_max;
		case boundary_type::none:
			assert(false && "This type is not a boundary and therefore has no associated condition");
		}
		std::unreachable();
		// Having the compiler optimize a switch statement is much less cumbersome than doing it by hand.
		// But constexpr won't allow it.
		// TODO(tuna): see why the below does not work
		// const auto underlying = std::to_underlying(type);
		// assert(underlying <= std::to_underlying(boundary_type::z_max) && type != boundary_type::none);
		// NOLINTNEXTLINE(*pointer-arithmetic)
		// return (&x_min)[underlying - 1];
	}

	static constexpr boundary_conditions_descriptor all(boundary_condition condition) noexcept {
		return {
			.x_min = condition,
			.y_min = condition,
			.z_min = condition,
			.x_max = condition,
			.y_max = condition,
			.z_max = condition
		};
	}

	constexpr bool operator==(const boundary_conditions_descriptor&) const noexcept = default;
};

struct sim_configuration {
	double delta_t;
	double cutoff_radius;
	double sigma;
	double epsilon;
	boundary_conditions_descriptor boundary_behavior;
	double end_time;
	unsigned write_frequency;
	std::string base_name;
	vec domain;
	double meshwidth;
};

class boundary_conditions_descriptor_constructor {
private:
	static constexpr boundary_condition get_enum(std::string_view condition_desc) {
		using enum boundary_condition;
		if (condition_desc == "outflow") return outflow;
		if (condition_desc == "reflecting") return reflecting;
		throw std::invalid_argument(fmt::format("Unknown boundary condition: {}", condition_desc));
	}

public:
	static constexpr boundary_conditions_descriptor operator()(
		std::string_view x_min, std::string_view y_min, std::string_view z_min, std::string_view x_max,
		std::string_view y_max, std::string_view z_max
	) {
		return {
			.x_min = get_enum(x_min),
			.y_min = get_enum(y_min),
			.z_min = get_enum(z_min),
			.x_max = get_enum(x_max),
			.y_max = get_enum(y_max),
			.z_max = get_enum(z_max)
		};
	};
};

class sim_configuration_constructor {
private:
	static constexpr double sixth_root_of_2 = 1.12246204830937298143353304967917951623241111061398;

public:
	template <
		fwd_reference_to<std::string> String, fwd_reference_to<vec> Vec,
		fwd_reference_to<boundary_conditions_descriptor> Descriptor>
	static constexpr sim_configuration operator()(
		double delta_t, double cutoff_radius, double sigma, double epsilon, Descriptor&& conditions, double end_time,
		unsigned write_frequency, String&& base_name, Vec&& domain
	) {
		return sim_configuration{
			delta_t,
			cutoff_radius,
			sigma,
			epsilon,
			std::forward<Descriptor>(conditions),
			end_time,
			write_frequency,
			std::forward<String>(base_name),
			std::forward<Vec>(domain),
			sixth_root_of_2 * sigma
		};
	}
};

template <typename T>
struct vec_2d {
	T x;
	T y;

	// NOLINTNEXTLINE(*explicit*)
	constexpr explicit(false) operator vec_3d<T>() const noexcept {
		return {x, y, 0};
	}
};

template <std::size_t N, typename T>
	requires(N == 2 || N == 3)
using vec_n = std::conditional_t<N == 2, vec_2d<T>, vec_3d<T>>;

template <std::size_t Dimensions>
struct cuboid_parameters {
	vec_n<Dimensions, double> origin;
	vec_n<Dimensions, ParticleContainer::size_type> scale;
	vec_n<Dimensions, double> velocity;
	double brownian_mean = -1;
	double particle_mass = -1;
};

struct disc_parameters {
	vec_2d<double> center;
	double radius;
	vec_2d<double> velocity;
	double particle_mass;
};

namespace daw::json {
template <typename T>
struct json_data_contract<vec_3d<T>> {
	using type = json_tuple_member_list<T, T, T>;
};

template <typename T>
struct json_data_contract<vec_2d<T>> {
	using type = json_tuple_member_list<T, T>;
};

template <>
struct json_data_contract<boundary_conditions_descriptor> {
	using constructor_t = boundary_conditions_descriptor_constructor;
	using type = json_member_list<
		json_string<"x_min">, json_string<"y_min">, json_string<"z_min">, json_string<"x_max">, json_string<"y_max">,
		json_string<"z_max">>;
};

template <>
struct json_data_contract<sim_configuration> {
	using constructor_t = sim_configuration_constructor;
	using type = json_member_list<
		json_number<"delta_t", double>, json_number<"cutoff_radius", double>, json_number<"sigma", double>,
		json_number<"epsilon", double>, json_class<"boundary_conditions", boundary_conditions_descriptor>,
		json_number<"end_time", double>, json_number<"write_frequency", unsigned>, json_string<"base_name">,
		json_class<"domain", vec>>;
};

template <std::size_t Dimensions>
struct json_data_contract<cuboid_parameters<Dimensions>> {
	using type = json_member_list<
		json_class<"origin", vec_n<Dimensions, double>>,
		json_class<"scale", vec_n<Dimensions, ParticleContainer::size_type>>,
		json_class<"velocity", vec_n<Dimensions, double>>, json_number<"brownian_mean", double>,
		json_number<"particle_mass", double>>;
};

template <>
struct json_data_contract<Particle> {
	using type = json_member_list<
		json_class<"position", vec>, json_class<"velocity", vec>, json_number<"mass", double>,
		json_number<"type", decltype(Particle::type)>>;
};

template <>
struct json_data_contract<disc_parameters> {
	using type = json_member_list<
		json_class<"center", vec_2d<double>>, json_number<"radius", double>, json_class<"velocity", vec_2d<double>>,
		json_number<"particle_mass", double>>;
};

}  // namespace daw::json

namespace FileReader {
constexpr std::pair<sim_configuration, ParticleContainer> parse(std::string_view json_data) noexcept(false) {
	SPDLOG_TRACE("Running from input content: {}", json_data);
	const auto json_val = daw::json::from_json<daw::json::json_value>(
		json_data, daw::json::options::parse_flags<daw::json::options::UseExactMappingsByDefault::yes>
	);
	auto sim_config = json_val["configuration"].as<sim_configuration>();
	const auto& domain = sim_config.domain;
	ParticleContainer particles(domain.x, domain.y, domain.z, sim_config.cutoff_radius);
	std::size_t seq_no = 0;
	for (const auto& [_, body] : json_val["bodies"]) {
		const std::string type = body["type"].as<std::string>();
		// TODO(tuna): this might not actually be winning us anything because
		// strcmp is simd so the whole check is proabably pretty cheap
		// the construction of the string_view and the subtraction to see where strcmp oughtta start
		// is not negligable work
		if (type.starts_with("cuboid_")) {
			if (type.ends_with("3d")) {
				const auto params = body["parameters"].as<cuboid_parameters<3>>();
				particles.add_cuboid<3>(
					params.origin, params.scale, sim_config.meshwidth, params.velocity, params.particle_mass,
					params.brownian_mean, seq_no
				);
				// TODO(tuna): rename to square
			} else if (type.ends_with("2d")) {
				const auto params = body["parameters"].as<cuboid_parameters<2>>();
				// TODO(tuna): figure out where on the z axis to place 2d objects.
				particles.add_cuboid<2>(
					{params.origin.x, params.origin.y, sim_config.domain.z / 2}, {params.scale.x, params.scale.y, 1},
					sim_config.meshwidth, params.velocity, params.particle_mass, params.brownian_mean, seq_no
				);
			}
		} else if (type == "particle") {
			const auto params = body["parameters"].as<Particle>();
			particles.place(MOVE_IF_DEBUG(params));
		} else if (type == "disc") {
			const auto params = body["parameters"].as<disc_parameters>();
			particles.add_disc(
				{params.center.x, params.center.y, sim_config.domain.z / 2}, params.radius, sim_config.meshwidth,
				params.velocity, params.particle_mass
			);
		} else {
			throw std::invalid_argument(fmt::format("Unknown body type: {}", type));
		}
	}
	return {std::move(sim_config), std::move(particles)};
}
}  // namespace FileReader
