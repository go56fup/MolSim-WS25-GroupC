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

#include "daw/json/daw_json_link.h"
#include "ParticleContainer.h"
#include "Vector.h"

struct sim_configuration {
	double delta_t;
	double cutoff_radius;
	double end_time;
	unsigned write_frequency;
	std::string base_name;
	vec domain;
};

struct cuboid_particle_parameters {
	double meshwidth;
	double mass;
	constexpr bool operator==(const cuboid_particle_parameters&) const = default;
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
	cuboid_particle_parameters particle_params{.meshwidth = -1, .mass = -1};
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
struct json_data_contract<sim_configuration> {
	using type = json_member_list<
		json_number<"delta_t", double>, json_number<"cutoff_radius", double>, json_number<"end_time", double>,
		json_number<"write_frequency", unsigned>, json_string<"base_name">, json_class<"domain", vec>>;
};

template <>
struct json_data_contract<cuboid_particle_parameters> {
	using type = json_member_list<json_number<"meshwidth", double>, json_number<"mass", double>>;
};

template <std::size_t Dimensions>
struct json_data_contract<cuboid_parameters<Dimensions>> {
	using type = json_member_list<
		json_class<"origin", vec_n<Dimensions, double>>,
		json_class<"scale", vec_n<Dimensions, ParticleContainer::size_type>>,
		json_class<"velocity", vec_n<Dimensions, double>>, json_number<"brownian_mean", double>,
		json_class<"particles", cuboid_particle_parameters>>;
};

template <>
struct json_data_contract<Particle> {
	using type = json_member_list<
		json_class<"position", vec>, json_class<"velocity", vec>, json_number<"mass", double>,
		json_number<"type", decltype(Particle::type)>>;
};
}  // namespace daw::json

namespace FileReader {
constexpr std::pair<sim_configuration, ParticleContainer> parse(std::string_view json_data) noexcept(false) {
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
					params.origin, params.scale, params.particle_params.meshwidth, params.velocity,
					params.particle_params.mass, params.brownian_mean, seq_no
				);
			} else if (type.ends_with("2d")) {
				const auto params = body["parameters"].as<cuboid_parameters<2>>();
				particles.add_cuboid<2>(
					params.origin, params.scale, params.particle_params.meshwidth, params.velocity,
					params.particle_params.mass, params.brownian_mean, seq_no
				);
			}
		} else if (type == "particle") {
			const auto params = body["parameters"].as<Particle>();
			particles.place(MOVE_IF_DEBUG(params));
		} else {
			throw std::invalid_argument(fmt::format("Unknown body type: {}", type));
		}
	}
	return {std::move(sim_config), std::move(particles)};
}
}  // namespace FileReader
