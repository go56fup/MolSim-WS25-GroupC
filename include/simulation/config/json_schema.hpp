#pragma once

#include <daw/json/daw_json_link.h>
#include <daw/json/daw_json_link_types.h>

#include "entities.hpp"
#include "simulation/config/entities.hpp"

namespace daw::json {
template <typename T>
struct json_data_contract<vec_3d<T>> {
	using type = json_tuple_member_list<T, T, T>;

	static constexpr auto to_json_data(const vec_3d<T>& v) noexcept {
		return std::forward_as_tuple(v.x, v.y, v.z);
	}
};

template <typename T>
struct json_data_contract<vec_2d<T>> {
	using type = json_tuple_member_list<T, T>;
};

template <>
struct json_data_contract<boundary_conditions_descriptor> {
	using constructor_t = boundary_conditions_descriptor_constructor;
	using type = json_member_list<
		json_string<"x_min">, json_string<"y_min">, json_string<"z_min">, json_string<"x_max">,
		json_string<"y_max">, json_string<"z_max">>;
};

template <>
struct json_data_contract<sim_configuration> {
	using constructor_t = sim_configuration_constructor;
	using type = json_member_list<
		json_number<"delta_t", double>, json_number<"cutoff_radius", double>,
		json_class<"boundary_conditions", boundary_conditions_descriptor>,
		json_class_null<"thermostat", std::optional<thermostat_parameters>>,
		json_number<"end_time", double>, json_number<"write_frequency", sim_iteration_t>,
		json_string<"base_name">, json_class<"domain", vec>, json_bool<"create_checkpoint">,
		json_number_null<"gravitational_constant", std::optional<double>>>;
};

template <>
struct json_data_contract<unprocessed_config> {
	using type =
		json_member_list<json_class<"configuration", sim_configuration>, json_raw<"bodies">>;
};

template <std::size_t Dimensions>
struct json_data_contract<cuboid_parameters<Dimensions>> {
	using type = json_member_list<
		json_class<"origin", vec_n<Dimensions, double>>,
		json_class<"scale", vec_n<Dimensions, particle_container::size_type>>,
		json_class<"velocity", vec_n<Dimensions, double>>, json_number<"meshwidth", double>,
		json_number<"brownian_mean", double>, json_number<"particle_mass", double>,
		json_number<"sigma", double>, json_number<"epsilon", double>>;
};

template <>
struct json_data_contract<particle_parameters> {
	using type = json_member_list<
		json_class<"position", vec>, json_class<"velocity", vec>, json_number<"mass", double>,
		json_number<"sigma", double>, json_number<"epsilon", double>>;
};

template <>
struct json_data_contract<particle> {
	using type = json_member_list<
		json_class<"position", vec>, json_class<"velocity", vec>, json_class<"force", vec>,
		json_class<"old_force", vec>, json_number<"mass", double>, json_number<"sigma", double>,
		json_number<"epsilon", double>, json_number<"type", decltype(particle::type)>>;

	static constexpr auto to_json_data(const particle& p) noexcept {
		return std::forward_as_tuple(p.x, p.v, p.f, p.old_f, p.m, p.sigma, p.epsilon, p.type);
	}
};

template <>
struct json_data_contract<disc_parameters> {
	using type = json_member_list<
		json_class<"center", vec_2d<double>>, json_number<"radius", double>,
		json_class<"velocity", vec_2d<double>>, json_number<"meshwidth", double>,
		json_number<"brownian_mean", double>, json_number<"particle_mass", double>,
		json_number<"sigma", double>, json_number<"epsilon", double>>;
};

template <typename Body>
struct json_data_contract<serialization_view<Body>> {
	using type = json_member_list<json_string<"type">, json_class<"parameters", Body>>;

	static constexpr auto to_json_data(const serialization_view<Body>& body) {
		return std::forward_as_tuple(body.type, body.parameters);
	}
};

template <>
struct json_data_contract<thermostat_parameters> {
	using constructor_t = thermostat_parameters_constructor;
	using type = json_member_list<
		json_number<"initial_temperature", double>,
		json_number<
			"application_frequency", sim_iteration_t>,
		json_number<"target_temperature", double>,
		json_number_null<"max_temperature_difference", std::optional<double>>
		>;
};
}  // namespace daw::json
