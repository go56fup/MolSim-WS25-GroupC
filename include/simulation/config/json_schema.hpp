#pragma once

#include <daw/json/daw_json_data_contract.h>
#include <daw/json/daw_json_link.h>
#include <daw/json/daw_json_link_types.h>

#include "grid/particle_container/fwd.hpp"
#include "physics/vec_3d.hpp"
#include "simulation/entities.hpp"

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

template <std::size_t Dimensions>
struct json_data_contract<cuboid_parameters<Dimensions>> {
	using type = json_member_list<
		json_class<"origin", vec_n<Dimensions, double>>,
		json_class<"scale", vec_n<Dimensions, particle_container::size_type>>>;
};

template <>
struct json_data_contract<particle_parameters> {
	using type = json_member_list<json_class<"position", vec>>;
};

template <>
struct json_data_contract<particle_state_parameters> {
	using type = json_member_list<
		json_class<"position", vec>, json_class<"force", vec>, json_class<"old_force", vec>>;

	static constexpr auto to_json_data(const particle_state_parameters& p) noexcept {
		return std::forward_as_tuple(p.position, p.force, p.old_force);
	}
};

template <>
struct json_data_contract<disc_parameters<2>> {
	using type =
		json_member_list<json_class<"center", vec_2d<double>>, json_number<"radius", double>>;
};

template <>
struct json_data_contract<thermostat_parameters> {
	using constructor_t = thermostat_parameters_constructor;
	using type = json_member_list<
		json_number<"initial_temperature", double>,
		json_number<"application_frequency", sim_iteration_t>,
		// TODO(tuna): can you also just say double?
		json_number_null<"target_temperature", std::optional<double>>,
		json_number_null<"max_temperature_difference", std::optional<double>>,
		json_bool<"enforce_initial_temperature">>;
};

template <>
struct json_data_contract<body_common_parameters> {
	using constructor_t = body_common_parameters_constructor;
	using type = json_member_list<
		json_number<"meshwidth", double>, json_number_null<"brownian_mean", std::optional<double>>>;

	static constexpr auto to_json_data(const body_common_parameters& params) noexcept {
		return std::forward_as_tuple(params.meshwidth, params.brownian_mean);
	};
};

template <>
struct json_data_contract<material_description> {
	using type = json_member_list<
		json_number<"mass", double>, json_number<"sigma", double>, json_number<"epsilon", double>>;

	static constexpr auto to_json_data(const material_description& material) noexcept {
		return std::forward_as_tuple(material.mass, material.sigma, material.epsilon);
	};
};

template <>
struct json_data_contract<body_entry> {
	using type = json_member_list<
		json_string<"type">, json_raw<"geometry">,
		json_class_null<"parameters", std::optional<body_common_parameters>>,
		json_class<"velocity", vec>, json_class<"material", material_description>>;

	static constexpr auto to_json_data(const body_entry& body) noexcept {
		return std::forward_as_tuple(
			body.type, body.geometry, body.parameters, body.velocity, body.material
		);
	}
};
}  // namespace daw::json
