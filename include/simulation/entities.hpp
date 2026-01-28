#pragma once

#include <limits>
#include <optional>

#include <daw/json/daw_json_link.h>
#include <utility>

#include "grid/enums.hpp"
#include "physics/vec_3d.hpp"
#include "utility/concepts.hpp"

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
		}
		assert(false && "This type is not a boundary and therefore has no associated condition");
		std::unreachable();
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

class boundary_conditions_descriptor_constructor {
private:
	static constexpr boundary_condition get_enum(std::string_view condition_desc) {
		using enum boundary_condition;
		if (condition_desc == "outflow") return outflow;
		if (condition_desc == "reflecting") return reflecting;
		if (condition_desc == "periodic") return periodic;
		throw std::invalid_argument(fmt::format("Unknown boundary condition: {}", condition_desc));
	}

public:
	static constexpr boundary_conditions_descriptor operator()(
		std::string_view x_min, std::string_view y_min, std::string_view z_min,
		std::string_view x_max, std::string_view y_max, std::string_view z_max
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

using sim_iteration_t = unsigned;

struct thermostat_parameters {
	double initial_temperature;
	sim_iteration_t application_frequency;
	double target_temperature;
	double max_temperature_difference = std::numeric_limits<double>::infinity();
	bool enforce_initial_temperature;
};

struct thermostat_parameters_constructor {
	static constexpr thermostat_parameters operator()(
		double initial_temperature, sim_iteration_t application_frequency,
		std::optional<double> target_temperature, std::optional<double> max_temperature_difference,
		bool enforce_initial_temperature
	) {
		return thermostat_parameters{
			.initial_temperature = initial_temperature,
			.application_frequency = application_frequency,
			.target_temperature = target_temperature.value_or(initial_temperature),
			.max_temperature_difference =
				max_temperature_difference.value_or(std::numeric_limits<double>::infinity()),
			.enforce_initial_temperature = enforce_initial_temperature
		};
	}
};

struct material_description {
	double mass;
	double sigma;
	double epsilon;

	constexpr bool operator==(const material_description&) const noexcept = default;
};

struct particle_properties {
	vec velocity;
	material_description material{};
};

struct sim_configuration {
	static constexpr auto max_base_name_len = 32;
	double delta_t{};
	double cutoff_radius{};
	boundary_conditions_descriptor boundary_behavior{};
	std::optional<thermostat_parameters> thermostat = std::nullopt;
	double end_time{};
	sim_iteration_t write_frequency = std::numeric_limits<sim_iteration_t>::max();
	p3094::fixed_string<max_base_name_len> base_name{std::from_range, "unused"};
	vec domain;
	bool create_checkpoint = false;
	std::uint8_t dimensions = 3;
	double gravitational_constant = 0;
};

struct sim_configuration_constructor {
	template <
		typename String, fwd_reference_to<vec> Vec,
		fwd_reference_to<boundary_conditions_descriptor> Descriptor,
		fwd_reference_to<std::optional<thermostat_parameters>> Thermostat>
	static constexpr sim_configuration operator()(
		double delta_t, double cutoff_radius, Descriptor&& conditions, Thermostat&& thermostat,
		double end_time, sim_iteration_t write_frequency, String&& base_name, Vec&& domain,
		bool create_checkpoint, std::optional<double> gravity
	) {
		const decltype(sim_configuration::dimensions) dims = domain.z > cutoff_radius ? 3 : 2;
		return sim_configuration{
			delta_t,
			cutoff_radius,
			std::forward<Descriptor>(conditions),
			std::forward<Thermostat>(thermostat),
			end_time,
			write_frequency,
			{std::from_range, std::forward<String>(base_name)},
			std::forward<Vec>(domain),
			create_checkpoint,
			dims,
			gravity.value_or(0)
		};
	}
};

template <typename T>
struct vec_2d {
	T x;
	T y;

	constexpr bool operator==(const vec_2d<T>&) const noexcept = default;
};

template <std::size_t N, typename T>
	requires(N == 2 || N == 3)
using vec_n = std::conditional_t<N == 2, vec_2d<T>, vec_3d<T>>;

struct body_common_parameters {
	double meshwidth;
	double brownian_mean = 0;
};

struct body_common_parameters_constructor {
	static constexpr body_common_parameters
	operator()(double meshwidth, std::optional<double> brownian_mean) {
		return {.meshwidth = meshwidth, .brownian_mean = brownian_mean.value_or(0)};
	}
};

template <std::size_t Dimensions>
struct cuboid_parameters {
	vec_n<Dimensions, double> origin;
	// TODO(tuna): link this back to particle_container without introducing circular deps
	vec_n<Dimensions, std::uint32_t> scale;

	constexpr cuboid_parameters<3> extend_to_3d(const vec& domain) const {
		return {.origin = {origin.x, origin.y, domain.z / 2}, .scale = {scale.x, scale.y, 1}};
	}
};

template <std::size_t N>
struct disc_parameters {
	vec_n<N, double> center;
	double radius;

	constexpr disc_parameters<3> extend_to_3d(const vec& domain) const {
		return {.center = {center.x, center.y, domain.z / 2}, .radius = radius};
	}
};

inline constexpr std::array<std::string_view, 1> three_d_objects{"cuboid"};
inline constexpr std::array<std::string_view, 2> two_d_objects{"disc", "rectangle"};
inline constexpr std::array<std::string_view, 2> one_d_objects{"particle", "particle_state"};

struct body_entry {
	std::string type;
	daw::json::json_value geometry;
	std::optional<body_common_parameters> parameters = std::nullopt;
	vec velocity;
	material_description material{};
};

struct particle_parameters {
	vec position;
};

struct particle_state_parameters {
	vec position;
	vec force;
	vec old_force;
};
