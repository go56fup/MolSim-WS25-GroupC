#pragma once

#include <daw/json/daw_json_link.h>

#include "grid/enums.hpp"
#include "grid/particle_container/fwd.hpp"
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

struct sim_configuration {
	double delta_t;
	double cutoff_radius;
	boundary_conditions_descriptor boundary_behavior;
	double end_time;
	unsigned write_frequency;
	p3094::fixed_string<32> base_name;
	vec domain;
	bool create_checkpoint;
};

struct sim_configuration_constructor {
	template <
		typename String, fwd_reference_to<vec> Vec,
		fwd_reference_to<boundary_conditions_descriptor> Descriptor>
	static constexpr sim_configuration operator()(
		double delta_t, double cutoff_radius, Descriptor&& conditions, double end_time,
		unsigned write_frequency, String&& base_name, Vec&& domain, bool create_checkpoint
	) {
		return sim_configuration{
			delta_t,
			cutoff_radius,
			std::forward<Descriptor>(conditions),
			end_time,
			write_frequency,
			{std::from_range, std::forward<String>(base_name)},
			std::forward<Vec>(domain),
			create_checkpoint
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

// TODO(tuna): factor out meshwidth-brownian etc. into subclass that is shared
template <std::size_t Dimensions>
struct cuboid_parameters {
	vec_n<Dimensions, double> origin;
	vec_n<Dimensions, particle_container::size_type> scale;
	vec_n<Dimensions, double> velocity;
	double meshwidth;
	double brownian_mean;
	double particle_mass;
	double sigma;
	double epsilon;
};

struct disc_parameters {
	vec_2d<double> center;
	double radius;
	vec_2d<double> velocity;
	double meshwidth;
	double brownian_mean;
	double particle_mass;
	double sigma;
	double epsilon;
};

struct unprocessed_config {
	sim_configuration config;
	daw::json::json_value bodies;
};

struct particle_parameters {
	vec position;
	vec velocity;
	double mass;
	double sigma;
	double epsilon;
};

template <typename Body>
struct serialization_view {
	std::string_view type;
	const Body& parameters;
};
