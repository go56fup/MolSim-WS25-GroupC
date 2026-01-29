#pragma once

#include <utility>
#include <vector>

#include "physics/vec_3d.hpp"
#include "simulation/entities.hpp"
#include "simulation/particle.hpp"

#define system_arithmetic(op, idx, value, x_comp, y_comp, z_comp)                                  \
	do {                                                                                           \
		(x_comp)[idx] op value.x;                                                                  \
		(y_comp)[idx] op value.y;                                                                  \
		(z_comp)[idx] op value.z;                                                                  \
	} while (0)

class particle_system {
private:
	constexpr void destructure_push_back(
		const vec& value, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z
	) {
		x.push_back(value.x);
		y.push_back(value.y);
		z.push_back(value.z);
	}

	constexpr void
	add_common(const vec& pos, const vec& velocity, const material_description& material) {
		destructure_push_back(pos, x, y, z);
		destructure_push_back(velocity, vx, vy, vz);
		mass.push_back(material.mass);
		sigma.push_back(material.sigma);
		epsilon.push_back(material.epsilon);
	}

public:
	std::vector<double> x, y, z;
	std::vector<double> vx, vy, vz;
	std::vector<double> fx, fy, fz;
	std::vector<double> old_fx, old_fy, old_fz;
	std::vector<double> mass;
	std::vector<double> epsilon;
	std::vector<double> sigma;

	enum class property : std::uint8_t { position, velocity, force, old_force };

	constexpr void
	add_particle(const vec& pos, const vec& velocity, const material_description& material) {
		add_common(pos, velocity, material);
		// TODO(tuna): add finalization pass for numa
		fx.push_back(0);
		fy.push_back(0);
		fz.push_back(0);

		old_fx.push_back(0);
		old_fy.push_back(0);
		old_fz.push_back(0);
	}

	constexpr void add_particle(
		const vec& pos, const vec& velocity, const vec& force, const vec& old_force,
		const material_description& material
	) {
		add_common(pos, velocity, material);
		destructure_push_back(force, fx, fy, fz);
		destructure_push_back(old_force, old_fx, old_fy, old_fz);
	}

	constexpr vec serialize_force(particle_id i) const noexcept {
		return {fx[i], fy[i], fz[i]};
	}

	constexpr vec serialize_velocity(particle_id i) const noexcept {
		return {vx[i], vy[i], vz[i]};
	}

	constexpr vec serialize_old_force(particle_id i) const noexcept {
		return {old_fx[i], old_fy[i], old_fz[i]};
	}

	constexpr vec serialize_position(particle_id i) const noexcept {
		return {x[i], y[i], z[i]};
	}

	constexpr std::span<const double> position_component(axis a) const noexcept {
		switch (a) {
		case axis::x:
			return x;
		case axis::y:
			return y;
		case axis::z:
			return z;
		}
		std::unreachable();
	}

	constexpr std::span<double> position_component(axis a) noexcept {
		switch (a) {
		case axis::x:
			return x;
		case axis::y:
			return y;
		case axis::z:
			return z;
		}
		std::unreachable();
	}

	constexpr std::span<double> force_component(axis a) noexcept {
		switch (a) {
		case axis::x:
			return fx;
		case axis::y:
			return fy;
		case axis::z:
			return fz;
		}
		std::unreachable();
	}

	constexpr std::size_t size() const noexcept {
		assert(
			x.size() == y.size() && y.size() == z.size() && z.size() == fx.size() &&
			fx.size() == fy.size() && fy.size() == fz.size() && fz.size() == old_fx.size() &&
			old_fx.size() == old_fy.size() && old_fy.size() == old_fz.size() &&
			old_fz.size() == mass.size() && mass.size() == sigma.size() &&
			sigma.size() == epsilon.size()
		);
		return x.size();
	}
};
