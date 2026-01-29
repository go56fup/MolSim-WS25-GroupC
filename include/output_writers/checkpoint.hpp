#pragma once

#include <fstream>
#include <vector>

#include <daw/json/daw_json_link.h>

#include "grid/particle_container/fwd.hpp"
#include "grid/particle_container/system.hpp"
#include "simulation/config/json_schema.hpp"
#include "simulation/entities.hpp"

constexpr std::string dump_state(particle_container& container) {
	const auto& system = container.system();

	// daw::json::json_value is non-owning, so the underlying string has to persist until the final
	// serialization.
	std::vector<std::string> geometries;
	geometries.resize(system.size());

	std::vector<body_entry> out;
	out.resize(system.size());

	for (particle_id i = 0; i < system.size(); ++i) {
		const particle_state_parameters state{
			.position = system.serialize_position(i),
			.force = system.serialize_force(i),
			.old_force = system.serialize_old_force(i)
		};
		TRACE_CHECKPOINT(
			"Serializing particle state: position={}, force={}, old_force={}", state.position,
			state.force, state.old_force
		);
		geometries[i] = daw::json::to_json(state);
		TRACE_CHECKPOINT("Serialized particle state at {}: {}", i, geometries[i]);
		const body_entry particle_entry{
			.type = "particle_state",
			// TODO(tuna): check if this is a view over a string that gets destructed at the end of
		    // the loop or is owning
			.geometry = daw::json::json_value(geometries[i]),
			.velocity = system.serialize_velocity(i),
			.material =
				{.mass = system.mass[i], .sigma = system.sigma[i], .epsilon = system.epsilon[i]}
		};
		out[i] = particle_entry;
	}
	return daw::json::
		to_json_array(out, daw::json::options::output_flags<daw::json::options::SerializationFormat::Pretty>);
};

inline void write_state_to_file(std::string_view state, std::string_view output_path) {
	std::ofstream out(output_path.data());
	out << state;
}
