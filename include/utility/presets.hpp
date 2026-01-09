#pragma once

#include <string_view>

#include "simulation/config/parse.hpp"
#include "simulation/molsim.hpp"

namespace presets {
namespace assignment4::task2 {
constexpr void big_experiment() {
#if HAS_EMBED
	constexpr char data[]{
#embed "assignment4/task2/big_experiment.json"
		, 0
	};
	constexpr auto json_data = std::string_view(data);
#else
	constexpr std::string_view json_data = R"({
  "configuration": {
    "delta_t": 0.0005,
    "cutoff_radius": 3,
    "boundary_conditions": {
      "x_min": "periodic",
      "y_min": "reflecting",
      "z_min": "reflecting",
      "x_max": "periodic",
      "y_max": "reflecting",
      "z_max": "reflecting"
    },
      "end_time": 0.5,
	  "thermostat": {
		  "initial_temperature": 40,
		  "application_frequency": 1000,
		  "enforce_initial_temperature": true
	  },
    "write_frequency": 10,
    "base_name": "MD_vtk",
      "domain": [
		  300,
		  54
          2.5
    ],
	  "gravitational_constant": -12.44,
    "create_checkpoint": false
  },
  "bodies": [
    {
      "type": "rectangle",
      "parameters": {
        "origin": [
          0.6,
          2
        ],
        "scale": [
          250,
          20
        ],
        "velocity": [
          0,
          0
        ],
        "brownian_mean": 0.1,
		"meshwidth": 1.2,
        "sigma": 1.2,
        "epsilon": 1,
        "particle_mass": 1
      }
    },
    {
      "type": "rectangle",
      "parameters": {
        "origin": [
			0.6,
			27
        ],
        "scale": [
			250,
            20
        ],
        "velocity": [
          0,
          0
        ],
        "brownian_mean": 0.1,
		"meshwidth": 1.2,
        "sigma": 1.1,
        "epsilon": 1,
        "particle_mass": 2
      }
    }
  ]
})";
#endif
	static constexpr auto parse_result = config::parse(json_data);
	particle_container container(parse_result.config.domain, parse_result.config.cutoff_radius);
	config::populate_simulation(container, parse_result.config, parse_result.bodies);
	run_simulation(container, parse_result.config, lennard_jones_force, "unused");
}
}  // namespace assignment4::task2
}  // namespace presets
