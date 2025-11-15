/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "MolSim.h"
#include "ParticleContainer.h"
#include "Vector.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace detail {
// TODO(tuna): add proper error handling during parse, maybe even just switch over to scnlib

void parse_particle(ParticleContainer& particles, std::istringstream& datastream) {
	double x_x{};
	double x_y{};
	double x_z{};
	double v_x{};
	double v_y{};
	double v_z{};

	double m{};

	datastream >> x_x >> x_y >> x_z >> v_x >> v_y >> v_z >> m;
	particles.emplace_back(
		std::piecewise_construct, std::forward_as_tuple(x_x, x_y, x_z), std::forward_as_tuple(v_x, v_y, v_z), m
	);
}

void parse_cuboid(ParticleContainer& particles, std::istringstream& datastream) {
	vec origin{};
	vec_3d<int> scale{};
	double distance{};
	vec velocity{};
	double mass{};
	double brownian_mean{};

	for (auto& i : origin) {
		datastream >> i;
	}

	for (auto& i : scale) {
		datastream >> i;
	}

	datastream >> distance;
	for (auto& i : velocity) {
		datastream >> i;
	}
	datastream >> mass;
	datastream >> brownian_mean;
	cuboid_generator(particles, origin, scale, distance, velocity, mass, brownian_mean);
	spdlog::trace("Read cuboid: {}, {}, {}, {}, {}, {}", origin, scale, distance, velocity, mass, brownian_mean);
}
}  // namespace detail

namespace FileReader {
void readFile(ParticleContainer& particles, std::string_view filename) {
	int num_particles = 0;

	// Comes from argv, will be null terminated.
	// NOLINTNEXTLINE(*suspicious-stringview-data-usage)
	std::ifstream input_file(filename.data());
	std::string tmp_string;

	if (!input_file.is_open()) {
		spdlog::error("Error: could not open file: {}", filename);
		exit(-1);  // NOLINT(*mt-unsafe)
	}
	while (tmp_string.empty() or tmp_string[0] == '#') {
		getline(input_file, tmp_string);
		spdlog::trace("Read line: {}", tmp_string);
	}

	std::istringstream numstream(tmp_string);
	numstream >> num_particles;
	spdlog::trace("Reading {}.", num_particles);
	getline(input_file, tmp_string);
	spdlog::trace("Read line: {}", tmp_string);

	if (num_particles <= 0) {
		spdlog::error("Error reading file: non-positive particle count: {}", num_particles);
		exit(-1);  // NOLINT(*mt-unsafe)
	}
	particles.reserve(static_cast<std::size_t>(num_particles));

	for (int i = 0; i < num_particles; i++) {
		std::istringstream datastream(tmp_string);
		std::string type;
		datastream >> type;
		if (type == "particle") {
			detail::parse_particle(particles, datastream);
		} else if (type == "cuboid") {
			detail::parse_cuboid(particles, datastream);
		} else {
			spdlog::error("Unrecognized entity: {}", type);
			exit(-1);  // NOLINT(*mt-unsafe)
		}

		getline(input_file, tmp_string);

		spdlog::trace("Read line: {}", tmp_string);
	}
}
}  // namespace FileReader
