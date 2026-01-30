#pragma once

#include <fstream>
#include <mutex>
#include <string>
#include <string_view>

#include "grid/particle_container/fwd.hpp"
#include "simulation/entities.hpp"
#include "simulation/statistics.hpp"

inline void write_statistics_output(
	particle_container& container, const sim_configuration& config, std::string_view output_prefix,
	sim_iteration_t iteration
) {
	static const std::string densities_file_path = std::string(output_prefix) + "_densities.csv";
	static const std::string variance_file_path = std::string(output_prefix) + "_variance.csv";
	static std::once_flag csv_headers_created;

	std::call_once(csv_headers_created, [&] {
		std::ofstream densities_begin(densities_file_path);
		std::ofstream variance_begin(variance_file_path);

		variance_begin << "iteration,variance\n";
		densities_begin << "iteration,interval,density\n";
	});
	std::ofstream densities_file(densities_file_path, std::ios::app);
	std::ofstream variance_file(variance_file_path, std::ios::app);

	const auto densities = radial_distribution_function(container, config);

	for (std::size_t i = 0; i < densities.size(); ++i) {
		densities_file << iteration << ',' << i << ',' << densities[i] << '\n';
	}
	const auto variance = diffusion(container, config);
	variance_file << iteration << ',' << variance << '\n';
};
