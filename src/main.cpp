#include <chrono>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ratio>
#include <sstream>
#include <string_view>
#include <vector>

#include <argparse/argparse.hpp>
#include <daw/json/daw_json_exception.h>
#include <fmt/base.h>
#include <spdlog/common.h>
#include <spdlog/spdlog.h>

#include "grid/particle_container/fwd.hpp"
#include "output_writers/checkpoint.hpp"
#include "simulation/config/parse.hpp"
#include "simulation/entities.hpp"
#include "simulation/molsim.hpp"

template <typename... Args>
[[noreturn]] void terminating_log(fmt::format_string<Args...> fmt_string, Args&&... args) {
	SPDLOG_CRITICAL(fmt_string, std::forward<Args>(args)...);
	std::exit(1);  // NOLINT(*mt-unsafe)
}

[[noreturn]] void terminating_log(std::string_view msg) {
	SPDLOG_CRITICAL(msg);
	std::exit(1);  // NOLINT(*mt-unsafe)
}

[[noreturn]] void unhandled_json(const daw::json::json_exception& ex) {
	terminating_log("Uncaught JSON exception: {}", ex.what());
}

/**
 * @brief Entry point of our simulation program.
 *
 * The program performs a molecular dynamics simulation:
 *  - Reads initial particle data from a file.
 *  - Calculates particle states over time.
 *  - Periodically writes simulation snapshots to output files.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit status code (0 on success).
 */
// NOLINTNEXTLINE(*avoid-c-arrays)
int usual_main(int argc, char* argv[]) {
	argparse::ArgumentParser program("MolSim");
	// TODO(tuna): mark all NOLINTS with explanations
	program.add_argument("-o", "--output")
		.help("output directory for resulting files")
		.default_value(".");
	program.add_argument("--log-level")
		.help("verbosity of log output")
		.choices("off", "critcal", "err", "warn", "info", "debug", "trace")
		.default_value("info");

	program.add_argument("-c", "--config")
		.help("input file to read global simulation parameters from");

	program.add_argument("-b", "--bodies")
		.help("input file to read initial configuration of particles and bodies from")
		.append();

	try {
		program.parse_args(argc, argv);
	} catch (const std::exception& err) {
		spdlog::error(err.what());
		std::ostringstream os;
		os << program;
		terminating_log(os.str());
	}
	spdlog::set_level(spdlog::level::from_str(program.get("--log-level")));

	SPDLOG_INFO("Hello from MolSim for PSE!");
	const std::string output_name = program.get("--output");
	if (!std::filesystem::exists(output_name)) {
		std::filesystem::create_directories(output_name);
	}

	const auto config_filepath = program.get<std::string>("--config");
	const std::ifstream t(config_filepath);
	std::stringstream buffer;
	buffer << t.rdbuf();
	const std::string_view config_file_contents = buffer.view();
	sim_configuration config;

	try {
		config = config::parse_config(config_file_contents);
	} catch (const daw::json::json_exception& ex) {
		unhandled_json(ex);
	}

	const auto body_file_paths = program.get<std::vector<std::string>>("--bodies");

	particle_container container(config.domain, config.cutoff_radius);

	for (const auto& body_file : body_file_paths) {
		const std::ifstream t(body_file);
		std::stringstream body_buffer;
		body_buffer << t.rdbuf();
		const std::string_view body_file_contents = body_buffer.view();
		try {
			auto tail = config::parse_bodies(body_file_contents);
			config::populate_simulation(container, config, tail);
		} catch (const daw::json::json_exception& ex) {
			unhandled_json(ex);
		}
	}
	const std::string output_prefix = std::string(output_name) + "/" + config.base_name.c_str();

	const auto start = std::chrono::steady_clock::now();
	const auto iteration = run_simulation(container, config, output_prefix);
	const auto finish = std::chrono::steady_clock::now();
	const std::chrono::duration<double, std::milli> elapsed = finish - start;

	SPDLOG_INFO(
		"Runtime: {} ms, iterations: {}, tick length: {} ms", elapsed.count(), iteration,
		elapsed.count() / iteration
	);

	if (config.create_checkpoint) {
		// TODO(tuna): specify both in terms of plot_particles, where iteration is used in the
		// filename Or just remove plot_particles
		write_state_to_file(dump_state(container), output_prefix + "_checkpoint.json");
	}

	SPDLOG_INFO("output written. Terminating...");

	return 0;
}

#ifdef BENCHMARK_PRESET
#include "utility/presets.hpp"
#endif

// NOLINTNEXTLINE(bugprone-exception-escape)
int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
	// TODO(tuna): make this work
#ifdef BENCHMARK_PRESET
	presets::BENCHMARK_PRESET();
#else
	omp_set_num_threads(2);
	return usual_main(argc, argv);
#endif
}
