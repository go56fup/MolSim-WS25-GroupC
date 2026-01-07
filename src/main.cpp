#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string_view>

#include <argparse/argparse.hpp>
#include <fmt/base.h>
#include <spdlog/common.h>
#include <spdlog/spdlog.h>

#include "physics/forces.hpp"
#include "physics/particle.hpp"
#include "simulation/config/entities.hpp"
#include "simulation/config/parse.hpp"
#include "simulation/molsim.hpp"

struct use_formatting_t {};

inline constexpr use_formatting_t use_formatting{};

// TODO(tuna): look into how spdlog implements its overloads for fmt and remove the tag type if
// possible
// TODO(tuna): grep for exit and replace with terminating_log -- look up docs for die
template <typename... Args>
[[noreturn]] void terminating_log(fmt::format_string<Args...> fmt_string, Args&&... args) {
	spdlog::critical(fmt_string, std::forward<Args>(args)...);
	std::exit(1);  // NOLINT(*mt-unsafe)
}

[[noreturn]] void terminating_log(std::string_view msg) {
	spdlog::critical(msg);
	std::exit(1);  // NOLINT(*mt-unsafe)
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

	program.add_argument("file").help(
		"input file to read simulation parameters and initial condition of bodies from"
	);

	program.add_argument("-c", "--checkpoint")
		.help("list of particle states to add to simulation, usually generated via a checkpoint");

	try {
		program.parse_args(argc, argv);
	} catch (const std::exception& err) {
		spdlog::error(err.what());
		std::ostringstream os;
		os << program;
		terminating_log(os.str());
	}
	spdlog::set_level(spdlog::level::from_str(program.get("--log-level")));

	spdlog::info("Hello from MolSim for PSE!");
	const std::string output_name = program.get("--output");
	if (!std::filesystem::exists(output_name)) {
		std::filesystem::create_directories(output_name);
	}

	const auto file = program.get<std::string>("file");
	const std::ifstream t(file);
	std::stringstream buffer;
	buffer << t.rdbuf();
	unprocessed_config parse_result = config::parse(buffer.view());
	particle_container container(parse_result.config.domain, parse_result.config.cutoff_radius);
	config::populate_simulation(container, parse_result.config, parse_result.bodies);
	if (auto checkpoint = program.present("--checkpoint")) {
		const std::ifstream t(*checkpoint);
		std::stringstream buffer;
		buffer << t.rdbuf();
		config::populate_simulation(
			container, parse_result.config,
			daw::json::from_json<daw::json::json_value>(buffer.view())
		);
	}

	// TODO(tuna): when we add gravitation as an alternative force, change this to not be hardcoded
	// TODO(tuna): when we move to the lut approach for the epsilons, see if just comparing the type
	// is faster than branchless compute

	run_simulation(container, parse_result.config, lennard_jones_force, output_name);

	spdlog::info("output written. Terminating...");
	return 0;
}

#ifdef BENCHMARK_PRESET
#include "Presets.h"
#endif

// NOLINTNEXTLINE(bugprone-exception-escape)
int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
	// TODO(tuna): make this work
#ifdef BENCHMARK_PRESET
	presets::BENCHMARK_PRESET<{.create_output = false}>("unused");
#else
	return usual_main(argc, argv);
#endif
}
