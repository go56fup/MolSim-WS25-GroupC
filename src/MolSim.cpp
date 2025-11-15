#include <cstdlib>
#include <exception>
#include <filesystem>
#include <sstream>
#include <vector>

#include <argparse/argparse.hpp>
#include <spdlog/common.h>
#include <spdlog/spdlog.h>

#include "FileReader.h"
#include "MolSim.h"
#include "ParticleContainer.h"

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
	// NOLINTBEGIN(*magic-numbers)

	program.add_argument("-d", "--delta-t")
		.help("the tick length used for the simulation")
		.default_value(0.014)
		.scan<'g', double>();
	program.add_argument("-e", "--end-time")
		.help("point in time to simulate until")
		.default_value(1000.0)
		.scan<'g', double>();

	// NOLINTEND(*magic-numbers)

	program.add_argument("-o", "--output").help("output directory for resulting files").default_value(".");
	program.add_argument("--log-level")
		.help("verbosity of log output")
		.choices("off", "critcal", "err", "warn", "info", "debug", "trace")
		.default_value("info");

	// TODO(tuna): look into making positionals meshable, so -f val file1 file2 -f2 val should work
	program.add_argument("files")
		.help("input file(s) to read initial condition of particles from")
		.nargs(argparse::nargs_pattern::at_least_one);

	try {
		program.parse_args(argc, argv);
	} catch (const std::exception& err) {
		spdlog::error(err.what());
		std::ostringstream os;
		os << program;
		spdlog::error(os.str());
		std::exit(1);  // NOLINT(*mt-unsafe)
	}
	spdlog::set_level(spdlog::level::from_str(program.get("--log-level")));

	spdlog::info("Hello from MolSim for PSE!");
	const std::string output_name = program.get("--output");
	if (!std::filesystem::exists(output_name)) {
		std::filesystem::create_directories(output_name);
	}

	ParticleContainer particles;
	const auto files = program.get<std::vector<std::string>>("files");
	for (const auto& file : files) {
		FileReader::readFile(particles, file);
	}
	run_simulation(
		particles, {.delta_t = program.get<double>("--delta-t"), .end_time = program.get<double>("--end-time")},
		output_name
	);

	spdlog::info("output written. Terminating...");
	return 0;
}

#ifdef BENCHMARK_PRESET
#include "Presets.h"
#endif

// NOLINTNEXTLINE(bugprone-exception-escape)
int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
#ifdef BENCHMARK_PRESET
	presets::BENCHMARK_PRESET<{.create_output = false}>("unused");
#else
	return usual_main(argc, argv);
#endif
}
