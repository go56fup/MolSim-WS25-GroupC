#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string_view>

#include <argparse/argparse.hpp>
#include <fmt/base.h>
#include <spdlog/common.h>
#include <spdlog/spdlog.h>

#include "FileReader.h"
#include "ForceCalculators.h"
#include "MolSim.h"
#include "Particle.h"

struct use_formatting_t {};

inline constexpr use_formatting_t use_formatting{};

// TODO(tuna): look into how spdlog implements its overloads for fmt and remove the tag type if possible
// TODO(tuna): grep for exit and replace with terminating_log -- look up docs for die
template <typename... Args>
[[noreturn]] void terminating_log(use_formatting_t, fmt::format_string<Args...> fmt_string, Args&&... args) {
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
	// Specifying defaults here.
	// NOLINTBEGIN(*magic-numbers)

	// NOLINTEND(*magic-numbers)

	program.add_argument("-o", "--output").help("output directory for resulting files").default_value(".");
	program.add_argument("--log-level")
		.help("verbosity of log output")
		.choices("off", "critcal", "err", "warn", "info", "debug", "trace")
		.default_value("info");

	// TODO(tuna): look into making positionals meshable, so -f val file1 file2 -f2 val should work
	// TODO(tuna): make exactly one
	program.add_argument("file").help("input file to read simulation parameters and initial condition of bodies from");

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
	auto [cfg, particles] = FileReader::parse(buffer.view());

	// TODO(tuna): when we add gravitation as an alternative force, change this to not be hardcoded
	auto calc = [sigma = cfg.sigma, eps = cfg.epsilon](const Particle& p1, const Particle& p2) noexcept {
		return lennard_jones_force(p1, p2, sigma, eps);
	};

	run_simulation(particles, cfg, calc, output_name);

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
