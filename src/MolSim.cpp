#if __cpp_lib_to_chars == 201611L
#include <charconv>
#endif
#include <cstdlib>
#include <exception>
#include <iostream>
#include <span>
#include <string>
#include <string_view>

#include "FileReader.h"
#include "MolSim.h"
#include "ParticleContainer.h"

/**
 * @brief Prints an error message and terminates the program.
 *
 * Used when command-line arguments are invalid or missing.
 *
 * @note This function never returns.
 */
[[noreturn]] void die() {
	std::cout << "Erroneous programme call!\n./molsym filename\n";
	std::exit(1);  // NOLINT(*mt-unsafe)
}

/**
 * @brief Parses an argument from argv into a double value.
 *
 * Falls back to `std::stod` if `std::from_chars` is unavailable.
 * Exits the program if parsing fails.
 *
 * @param str View of the passed command-line argument.
 * @return The parsed double value.
 */
constexpr double get_double_from_argv(std::string_view str) {
	auto invalid_input = [] {
		std::cout << "Given decimal number is invalid.\n";
		die();
	};
	double result = 0;
#if __cpp_lib_to_chars == 201611L
	const auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.length(), result);
	if (*ptr != '\0' || ec == std::errc::invalid_argument || ec == std::errc::result_out_of_range) {
		invalid_input();
	}
#else
	std::size_t pos = 0;
	try {
		// Because the string_view points to argv entries, which are guaranteed to be null-terminated,
		// this is fine.
		// NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
		result = std::stod(str.data(), &pos);
	} catch (std::exception&) {
		invalid_input();
	}
	if (pos != str.length()) invalid_input();
#endif
	return result;
}

/**
 * @brief Returns a simulation parameter from either the default or a command-line argument.
 *
 * @param use_default Whether to use the default value.
 * @param input Pointer to a command-line argument to parse (ignored if @p use_default is `true`).
 * @param default_ The default value to use.
 * @return The parsed or default parameter.
 */
constexpr double get_sim_parameter(bool use_default, const char* input, double default_) {
	if (use_default) return default_;
	return get_double_from_argv(input);
}

/**
 * @brief Entry point of our simulation program.
 *
 * The program performs a basic molecular dynamics simulation:
 *  - Reads initial particle data from a file.
 *  - Calculates particle states over time.
 *  - Periodically writes simulation snapshots to VTK files.
 *
 * Command-line usage:
 * ```
 * ./MolSim input_file [delta_t end_time]
 * ```
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit status code (0 on success).
 */
int main(int argc, char* argv[]) {
	static constexpr int argc_with_everything_defaulted = 2;
	static constexpr int argc_with_everything_explicit = 5;
	std::cout << "Hello from MolSim for PSE!\n";
	if (argc < argc_with_everything_defaulted || argc > argc_with_everything_explicit) die();
	const bool defaulted = argc == 2;
	// NOLINTBEGIN(*pointer-arithmetic)
	const double delta_t = get_sim_parameter(defaulted, argv[2], 0.014);
	const double end_time = get_sim_parameter(defaulted, argv[3], 1000);

	ParticleContainer particles;
	FileReader::readFile(particles, argv[1]);
	// NOLINTEND(*pointer-arithmetic)

	run_simulation(particles, delta_t, end_time);
	std::cout << "output written. Terminating...\n";
	return 0;
}
