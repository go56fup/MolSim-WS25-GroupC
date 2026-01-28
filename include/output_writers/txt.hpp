//
// Created by Gabriel Ribeiro Fernandes on 04.11.25.
//

#include <fstream>
#include <iomanip>
#include <sstream>

#include "physics/particle.hpp"

namespace output_writer::txt {
inline void
plotParticles(std::span<const Particle> particles, std::string_view filename, unsigned iteration) {
	std::ofstream file;
	std::stringstream strstr;
	strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".txt";

	file.open(strstr.str().c_str());

	for (const auto& p : particles) {
		file << p << '\n';
	}

	file.close();
}
}  // namespace output_writer::txt
