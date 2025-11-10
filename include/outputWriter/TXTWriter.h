//
// Created by Gabriel Ribeiro Fernandes on 04.11.25.
//

#include <fstream>
#include <iomanip>
#include <sstream>

#include "Particle.h"

namespace outputWriter::TXTWriter {
inline void plotParticles(std::span<const Particle> particles, std::string_view filename, int iteration) {
	std::ofstream file;
	std::stringstream strstr;
	strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".txt";

	file.open(strstr.str().c_str());

	for (const auto& p : particles) {
		file << p.toString() << '\n';
	}
	
	file.close();
}
}  // namespace outputWriter::TXTWriter
